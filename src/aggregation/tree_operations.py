"""
Tree Operations Module for TracerX Marker Selection Pipeline

This module provides tree manipulation operations including node collapsing
and tree modification for phylogenetic analysis.

Purpose:
- Tree node collapsing based on frequency thresholds
- Tree structure modification and validation
- Supporting utilities for tree analysis
"""

import numpy as np


class ModifyTree:
    """
    Tree modification class for manipulating phylogenetic tree structures.
    
    Provides functionality to delete nodes and maintain tree relationships
    while preserving the edge matrix structure.
    """
    
    def __init__(self, E):
        """
        Initialize tree from edge matrix.
        
        Args:
            E (np.ndarray): Edge matrix representing parent-child relationships
        """
        self.cp_tree = {}  # child -> parent mapping
        self.tree = {}     # parent -> children mapping
        self.E = E
        self.N = len(E)
        
        # Build tree relationships from edge matrix
        for i in range(self.N - 1, -1, -1):
            for j in range(self.N - 1, -1, -1):
                if int(E[i, j]) == 1:
                    self.cp_tree[j] = i
                    if i not in self.tree.keys():
                        self.tree[i] = [j]
                    else:
                        self.tree[i].append(j)

    def delete_node(self, idx):
        """
        Delete a node from the tree and maintain connectivity.
        
        Args:
            idx (int): Index of node to delete
            
        Raises:
            Exception: If attempting to delete root with multiple children
        """
        if self.is_root(idx):
            if self.num_children(idx) > 1:
                raise Exception('Cannot delete root node with more than one child!')
            child = self.tree[idx][0]
            del self.cp_tree[child]
            del self.tree[idx]
        elif self.is_leaf(idx):
            parent = self.cp_tree[idx]
            del self.cp_tree[idx]
            if self.num_children(parent) == 1:
                del self.tree[parent]
            else:
                self.tree[parent].remove(idx)
        else:
            # Internal node: reconnect children to parent
            parent = self.cp_tree[idx]
            children = self.tree[idx]
            self.tree[parent].remove(idx)
            del self.cp_tree[idx]
            for child in children:
                self.cp_tree[child] = parent
                self.tree[parent].append(child)
                self.E[parent, child] = 1
            del self.tree[idx]

    def is_leaf(self, idx):
        """Check if node is a leaf (no children)."""
        return idx not in self.tree.keys()

    def is_root(self, idx):
        """Check if node is root (no parent)."""
        return idx in self.tree.keys() and idx not in self.cp_tree.keys()

    def num_children(self, idx):
        """Get number of children for a node."""
        if self.is_leaf(idx):
            return 0
        else:
            return len(self.tree[idx])


def collapse_nodes(U, C, E, A, W, threshold=0.0, only_leaf=False):
    """
    Collapse tree nodes based on frequency and structural criteria.
    
    This function removes nodes with zero branch lengths or low frequencies
    to simplify the tree structure while preserving essential relationships.
    
    Args:
        U (np.ndarray): Node frequency matrix
        C (np.ndarray): Clone matrix
        E (np.ndarray): Edge matrix
        A (np.ndarray): Adjacency matrix
        W (np.ndarray): Weight matrix
        threshold (float): Frequency threshold for node removal
        only_leaf (bool): If True, only collapse leaf nodes
        
    Returns:
        tuple: Updated matrices (U_new, C_new, E_new, A_new, W_new)
    """
    print("Loading collapse nodes")
    
    # Generate the tree structure
    tree = ModifyTree(E)
    
    if not only_leaf:
        # Step 1: Collapse branches with zero length
        branch_remove_idx = []
        for i in range(tree.N-1, -1, -1):
            for j in range(tree.N-1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0:
                    branch_remove_idx.append(j)
        
        # Remove zero-length branches
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            tree.delete_node(node)

        # Step 2: Collapse nodes with low frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        for i in range(tree.N-1, -1, -1):
            if i in branch_remove_idx:
                continue
            if tree.is_root(i):
                continue
            if np.mean(U[:, i]) <= threshold:
                if tree.num_children(i) == 1:
                    freq_remove_idx.append(i)
                elif tree.is_leaf(i):
                    freq_leaf_remove_idx.append(i)
        
        print(f"Frequency-based removal: {freq_remove_idx}")
        
        # Remove low-frequency internal nodes
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]
        
        # Remove low-frequency leaf nodes
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
    else:
        # Only leaf mode: more conservative collapsing
        branch_remove_idx = []
        print("Collapsing branch with length 0")
        
        # Only collapse zero-length branches that don't lead to leaf nodes
        for i in range(tree.N - 1, -1, -1):
            for j in range(tree.N - 1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and not tree.is_leaf(j):
                    branch_remove_idx.append(j)
        
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            tree.delete_node(node)

        # Collapse leaf nodes with zero frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        print('Collapsing leaf nodes with frequency 0')
        
        for i in range(tree.N - 1, -1, -1):
            if i in branch_remove_idx:
                continue
            if np.mean(U[:, i]) <= threshold and tree.is_leaf(i):
                freq_leaf_remove_idx.append(i)
        
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
        
        # Handle single-child internal nodes
        for i in range(tree.N - 1, -1, -1):
            if tree.num_children(i) == 1:
                freq_remove_idx.append(i)
        
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]

    # Combine all nodes to be removed
    remove_idx = branch_remove_idx + freq_remove_idx + freq_leaf_remove_idx
    print(f'Nodes {remove_idx} will be collapsed.')
    
    # Update all matrices by removing collapsed nodes
    U_new = np.delete(U, remove_idx, axis=1)
    C_new = np.delete(C, remove_idx, axis=0)
    A_new = np.delete(A, remove_idx, axis=0)
    A_new = np.delete(A_new, remove_idx, axis=1)
    E_new = np.delete(tree.E, remove_idx, axis=0)
    E_new = np.delete(E_new, remove_idx, axis=1)
    W_new = np.delete(W, remove_idx, axis=0)
    
    print(f"Collapse complete: U shape {U_new.shape}, C shape {C_new.shape}")
    return U_new, C_new, E_new, A_new, W_new