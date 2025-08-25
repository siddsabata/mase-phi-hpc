"""
Tree Rendering Module for TracerX Marker Selection Pipeline

This module provides tree visualization and rendering functionality
using Graphviz for phylogenetic tree display.

Purpose:
- Graphviz-based tree visualization
- Tree structure utilities and conversions
- Data format conversions for visualization
"""

import numpy as np
from graphviz import Digraph


def render_tumor_tree(tree_structure, node_dict):
    """
    Render phylogenetic tree using Graphviz with mutation annotations.
    
    Creates a directed graph visualization of the tumor evolution tree
    with nodes labeled by their associated mutations.
    
    Args:
        tree_structure (dict): Parent -> children tree structure
        node_dict (dict): Node -> mutation list mapping
        
    Returns:
        graphviz.Digraph: Rendered tree visualization
    """
    w = Digraph(format='png')
    edge_idx = 0
    root = root_searching(tree_structure)
    
    for p, c_list in tree_structure.items():
        # Handle root node specially (labeled as 'normal')
        p_list = ['normal'] if p == root else node_dict[p]
        p_final = tuple(set(p_list))
        ne = len(p_final)
        
        for c in c_list:
            edge_idx += 1
            
            # Get child node mutations
            if c not in node_dict:
                c_type = str(c)
            else:
                c_type = c
            c_final = tuple(set(node_dict[c_type]))
            ne_ = len(c_final)
            
            # Format node labels based on mutation count
            # Split long lists across multiple lines for readability
            if ne >= 10:
                p_str = str(p) + ' ' + ' '.join((str(e) for e in p_final[:ne//3])) + '\n' + \
                       ' '.join((str(e) for e in p_final[ne//3:ne//3*2])) + '\n' + \
                       ' '.join((str(e) for e in p_final[ne//3*2:]))
                if ne_ >= 10:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final[:ne_//3])) + '\n' + \
                           ' '.join((str(e) for e in c_final[ne_//3:ne_//3*2])) + '\n' + \
                           ' '.join((str(e) for e in c_final[ne_//3*2:]))
                else:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final))
            else:
                p_str = str(p) + ' ' + ' '.join((str(e) for e in p_final))
                if ne_ >= 10:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final[:ne_//3])) + '\n' + \
                           ' '.join((str(e) for e in c_final[ne_//3:ne_//3*2])) + '\n' + \
                           ' '.join((str(e) for e in c_final[ne_//3*2:]))
                else:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final))
            
            # Add edge with branch label
            w.edge(p_str, c_str, "b" + str(edge_idx))
    
    return w


def add_prefix_tree(mutation):
    """
    Add 'mut_' prefix to mutation identifiers.
    
    Args:
        mutation (dict): Node -> mutation list mapping
        
    Returns:
        dict: Node -> prefixed mutation list mapping
    """
    mutation_prefixed = {}
    for node, mut_list in mutation.items():
        mutation_prefixed[node] = ["mut_" + str(mut) for mut in mut_list]
    return mutation_prefixed


def df2dict(df):
    """
    Convert DataFrame to mutation index dictionary.
    
    Extracts gene names or genomic coordinates for mutation identification.
    
    Args:
        df (pd.DataFrame): DataFrame with gene and position information
        
    Returns:
        dict: Index -> mutation name mapping
    """
    idx2name = {}
    for i in range(len(df)):
        if isinstance(df.loc[i]["Gene"], str):
            name = df.loc[i]["Gene"]
        else:
            name = str(df.loc[i]["Chromosome"]) + '_' + str(df.loc[i]["Genomic Position"])
        idx2name.setdefault(i, name)
    return idx2name


def W2node_dict(W_node, idx2name=None):
    """
    Convert weight matrix to node-mutation dictionary.
    
    Creates mapping from tree nodes to their associated mutations
    based on the weight matrix structure.
    
    Args:
        W_node (np.ndarray): Weight matrix (nodes x mutations)
        idx2name (dict, optional): Index to name mapping
        
    Returns:
        dict: Node -> mutation list mapping
    """
    node_dict = {}
    N, m = W_node.shape
    
    for i in range(N):
        node_dict.setdefault(i, [])
        for j in range(m):
            if W_node[i, j] == 1:
                if idx2name is not None:
                    node_dict[i].append(idx2name[j])
                else:
                    node_dict[i].append(j)
    
    return node_dict


def generate_tree(cp_tree):
    """
    Convert child-parent mapping to parent-children tree structure.
    
    Args:
        cp_tree (dict): Child -> parent mapping
        
    Returns:
        dict: Parent -> children list mapping
    """
    tree = {}
    for child, parent in cp_tree.items():
        if parent in tree.keys():
            tree[parent].append(child)
        else:
            tree[parent] = [child]
    return tree


def generate_cp(tree):
    """
    Convert parent-children mapping to child-parent structure.
    
    Args:
        tree (dict): Parent -> children list mapping
        
    Returns:
        dict: Child -> parent mapping
    """
    return {c: p for p in tree.keys() for c in tree[p]}


def root_searching(tree):
    """
    Find root node of the tree structure.
    
    Traverses upward from any node until finding a node with no parent.
    Includes loop detection for safety.
    
    Args:
        tree (dict): Tree structure (parent -> children)
        
    Returns:
        int or None: Root node index, or None if loop detected
    """
    tree_cp = generate_cp(tree)
    
    if not tree_cp:  # Empty tree
        return None
        
    start_node = list(tree_cp.keys())[0]
    iter_count = 0
    
    while True:
        iter_count += 1
        start_node = tree_cp[start_node]
        
        if start_node not in tree_cp.keys():
            break
            
        if iter_count >= 100:
            print("The directed tree exists self-loop.")
            return None
    
    return start_node


def E2tree(E):
    """
    Convert edge matrix to tree structure.
    
    Args:
        E (np.ndarray): Edge matrix (parent x child)
        
    Returns:
        dict: Parent -> children list mapping
    """
    tree = {}
    N = E.shape[0]
    
    for i in range(N):
        for j in range(N):
            if E[i, j] == 1:
                tree.setdefault(i, [])
                tree[i].append(j)
    
    return tree


def validate_sample_consistency(clonal_freq_data):
    """
    Validate consistency of sample counts across tree nodes.
    
    Ensures all nodes have the same number of samples for proper
    multi-sample analysis.
    
    Args:
        clonal_freq_data (dict): Node -> frequency data mapping
        
    Returns:
        int: Number of samples if consistent, -1 if inconsistent
    """
    sample_counts = []
    
    for node, freqs in clonal_freq_data.items():
        for freq_sample in freqs:
            sample_counts.append(len(freq_sample))
    
    if len(set(sample_counts)) > 1:
        print(f"Warning: Inconsistent sample counts across nodes: {set(sample_counts)}")
        return -1
    
    return sample_counts[0] if sample_counts else 0