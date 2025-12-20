#!/usr/bin/env python3
"""
Unit tests for clone frequency computation.

Tests the paper's algorithm:
    clone_freq = mean(VAF of markers) - sum(child clone frequencies)

Run with: python -m pytest test_clone_frequency.py -v
Or directly: python test_clone_frequency.py
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from clone_frequency import (
    compute_clone_frequencies_from_ddpcr,
    get_post_order_traversal,
    get_children,
    find_root_node,
    get_all_nodes
)


def create_mock_tree_distribution():
    """
    Create a mock tree distribution for testing.

    Tree structure:
           0 (root, no mutations)
          / \
         1   2
        /
       3

    Node mutations:
    - Node 0: [] (root, no mutations)
    - Node 1: [GENE_A, GENE_B]
    - Node 2: [GENE_C]
    - Node 3: [GENE_D]
    """
    tree_structure = {
        0: [1, 2],  # Root has children 1 and 2
        1: [3],     # Node 1 has child 3
        2: [],      # Node 2 is a leaf
        3: []       # Node 3 is a leaf
    }

    node_dict_name = {
        0: [],                    # Root has no mutations
        1: ['GENE_A', 'GENE_B'],  # Node 1 has 2 mutations
        2: ['GENE_C'],            # Node 2 has 1 mutation
        3: ['GENE_D']             # Node 3 has 1 mutation
    }

    # Single tree with 100% frequency
    return {
        'freq': np.array([100.0]),
        'tree_structure': [tree_structure],
        'node_dict_name': [node_dict_name]
    }


def create_mock_ddpcr_data():
    """
    Create mock ddPCR data.

    VAFs:
    - GENE_A: 0.5 (50%)
    - GENE_B: 0.4 (40%)
    - GENE_C: 0.3 (30%)
    - GENE_D: 0.2 (20%)
    """
    return pd.DataFrame({
        'MutDOR': [500, 400, 300, 200],
        'DOR': [1000, 1000, 1000, 1000]
    }, index=['GENE_A', 'GENE_B', 'GENE_C', 'GENE_D'])


def test_get_post_order_traversal():
    """Test post-order traversal returns leaves first."""
    tree_dist = create_mock_tree_distribution()
    tree_structure = tree_dist['tree_structure'][0]
    node_dict_name = tree_dist['node_dict_name'][0]

    post_order = get_post_order_traversal(tree_structure, node_dict_name)

    # Leaves (2, 3) should come before their parents (1, 0)
    assert post_order.index(3) < post_order.index(1), "Node 3 should come before Node 1"
    assert post_order.index(2) < post_order.index(0), "Node 2 should come before Node 0"
    assert post_order.index(1) < post_order.index(0), "Node 1 should come before Node 0"

    print(f"Post-order traversal: {post_order}")
    print("PASS: Post-order traversal correct")


def test_clone_frequency_computation():
    """
    Test clone frequency computation using paper's algorithm.

    Expected calculations:
    - Node 3 (leaf): mean([0.2]) = 0.2, no children -> freq = 0.2
    - Node 1: mean([0.5, 0.4]) = 0.45, child sum = 0.2 -> freq = 0.45 - 0.2 = 0.25
    - Node 2 (leaf): mean([0.3]) = 0.3, no children -> freq = 0.3
    - Node 0: mean([]) = 0.0, child sum = 0.25 + 0.3 = 0.55 -> freq = 0.0 - 0.55 = -0.55
    """
    tree_dist = create_mock_tree_distribution()
    ddpcr_data = create_mock_ddpcr_data()

    clone_freqs = compute_clone_frequencies_from_ddpcr(
        tree_dist, ddpcr_data, "test_timepoint"
    )

    print(f"\nComputed clone frequencies: {clone_freqs}")

    # Check Node 3 (leaf)
    expected_3 = 0.2
    assert abs(clone_freqs[3] - expected_3) < 0.001, \
        f"Node 3: expected {expected_3}, got {clone_freqs[3]}"

    # Check Node 1
    expected_1 = 0.45 - 0.2  # mean([0.5, 0.4]) - child_freq[3]
    assert abs(clone_freqs[1] - expected_1) < 0.001, \
        f"Node 1: expected {expected_1}, got {clone_freqs[1]}"

    # Check Node 2 (leaf)
    expected_2 = 0.3
    assert abs(clone_freqs[2] - expected_2) < 0.001, \
        f"Node 2: expected {expected_2}, got {clone_freqs[2]}"

    # Check Node 0 (root with no mutations)
    # raw_vaf = 0.0 (no mutations), child_sum = 0.25 + 0.3 = 0.55
    expected_0 = 0.0 - (expected_1 + expected_2)
    assert abs(clone_freqs[0] - expected_0) < 0.001, \
        f"Node 0: expected {expected_0}, got {clone_freqs[0]}"

    print("PASS: Clone frequency computation correct")
    print(f"  Node 0 (root): {clone_freqs[0]:.4f} (negative expected - no mutations)")
    print(f"  Node 1: {clone_freqs[1]:.4f}")
    print(f"  Node 2: {clone_freqs[2]:.4f}")
    print(f"  Node 3: {clone_freqs[3]:.4f}")


def test_negative_frequencies_allowed():
    """
    Test that negative frequencies are allowed (per paper).

    Paper says negative values indicate measurement imprecision when
    children have higher VAF than parent.
    """
    tree_dist = create_mock_tree_distribution()
    ddpcr_data = create_mock_ddpcr_data()

    clone_freqs = compute_clone_frequencies_from_ddpcr(
        tree_dist, ddpcr_data, "test_timepoint"
    )

    # Node 0 should be negative (no mutations, but children have mutations)
    assert clone_freqs[0] < 0, "Node 0 should have negative frequency"
    print(f"\nNegative frequency for root: {clone_freqs[0]:.4f}")
    print("PASS: Negative frequencies correctly allowed")


def test_find_root_node():
    """Test root node detection."""
    tree_structure = {
        0: [1, 2],
        1: [3],
        2: [],
        3: []
    }

    root = find_root_node(tree_structure)
    assert root == 0, f"Expected root 0, got {root}"
    print("PASS: Root node detection correct")


def test_get_children():
    """Test child retrieval."""
    tree_structure = {
        0: [1, 2],
        1: [3],
        2: [],
        3: []
    }

    assert get_children(tree_structure, 0) == [1, 2]
    assert get_children(tree_structure, 1) == [3]
    assert get_children(tree_structure, 2) == []
    assert get_children(tree_structure, 3) == []
    assert get_children(tree_structure, 99) == []  # Non-existent node

    print("PASS: Child retrieval correct")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Testing clone frequency computation (paper's algorithm)")
    print("=" * 60)

    test_find_root_node()
    test_get_children()
    test_get_post_order_traversal()
    test_clone_frequency_computation()
    test_negative_frequencies_allowed()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
