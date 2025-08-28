"""
Backward Compatibility Module for TracerX Marker Selection Pipeline

This module provides backward compatibility by importing from the new modular structure.
All functions are now available through their respective specialized modules.

DEPRECATED: This file is maintained for backward compatibility only.
New code should import directly from the specific modules:
- tree_operations: For tree manipulation and node collapsing
- tree_rendering: For Graphviz visualization and tree utilities  
- visualization: For main visualization and analysis functions
"""

# Import all functions from the new modular structure for backward compatibility
from tree_operations import collapse_nodes, ModifyTree
from tree_rendering import (
    render_tumor_tree, add_prefix_tree, df2dict, W2node_dict,
    generate_tree, generate_cp, root_searching, E2tree, validate_sample_consistency
)
from visualization import (
    combine_existing_tree_and_frequency_plots, analyze_tree_distribution, 
    plot_mut_profile_comparison
)

# Re-export all functions to maintain compatibility
__all__ = [
    'collapse_nodes', 'ModifyTree',
    'render_tumor_tree', 'add_prefix_tree', 'df2dict', 'W2node_dict',
    'generate_tree', 'generate_cp', 'root_searching', 'E2tree', 'validate_sample_consistency',
    'combine_existing_tree_and_frequency_plots', 'analyze_tree_distribution', 
    'plot_mut_profile_comparison'
]