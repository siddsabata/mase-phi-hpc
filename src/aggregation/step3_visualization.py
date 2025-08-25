"""
Step 3 Visualization Module for TracerX Marker Selection Pipeline

Main visualization module that combines tree rendering, analysis, and plotting
functionality for the aggregation stage of the pipeline.

Purpose:
- Combined tree and frequency visualizations
- Multi-sample analysis and plotting
- Tree distribution analysis with best tree identification
- Image combination and output management
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
import shutil
import os

# Import from our modular components
from tree_operations import collapse_nodes, ModifyTree
from tree_rendering import (
    render_tumor_tree, add_prefix_tree, df2dict, W2node_dict,
    generate_tree, generate_cp, root_searching, E2tree, validate_sample_consistency
)


def combine_existing_tree_and_frequency_plots(tree_png_path, freq_png_path, output_path, 
                                           patient_num, tree_idx, freq, type_name, actual_num_samples):
    """
    Combine existing tree and frequency PNG files into a side-by-side visualization.
    
    Creates publication-quality combined visualizations showing both phylogenetic
    tree structure and clonal frequency data in a single figure.
    
    Args:
        tree_png_path (str): Path to the Graphviz tree PNG file
        freq_png_path (str): Path to the matplotlib frequency PNG file
        output_path (str): Path for the combined output file
        patient_num (str): Patient identifier
        tree_idx (int): Tree index for title
        freq (float): Tree frequency for title
        type_name (str): Analysis type
        actual_num_samples (int): Number of samples
        
    Returns:
        str or None: Path to combined image if successful, None if failed
    """
    # Check if both input files exist
    if not Path(tree_png_path).exists():
        print(f"Warning: Tree PNG not found: {tree_png_path}")
        return None
    
    if not Path(freq_png_path).exists():
        print(f"Warning: Frequency PNG not found: {freq_png_path}")
        return None
    
    try:
        # Load the images
        tree_img = mpimg.imread(tree_png_path)
        freq_img = mpimg.imread(freq_png_path)
        
        # Create figure with side-by-side subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Display tree image on the left
        ax1.imshow(tree_img)
        ax1.axis('off')
        ax1.set_title('Phylogenetic Tree Structure', fontsize=14, fontweight='bold', pad=20)
        
        # Display frequency image on the right
        ax2.imshow(freq_img)
        ax2.axis('off')
        ax2.set_title('Clonal Frequencies', fontsize=14, fontweight='bold', pad=20)
        
        # Add main title for the entire figure
        fig.suptitle(f'{patient_num} - Tree {tree_idx} - Combined Analysis\n'
                    f'Bootstrap frequency: {freq} | Samples: {actual_num_samples} | Analysis: {type_name}', 
                    fontsize=16, fontweight='bold', y=0.95)
        
        # Adjust layout
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        
        # Save the combined visualization
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"Saved combined visualization: {output_path}")
        return output_path
        
    except Exception as e:
        print(f"Error combining images: {e}")
        return None


def analyze_tree_distribution(tree_distribution, directory, patient_num, type, fig=False, 
                            sample_prefix='Region', custom_sample_names=None):
    """
    Analyze tree distribution with multi-sample support and create combined tree-frequency visualizations.
    
    This function generates all tree visualizations in a subdirectory and identifies the best tree
    (highest frequency) to copy to the main results directory.
    
    Args:
        tree_distribution (dict): Tree distribution data from aggregation
        directory (Path): Output directory for files
        patient_num (str): Patient identifier
        type (str): Analysis type (e.g., 'initial')
        fig (bool): Whether to generate figures
        sample_prefix (str): Prefix for sample names (default: 'Region')
        custom_sample_names (list, optional): Custom sample names
    """
    if not fig:
        return
        
    # Create subdirectory for all tree visualizations
    all_trees_dir = directory / f'all_trees_{type}'
    all_trees_dir.mkdir(exist_ok=True, parents=True)
    
    # Track the best tree (highest frequency)
    best_tree_idx = None
    best_frequency = 0
    best_tree_combined_path = None
    
    for idx in range(len(tree_distribution['freq'])):
        tree_structure = tree_distribution['tree_structure'][idx]
        cp_tree = tree_distribution['cp_tree'][idx]
        node_dict = tree_distribution['node_dict'][idx]
        node_dict_name = tree_distribution['node_dict_name'][idx]
        freq = tree_distribution['freq'][idx]
        clonal_freq = tree_distribution['clonal_freq'][idx]

        # Track the best tree (highest frequency)
        if freq > best_frequency:
            best_frequency = freq
            best_tree_idx = idx

        # Validate sample consistency and detect number of samples
        num_samples = validate_sample_consistency(clonal_freq)
        if num_samples <= 0:
            print(f"Error: Invalid or inconsistent sample data for tree {idx}. Skipping visualization.")
            continue
        
        print(f"Processing tree {idx} with {num_samples} samples (frequency: {freq})")

        # Multi-sample clonal prevalence processing
        prev_mat = []
        for node, freqs in clonal_freq.items():
            for freq_sample in freqs:
                # Dynamic multi-sample processing
                actual_num_samples = len(freq_sample)
                for sample_idx in range(actual_num_samples):
                    # Generate sample names
                    if custom_sample_names and sample_idx < len(custom_sample_names):
                        sample_name = custom_sample_names[sample_idx]
                    else:
                        sample_name = f'{sample_prefix}_{sample_idx}'
                    
                    prev_mat.append({
                        'fraction': freq_sample[sample_idx], 
                        'sample': sample_name, 
                        'clone': node
                    })
        
        df_prev = pd.DataFrame(prev_mat)
        
        # Step 1: Generate phylogenetic tree (Graphviz) - temporary for combining
        g = render_tumor_tree(tree_structure, node_dict_name)
        tree_filename = all_trees_dir / f'{patient_num}_tree_dist{idx}_{type}_temp'
        g.render(filename=tree_filename, cleanup=True)
        tree_png_path = f"{tree_filename}.png"

        # Step 2: Generate frequency plot (matplotlib) - temporary for combining
        plt.figure(figsize=(12, 8))  # Slightly taller for better readability in combined view
        
        # Use seaborn color palette for better multi-sample visualization
        actual_num_samples = len(df_prev['sample'].unique())
        colors = sns.color_palette("Set2", actual_num_samples)
        
        # Create the bar plot with improved aesthetics
        sns.barplot(data=df_prev, x='clone', y='fraction', hue='sample', palette=colors)
        
        # Improved title with sample count information
        plt.title(f'{patient_num}_tree_{idx}_freq{freq} ({actual_num_samples} samples)', 
                 fontsize=14, fontweight='bold')
        
        # Better legend positioning for multiple samples
        plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Improve axis labels
        plt.xlabel('Clone', fontweight='bold', fontsize=12)
        plt.ylabel('Clonal Frequency', fontweight='bold', fontsize=12)
        
        # Rotate x-axis labels if many clones
        if len(df_prev['clone'].unique()) > 8:
            plt.xticks(rotation=45)
        
        # Add grid for better readability
        plt.grid(True, alpha=0.3, axis='y')
        
        # Save frequency plot temporarily
        freq_filename = all_trees_dir / f'{patient_num}_freq_dist{idx}_{type}_temp.png'
        plt.tight_layout()
        plt.savefig(freq_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Step 3: Combine both plots side-by-side and save to subdirectory
        combined_filename = all_trees_dir / f'{patient_num}_combined_tree_freq_{idx}_{type}.png'
        combined_result = combine_existing_tree_and_frequency_plots(
            tree_png_path=tree_png_path,
            freq_png_path=freq_filename,
            output_path=combined_filename,
            patient_num=patient_num,
            tree_idx=idx,
            freq=freq,
            type_name=type,
            actual_num_samples=actual_num_samples
        )
        
        # Track the best tree's combined visualization path
        if idx == best_tree_idx:
            best_tree_combined_path = combined_filename
        
        # Clean up temporary files
        try:
            if os.path.exists(tree_png_path):
                os.remove(tree_png_path)
            if os.path.exists(freq_filename):
                os.remove(freq_filename)
        except Exception as e:
            print(f"Warning: Could not clean up temporary files: {e}")
        
        if combined_result:
            print(f"Successfully created combined visualization for tree {idx}")
        else:
            print(f"Failed to create combined visualization for tree {idx}")
        
        print(f"Saved combined visualization for tree {idx} with {actual_num_samples} samples")

    # Copy the best tree visualization to the main aggregation results directory
    if best_tree_combined_path and best_tree_combined_path.exists():
        best_tree_main_path = directory / f'{patient_num}_combined_best_tree_{type}.png'
        try:
            shutil.copy2(best_tree_combined_path, best_tree_main_path)
            print(f"\nBest tree (index {best_tree_idx}, frequency {best_frequency}) copied to main directory:")
            print(f"  {best_tree_main_path}")
            print(f"All tree visualizations available in: {all_trees_dir}")
        except Exception as e:
            print(f"Warning: Failed to copy best tree visualization to main directory: {e}")
    else:
        print(f"Warning: Best tree visualization not found for copying to main directory")


def plot_mut_profile_comparison(tree_distribution, aggregated_results_file, directory, patient):
    """
    Plot mutation profile comparison across trees.
    
    This is a placeholder function for future mutation profile analysis.
    Currently provides basic framework for comparison visualizations.
    
    Args:
        tree_distribution (dict): Tree distribution data
        aggregated_results_file (str): Path to aggregated results file
        directory (Path): Output directory
        patient (str): Patient identifier
    """
    # Placeholder implementation for mutation profile comparison
    print(f"Mutation profile comparison for patient {patient} - functionality to be implemented")
    
    # Future implementation would include:
    # - Cross-tree mutation frequency analysis
    # - Mutation co-occurrence patterns
    # - Clonal evolution trajectory visualization
    # - Statistical comparison of mutation profiles
    pass