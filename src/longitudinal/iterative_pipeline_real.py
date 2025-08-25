from optimize import *                    # Import optimization functions
from optimize_fraction import *           # Import fraction optimization functions
import pandas as pd                       # For data manipulation
from zipfile import ZipFile               # For handling zip files
import json                               # For JSON data handling
import gzip                               # For handling gzip compressed files
import pickle                             # For serializing and deserializing Python objects
from analyze import *                     # Import analysis functions
from optimize import *                    # Import optimization functions (duplicate)
from adjust_tree_distribution import *    # Import tree distribution adjustment functions
import matplotlib.pyplot as plt           # For plotting


# Set patient and analysis parameters
patient_name='CRUK0041'                   # Patient identifier
used_num_blood=0                          # Number of blood samples used in bootstrapping
used_num_tissue=5                         # Number of tissue samples used in bootstrapping
type='common'                             # Type of analysis (common mutations)
directory = Path(f'/home/xuecongf/liquid_biopsy/abbosh')  # Main data directory
liquid_biopsy_directory = Path("/home/xuecongf/liquid_biopsy/abbosh/DateSample")  # Liquid biopsy data directory

# Alternative paths (commented out)
#directory = Path(f'Abbosh_et_al/data/bootstrap')
#liquid_biopsy_directory = Path("Abbosh_et_al/DateSample")

# Set analysis parameters
n_markers = 2                             # Number of markers to select
algo = "struct"                           # Algorithm type (structural)
file = directory / f'bootstrap/{patient_name}/{used_num_blood}_{used_num_tissue}/{type}'  # Path to bootstrap results
bootstrap_num = 100                       # Number of bootstrap samples
method = 'phylowgs'                       # Phylogenetic method used
num_chain=5                               # Number of chains for PhyloWGS
directory_xlsx = directory / f'{patient_name}_subset.xlsx'  # Path to patient data

# Load patient tissue data
inter = pd.read_excel(directory_xlsx, sheet_name='tissue_no_germline_0', index_col=0)
filtered = inter.copy(deep=True)          # Create a deep copy of the input data
calls = filtered                          # Variable for mutation calls

# Create mapping between gene IDs and indices
gene2idx = {'s' + str(i): i for i in range(len(inter))}  # Map gene IDs to indices
print(gene2idx)                           # Print the mapping
gene_list = list(gene2idx.keys())         # List of gene IDs
gene_name_list = []                       # Will store human-readable gene names
gene_count = {}                           # Counter for duplicate gene names

# Process gene names, handling duplicates and missing values
for i in range(inter.shape[0]):
    gene = calls["Gene"].loc[i]           # Get gene name
    if gene in gene_name_list:            # Check if gene name already exists
        gene_count[gene] += 1             # Increment count for duplicate genes
        gene = gene + '_' + str(gene_count[gene])  # Append counter to duplicate gene names
    else:
        gene_count[gene] = 1              # Initialize counter for new gene
    if not isinstance(gene, str):         # Handle non-string gene names
        gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])  # Use location as name
    gene_name_list.append(gene)           # Add processed gene name to list

num_marker = len(gene_list)               # Count total number of markers

# Load ddPCR data (liquid biopsy data)
directory_ddpcr = Path("/home/xuecongf/liquid_biopsy/abbosh/DateSample")
file_ddpcr = directory_ddpcr / f"{patient_name}__DateSample.xlsx"
ddpcr_raw = pd.read_excel(file_ddpcr, sheet_name=None, index_col=0)

# Process dates from ddPCR data and sort chronologically
date_list = []
for date in ddpcr_raw.keys():
    date_list.append(pd.to_datetime(date, format="%Y-%m-%d"))
date_keys_sorted = [ts.strftime("%Y-%m-%d") for ts in sorted(date_list)]  # Sort dates chronologically

# Iterate through each time point (ordered by date)
for order_idx in range(0, len(date_keys_sorted)):
    # For the first timepoint (initial analysis)
    if order_idx == 0:
        # Load initial tree distribution summary
        tree_distribution_file_summary = file / f'{method}_bootstrap_summary.pkl'
        with open(tree_distribution_file_summary, 'rb') as f:
            tree_distribution_summary = pickle.load(f)
        
        # Extract tree structures, node information, and frequencies
        tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
            tree_distribution_summary['tree_structure'], tree_distribution_summary[
                'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
        
        # Load full tree distribution data
        tree_distribution_file = file / f'{method}_bootstrap_aggregation.pkl'
        with open(file / tree_distribution_file, 'rb') as f:
            tree_distribution = pickle.load(f)
        
        # Extract tree structures, node dictionary, clonal frequencies, and tree frequencies
        tree_list, node_list, clonal_freq_list, tree_freq_list = tree_distribution['tree_structure'], tree_distribution[
            'node_dict'], tree_distribution['vaf_frac'], tree_distribution['freq']
    
    # For subsequent timepoints (using updated tree distributions)
    else:
        # Load the updated tree distribution from previous timepoint
        tree_distribution_file_summary = file / f'{method}_bootstrap_summary_updated_{algo}_{n_markers}_{order_idx-1}_bayesian.pkl'
        with open(tree_distribution_file_summary, 'rb') as f:
            tree_distribution_summary = pickle.load(f)
        
        # Extract tree structures, node information, and frequencies
        tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
            tree_distribution_summary['tree_structure'], tree_distribution_summary[
                'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
        tree_list, node_list, tree_freq_list = tree_distribution_summary['tree_structure'], tree_distribution_summary[
            'node_dict'], tree_distribution_summary['freq']
        
        # Recalculate clonal frequencies by averaging across samples
        clonal_freq_list = []
        for idx in range(len(tree_distribution_summary['vaf_frac'])):
            clonal_freq_dict = tree_distribution_summary['vaf_frac'][idx]
            clonal_freq_dict_new = {}
            for node, freqs in clonal_freq_dict.items():
                clonal_freq_dict_new[node] = [list(np.array(freqs).mean(axis=0))]
            clonal_freq_list.append(clonal_freq_dict_new)

    # Clean up node list data: convert string keys to integers
    node_list_scrub = []
    for node_dict in node_list:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_scrub.append(temp)

    # Clean up clonal frequency data: convert string keys to integers
    clonal_freq_list_scrub = []
    for clonal_freq_dict in clonal_freq_list:
        temp = {}
        for key, values in clonal_freq_dict.items():
            temp.setdefault(int(key), values[0])
        clonal_freq_list_scrub.append(temp)

    # Load ddPCR data for current timepoint
    sample_df = pd.read_excel(liquid_biopsy_directory / f"{patient_name}__DateSample.xlsx", sheet_name=None, index_col=0)
    subset_markers = list(sample_df[list(sample_df.keys())[0]].index)  # Get markers from first sheet
    subset_list = list(inter[inter.Gene.isin(subset_markers)].index)  # Get indices of markers in original data
    subset_markers_s = list([f"s{i}" for i in subset_list])  # Format marker IDs
    gene2idx_sub = {subset_markers_s[i]: i for i in range(len(subset_markers_s))}  # Create mapping for subset

    # Parameters for marker selection
    read_depth=90000                     # Sequencing read depth
    lam1 = 0                             # Weight for tree fractions (not used)
    lam2 = 1                             # Weight for tree distributions (fully weighted)

    # Select optimal markers based on tree structure
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(
        gene_list, n_markers, tree_list, node_list, clonal_freq_list, 
        gene2idx, tree_freq_list, read_depth=read_depth, 
        lam1=lam1, lam2=lam2, focus_sample_idx=0)
    
    # Convert selected marker IDs to gene names
    selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
    ddpcr = []  # Will store ddPCR data for selected markers

    # Get current blood sample timepoint
    blood_sample_idx = date_keys_sorted[order_idx]
    ddpcr_raw_sample = ddpcr_raw[blood_sample_idx]
    
    # Extract ddPCR data for selected markers
    for gene in selected_markers2_genename:
        ddpcr.append({
            'gene': gene, 
            'mut': ddpcr_raw_sample["MutDOR"].loc[gene],  # Mutant reads
            'WT': ddpcr_raw_sample["DOR"].loc[gene],      # Total depth
            'liquid_biopsy_sample': blood_sample_idx      # Sample date
        })

    # Create DataFrame with ddPCR data
    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}  # Map indices to genes

    # Extract tree information from summary
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
    
    # Extract ddPCR counts for markers
    ddpcr_marker_counts = list(df_ddpcr_2["mut"])  # Mutant read counts
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])  # Total read depths

    # Update tree distributions using Bayesian approach based on new ddPCR data
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(
        tree_list_summary, node_name_list_summary,
        tree_freq_list_summary, read_depth_list,
        ddpcr_marker_counts, marker_idx2gene)
    
    # Create updated tree distribution summary
    updated_tree_distribution_summary = update_tree_distribution_bayesian(
        tree_distribution_summary, updated_tree_freq_list)

    # Save updated tree distribution for next iteration
    tree_distribution_file_summary_updated = file / f'{method}_bootstrap_summary_updated_{algo}_{n_markers}_{order_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)