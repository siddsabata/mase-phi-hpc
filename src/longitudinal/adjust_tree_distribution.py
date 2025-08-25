import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations, permutations
from collections import deque

import pandas as pd
from scipy.stats import norm, chi2

from pathlib import Path
import pickle
from analyze import *
import math
from scipy.integrate import quad

def adjust_tree_distribution_struct(tree_list, node_dict_list, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha):
    accepted_tree_indices = []
    for tree_idx in range(len(tree_list)):
        tree_structure = tree_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        bool_reject = test_single_tree_struct(tree_structure, node_dict, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha)
        if not bool_reject:
            accepted_tree_indices.append(tree_idx)
    return accepted_tree_indices

def adjust_tree_distribution_frac(clonal_freq_list,node_dict_list, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha):
    accepted_tree_indices = []
    for tree_idx in range(len(clonal_freq_list)):
        clonal_freq_dict = clonal_freq_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        bool_reject = test_single_tree_frac(clonal_freq_dict, node_dict, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha)
        if not bool_reject:
            accepted_tree_indices.append(tree_idx)
    return accepted_tree_indices

def adjust_tree_distribution_struct_bayesian(tree_list, node_dict_list, tree_freq_list, read_depth, ddpcr_marker_counts, marker_idx2gene):
    """
    Update tree distribution using Bayesian inference with ALL available markers.
    
    MAJOR FIX: This version uses all marker pairs instead of randomly selecting one,
    providing deterministic and more robust tree updates. This resolves the 
    non-deterministic behavior where identical inputs produced different results.
    
    For n markers, this processes C(n,2) = n*(n-1)/2 marker pairs and combines
    their likelihood evidence to update tree frequencies.
    """
    updated_tree_freq_list = []
    n_markers = len(ddpcr_marker_counts)
    
    # Handle edge case: only 1 marker (can't form pairs)
    if n_markers < 2:
        print(f"Warning: Only {n_markers} marker(s) available. Need at least 2 for pairwise analysis.")
        return tree_freq_list  # Return unchanged frequencies
    
    # Generate ALL marker pairs instead of randomly selecting one
    all_marker_pairs = list(combinations(range(n_markers), 2))
    print(f"Using all {len(all_marker_pairs)} marker pairs: {all_marker_pairs}")
    
    for tree_idx in range(len(tree_list)):
        tree_structure = tree_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        tree_freq = tree_freq_list[tree_idx]
        
        # Calculate combined likelihood from ALL marker pairs
        combined_likelihood = 1.0  # Start with neutral likelihood
        
        for marker_pair in all_marker_pairs:
            pair_likelihood = calculate_single_pair_likelihood(
                tree_structure, node_dict, read_depth, ddpcr_marker_counts, 
                marker_idx2gene, marker_pair)
            combined_likelihood *= pair_likelihood
            
        # Update tree frequency with combined evidence
        updated_tree_freq = tree_freq * combined_likelihood
        updated_tree_freq_list.append(updated_tree_freq)
    
    # Normalize frequencies to sum to 100
    updated_tree_freq_list_np = np.array(updated_tree_freq_list)
    total_freq = np.sum(updated_tree_freq_list_np)
    if total_freq > 0:
        updated_tree_freq_list_np = (updated_tree_freq_list_np / total_freq) * 100
    else:
        # Fallback: uniform distribution if all likelihoods are zero
        updated_tree_freq_list_np = np.ones(len(tree_list)) * (100.0 / len(tree_list))
        print("Warning: All tree likelihoods were zero, using uniform distribution")
    
    return updated_tree_freq_list_np.tolist()

def calculate_single_pair_likelihood(tree_structure, node_dict, read_depth_list, ddpcr_marker_counts, marker_idx2gene, marker_pair, lower_bound=0, upper_bound=1):
    """
    Calculate the likelihood for a single pair of markers given a tree structure.
    
    Args:
        tree_structure: Phylogenetic tree structure
        node_dict: Mapping of nodes to mutations
        read_depth_list: List of read depths for each marker
        ddpcr_marker_counts: List of mutant counts for each marker
        marker_idx2gene: Mapping from marker index to gene name
        marker_pair: Tuple of (marker_idx1, marker_idx2)
        lower_bound: Lower bound for VAF integration
        upper_bound: Upper bound for VAF integration
        
    Returns:
        Conditional probability (likelihood) for this marker pair
    """
    marker_idx1, marker_idx2 = marker_pair
    a2d_matrix = ancestor2descendant(tree_structure)
    mut2node_dict = mut2node(node_dict)
    
    # Get nodes for each marker
    gene1 = marker_idx2gene[marker_idx1]
    gene2 = marker_idx2gene[marker_idx2]
    
    if gene1 not in mut2node_dict or gene2 not in mut2node_dict:
        print(f"Warning: Marker {gene1} or {gene2} not found in tree. Skipping this pair.")
        return 1.0  # Neutral likelihood
    
    node1, node2 = mut2node_dict[gene1], mut2node_dict[gene2]
    
    # Prepare data for likelihood calculation
    read_depth_sublist = [int(np.round(read_depth_list[marker_idx1])), int(np.round(read_depth_list[marker_idx2]))]
    read_counts_sublist = [int(np.round(ddpcr_marker_counts[marker_idx1])), int(np.round(ddpcr_marker_counts[marker_idx2]))]
    
    print(f"Marker pair ({gene1}, {gene2}): depths={read_depth_sublist}, counts={read_counts_sublist}")
    
    # Determine phylogenetic relationship
    if node1 == node2:
        relation = 'same'
    elif a2d_matrix[node1, node2] == 1:
        relation = 'ancestor'
    elif a2d_matrix[node1, node2] == 0 and a2d_matrix[node2, node1] == 1:
        relation = 'descendant'
    else:
        relation = 'null'
    
    print(f"Relationship: {relation}")
    
    # Calculate conditional probability based on relationship
    if relation == 'same':
        conditional_prob = outer_integral_single(read_depth_sublist[0], read_depth_sublist[1], 
                                               read_counts_sublist[0], read_counts_sublist[1], 
                                               lower_bound, upper_bound)
    elif relation == 'null':
        # For unrelated markers, use neutral likelihood
        conditional_prob = 1.0
        print(f"Unrelated markers, using neutral likelihood")
    else:
        conditional_prob = outer_integral(read_depth_sublist[0], read_depth_sublist[1], 
                                        read_counts_sublist[0], read_counts_sublist[1], 
                                        relation, lower_bound, upper_bound)
    
    print(f"Conditional probability: {conditional_prob}")
    return conditional_prob


def update_single_tree_fractions(tree_structure, node_dict, tree_freq, read_depth_list, ddpcr_marker_counts, marker_idx2gene, marker_idx_list, lower_bound=0, upper_bound=1):
    """
    Legacy function - kept for backward compatibility.
    Use calculate_single_pair_likelihood for new code.
    """
    marker_pair = tuple(marker_idx_list)
    conditional_prob = calculate_single_pair_likelihood(
        tree_structure, node_dict, read_depth_list, ddpcr_marker_counts, 
        marker_idx2gene, marker_pair, lower_bound, upper_bound)
    updated_tree_freq = tree_freq * conditional_prob
    return updated_tree_freq

def ncr(n, r): ##for python 3.7
    return math.factorial(n) // math.factorial(r) // math.factorial(n-r)
def integrand(f2, f1, d1, d2, r1, r2):
    try:
        # return math.comb(d1, r1) * (f1 ** r1) * ((1 - f1) ** (d1 - r1)) * math.comb(d2, r2) * (f2 ** r2) * ((1 - f2) ** (d2 - r2))
        return math.exp(math.log(math.comb(d1, r1)) + r1 * math.log(f1) + (d1 - r1) * math.log(1 - f1) +
               math.log(math.comb(d2, r2)) + r2 * math.log(f2) + (d2 - r2) * math.log(1 - f2))
    except AttributeError:
        return ncr(d1, r1) * (f1 ** r1) *((1 - f1)**(d1-r1)) * ncr(d2, r2) * (f2 ** r2) * ((1-f2) ** (d2-r2))

def recursive_integral(f1, d1, d2, r1, r2, relation,lower_bound, upper_bound):
    integral, error = quad(integrand, f2_start(f1, relation, lower_bound), f2_end(f1, relation, upper_bound), args=(f1, d1, d2, r1, r2))
    return integral

def outer_integral(d1, d2, r1, r2, relation, lower_bound, upper_bound):
    integral, error = quad(recursive_integral, f1_start(lower_bound), f1_end(upper_bound), args=(d1, d2, r1, r2,relation,lower_bound, upper_bound))
    return integral

def integrand_single(f, d1, d2, r1, r2):
    """
    Compute the integrand using log-space arithmetic to avoid numerical overflow.
    This calculates the binomial likelihood for two markers with VAF f.
    """
    # Clamp f to avoid log(0) issues
    f = max(1e-16, min(1 - 1e-16, f))
    
    try:
        # Use log-space arithmetic to prevent overflow
        # log(C(n,k)) = log(n!) - log(k!) - log((n-k)!)
        from scipy.special import gammaln
        
        # Log binomial coefficients using gammaln (more stable than math.lgamma)
        log_comb1 = gammaln(d1 + 1) - gammaln(r1 + 1) - gammaln(d1 - r1 + 1)
        log_comb2 = gammaln(d2 + 1) - gammaln(r2 + 1) - gammaln(d2 - r2 + 1)
        
        # Log likelihood terms
        log_likelihood = (log_comb1 + log_comb2 + 
                         r1 * np.log(f) + (d1 - r1) * np.log(1 - f) +
                         r2 * np.log(f) + (d2 - r2) * np.log(1 - f))
        
        # Convert back from log space, with overflow protection
        if log_likelihood > 700:  # exp(700) is near float64 limit
            return 1e300  # Large but finite value
        elif log_likelihood < -700:
            return 1e-300  # Small but non-zero value
        else:
            return np.exp(log_likelihood)
            
    except (AttributeError, OverflowError, ValueError):
        # Fallback to simpler calculation with overflow protection
        try:
            # Try to compute directly with smaller numbers
            if d1 > 1000 or d2 > 1000:
                # For very large depths, use approximation
                # Binomial -> Normal approximation when n is large
                from scipy.stats import norm
                
                # VAF expected value and variance for each marker
                mu1, var1 = f * d1, f * (1 - f) * d1
                mu2, var2 = f * d2, f * (1 - f) * d2
                
                # Use normal approximation to binomial
                prob1 = norm.pdf(r1, mu1, np.sqrt(var1)) if var1 > 0 else 1e-300
                prob2 = norm.pdf(r2, mu2, np.sqrt(var2)) if var2 > 0 else 1e-300
                
                return prob1 * prob2
            else:
                # For smaller depths, use original calculation
                return ncr(d1, r1) * (f ** r1) * ((1 - f) ** (d1 - r1)) * ncr(d2, r2) * (f ** r2) * ((1 - f) ** (d2 - r2))
                
        except (OverflowError, ValueError):
            # Ultimate fallback
            return 1e-300

def outer_integral_single(d1, d2, r1, r2, lower_bound, upper_bound):
    """
    Compute the integral with robust error handling for numerical issues.
    """
    try:
        integral, error = quad(integrand_single, f1_start(lower_bound), f1_end(upper_bound), 
                              args=(d1, d2, r1, r2), limit=50, epsabs=1e-8, epsrel=1e-6)
        
        # Check for reasonable result
        if np.isnan(integral) or np.isinf(integral) or integral <= 0:
            # Fallback to simpler approximation
            integral = 1e-10
            
        return integral
        
    except (OverflowError, ValueError, RuntimeError) as e:
        # If integration fails, return a small positive value
        # This allows the algorithm to continue with reduced confidence in this tree
        print(f"Integration warning: {e}. Using fallback value.")
        return 1e-10

def f2_start(f1, relation, lower_bound):
    if relation == 'descendant':
        return f1
    elif relation == 'ancestor' or 'null':
        return lower_bound
    elif relation == 'same':
        raise AttributeError
def f2_end(f1, relation, upper_bound):
    if relation == 'descendant':
        return upper_bound
    elif relation == 'ancestor':
        return f1
    elif relation == 'null':
        return 1 - f1
    elif relation == 'same':
        raise AttributeError


def f1_start(lower_bound=0):
    return lower_bound

def f1_end(upper_bound=1):
    return upper_bound


def update_tree_distribution(tree_distribution, accepted_tree_indices):
    updated_tree_distribution = {}
    for key, item_list in tree_distribution.items():
        print(key, len(item_list))
        updated_item_list = [item_list[idx] for idx in accepted_tree_indices]
        updated_tree_distribution[key] = updated_item_list
    return updated_tree_distribution
def update_tree_distribution_bayesian(tree_distribution, update_tree_freq_list):
    updated_tree_distribution = {}
    for key, item_list in tree_distribution.items():
        print(key, len(item_list))
        if key != 'freq':
            updated_tree_distribution[key] = item_list
        else:
            updated_tree_distribution['freq'] = update_tree_freq_list
    return updated_tree_distribution

def mut2node(node_dict):
    mut2node_dict = {}
    for node, mut_list in node_dict.items():
        for mut in mut_list:
            mut2node_dict[mut] = int(node)
    return mut2node_dict

def find_root(tree):
    non_root = []
    for item in tree.values():
        non_root += list(item)
    for node in tree.keys():
        if node not in non_root:
            return node

def bfs_structure(tree):  # O(k)
    order = []
    root = find_root(tree)
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order

def ancestor2descendant(tree):
    order = bfs_structure(tree)
    a2d = np.zeros((len(order), len(order)))
    for node in order[::-1]:
        if node in tree.keys():
            for child in tree[node]:
                a2d[int(node)][int(child)] = 1
                a2d[int(node)] += a2d[int(child)]
    return a2d

def test_single_tree_struct(tree_structure, node_dict, read_depth_list, ddpcr_marker_counts, marker_idx2gene, alpha):
    n_markers = len(ddpcr_marker_counts)
    a2d_matrix = ancestor2descendant(tree_structure)
    mut2node_dict = mut2node(node_dict)
    for (marker_idx1, marker_idx2) in list(combinations(list(range(len(ddpcr_marker_counts))), 2)):
        node1, node2 = mut2node_dict[marker_idx2gene[marker_idx1]], mut2node_dict[marker_idx2gene[marker_idx2]]
        read_depth_sublist = [read_depth_list[marker_idx1], read_depth_list[marker_idx1]]

        freq_hat_1, freq_hat_2 = ddpcr_marker_counts[marker_idx1]/read_depth_sublist[0], ddpcr_marker_counts[marker_idx2]/read_depth_sublist[1]
        if node1 == node2:
            relation = 'same'
        elif a2d_matrix[node1, node2] == 1:
            relation = 'ancestor'
        elif a2d_matrix[node1, node2] == 0 and a2d_matrix[node2, node1] == 1:
            relation = 'descendant'
        else:
            relation = 'null'
        if relation == 'null':
            bool_reject = False
            continue

        # node_correct - 1 is due to the fact that the normal node 0 contains no mutations, thus being omitted
        bool_reject, W, z = wald_test(freq_hat_1, freq_hat_2, n_markers * (n_markers - 1) / 2, relation, read_depth_sublist, alpha=alpha)
        if bool_reject:
            print(relation, bool_reject)
            return bool_reject
        print(relation, bool_reject, W, -z)
    return bool_reject


def test_single_tree_frac(clonal_freq_dict, node_dict, read_depth_list, ddpcr_marker_counts, marker_idx2gene, alpha):
    n_markers = len(ddpcr_marker_counts)
    mut2node_dict = mut2node(node_dict)
    for marker_idx in range(n_markers):
        marker_count = ddpcr_marker_counts[marker_idx]
        node = mut2node_dict[marker_idx2gene[marker_idx]]
        read_depth = read_depth_list[marker_idx]
        clonal_freq = clonal_freq_dict[node][0][0]
        # node_correct - 1 is due to the fact that the normal node 0 contains no mutations, thus being omitted
        bool_reject, T, chi_square = chi_square_test(clonal_freq, read_depth, marker_count, n_markers, alpha=alpha)
        if bool_reject:
            print(bool_reject, T, chi_square)
            return bool_reject
        print(bool_reject, T, chi_square)
    return bool_reject

def chi_square_test(clonal_freq, read_depth, marker_count, correction_rate,alpha=0.05):
    expected_num = [clonal_freq * read_depth, read_depth *(1-clonal_freq)]
    T = (marker_count - expected_num[0])**2 / expected_num[0] + (read_depth - marker_count - expected_num[1])**2 / expected_num[1]
    chi_square = chi2.ppf(1-alpha/correction_rate, 1)
    if T > chi_square:
        return True, T, chi_square
    else:
        return False, T, chi_square

def wald_test(freq_hat_1, freq_hat_2, correction_rate, relation='ancestor', depth=[100,100], alpha=0.05):
    '''
    return True if reject
    '''
    print(freq_hat_1, freq_hat_2)
    assert relation in ['ancestor', 'descendant', 'same']
    if relation == 'ancestor':
        W = (freq_hat_2 - freq_hat_1) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2)/ depth[1]))
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'descendant':
        W = (freq_hat_1 - freq_hat_2) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2) / depth[1]))
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'same':
        W = np.abs((freq_hat_1 - freq_hat_2) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2) / depth[1])))
        z = norm.ppf(alpha / correction_rate / 2)
    #print(W)
    if W > - z:
        return True, W, - z
    else:
        return False, W, z


if __name__ == '__main__':
    ddpcr = [{
        'gene': 'MLH1', 'mut': 100.21, 'WT': 1080.75, 'liquid_biopsy_sample': 'i', 'patient':256
    },
        {
        'gene': 'TP53', 'mut': 204.6, 'WT': 837.65, 'liquid_biopsy_sample': 'i', 'patient':256
        }]
    n_markers = 2
    blood_sample_idx = 0
    df_ddpcr = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr["gene"][i] for i in range(len(df_ddpcr)) }
    directory = Path('/data/liquid_biopsy/MOONSHOT2/bootstrap')
    patient_num = 256
    method = 'phylowgs'
    tree_distribution_file = directory / f'{patient_num}_bootstrap_{method}/{method}_bootstrap_summary.pkl'
    with open(directory / tree_distribution_file, 'rb') as f:
        tree_distribution_summary = pickle.load(f)

    gene_list = []
    gene2idx = {}
    inter = pd.read_excel(directory / f"xlsx/Patient_{patient_num}_bootstrap.xlsx",
                          sheet_name='common_blood_tissue_no_germline', index_col=0)
    inter = inter[inter["Allele Frequency_x"] < 0.9]
    inter = inter[inter["Allele Frequency_y"] < 0.9]
    calls = inter
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    genename2idx = {inter["Gene"][i] : i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    gene_name_list = []
    gene_count = {}
    for i in range(inter.shape[0]):
        gene = calls.iloc[i, 0]
        if gene in gene_name_list:
            gene_count[gene] += 1
            gene = gene + '_' + str(gene_count[gene])
        else:
            gene_count[gene] = 1
        if not isinstance(gene, str):
            gene = calls.iloc[i, 1] + '_' + str(calls.iloc[i, 2])
        gene_name_list.append(gene)

    tree_list, node_list, node_name_list, tree_freq_list = tree_distribution_summary['tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']

    # scrub node_list
    node_list_scrub = []
    for node_dict in node_list:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_scrub.append(temp)

    adjust_algo = 'bayesian'
    ddpcr_marker_counts = list(df_ddpcr["mut"])
    read_depth_list = list(df_ddpcr["mut"] + df_ddpcr["WT"])
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(tree_list, node_name_list,tree_freq_list, read_depth_list, ddpcr_marker_counts, marker_idx2gene)
    updated_tree_distribution_summary = update_tree_distribution_bayesian(tree_distribution_summary,
                                                                          updated_tree_freq_list)

    tree_distribution_file_summary_updated = directory / f'{patient_num}_bootstrap_{method}/{method}_bootstrap_summary_updated_{adjust_algo}_{n_markers}_{blood_sample_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)
    print(tree_freq_list, updated_tree_freq_list)

    plt.figure(figsize=(6, 18))
    timepoint = [0, 1]
    for i in range(len(tree_freq_list)):
        if tree_freq_list[i] > updated_tree_freq_list[i]:
            color = "tab:blue"
        else:
            print(i, tree_freq_list[i], updated_tree_freq_list[i])
            color = "tab:orange"
        plt.plot(timepoint, [tree_freq_list[i], updated_tree_freq_list[i]], marker='o', linestyle='-', color=color, ms=5)

    # Add labels and a legend
    plt.xlabel('time points')
    plt.ylabel('tree weights')
    plt.xticks([0,1])
    plt.savefig(f"Patient_{patient_num}_weights_change.eps",bbox_inches='tight')