"""
Similarity calculations for GO terms.

This module provides functions for calculating Jaccard similarity between GO
terms based on their gene annotations.
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import os


def calculate_jaccard_similarity_optimized(go_term_genes_dict):
    """
    Calculate Jaccard similarity matrix for GO terms.

    Uses scipy's pdist for efficient computation of pairwise similarities.

    Parameters
    ----------
    go_term_genes_dict : dict
        Dictionary mapping GO ID to set of gene IDs

    Returns
    -------
    pd.DataFrame
        Symmetric matrix with Jaccard similarities between all GO term pairs
    """
    go_terms = sorted(go_term_genes_dict.keys())
    all_genes = sorted(set().union(*go_term_genes_dict.values()))

    n_terms = len(go_terms)
    n_genes = len(all_genes)

    gene_to_idx = {gene: idx for idx, gene in enumerate(all_genes)}

    feature_matrix = np.zeros((n_terms, n_genes), dtype=np.bool_)
    for i, go_id in enumerate(go_terms):
        genes = go_term_genes_dict[go_id]
        for gene in genes:
            feature_matrix[i, gene_to_idx[gene]] = True

    condensed_distances = pdist(feature_matrix, metric='jaccard')
    distance_matrix = squareform(condensed_distances)
    similarities = 1 - distance_matrix

    return pd.DataFrame(similarities, index=go_terms, columns=go_terms)


def load_or_calculate_jaccard(go_term_genes_dict, cache_file):
    """
    Load cached Jaccard matrix or calculate and cache if not exists.

    Parameters
    ----------
    go_term_genes_dict : dict
        Dictionary mapping GO ID to set of gene IDs
    cache_file : str or Path
        Path to cache file for storing/loading Jaccard matrix

    Returns
    -------
    pd.DataFrame
        Jaccard similarity matrix
    """
    if os.path.exists(cache_file):
        print(f"Loading cached Jaccard matrix from {cache_file}")
        return pd.read_csv(cache_file, index_col=0)
    else:
        print(f"Computing Jaccard similarity matrix...")
        similarity_matrix = calculate_jaccard_similarity_optimized(
            go_term_genes_dict
        )
        print(f"Saving to {cache_file}")
        similarity_matrix.to_csv(cache_file)
        return similarity_matrix
