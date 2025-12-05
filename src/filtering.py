"""
GO term filtering functions.

This module provides functions for filtering GO terms based on IQR thresholds
and removing redundant terms with high Jaccard similarity.

"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def calculate_iqr_thresholds(df, gene_col='no_of_genes_in_hetio_GO_2016',
                             pct_col='pct_change_genes'):
    """
    Calculate IQR-based thresholds for filtering GO terms.

    Uses interquartile range to identify reasonable bounds for gene counts
    and percent change values.

    Parameters
    ----------
    df : pd.DataFrame
        GO term data with gene counts and percent changes
    gene_col : str
        Column name for gene counts
    pct_col : str
        Column name for percent change values

    Returns
    -------
    dict
        Dictionary with keys: min_genes, max_genes, median_genes,
        min_pct_change, max_pct_change, median_pct_change, total_terms
    """
    q1_genes = df[gene_col].quantile(0.25)
    q3_genes = df[gene_col].quantile(0.75)
    iqr_genes = q3_genes - q1_genes
    median_genes = df[gene_col].median()

    q1_pct = df[pct_col].quantile(0.25)
    q3_pct = df[pct_col].quantile(0.75)
    iqr_pct = q3_pct - q1_pct
    median_pct = df[pct_col].median()

    thresholds = {
        'min_genes': q1_genes - 1.5 * iqr_genes,
        'max_genes': q3_genes + 1.5 * iqr_genes,
        'median_genes': median_genes,
        'min_pct_change': q1_pct - 1.5 * iqr_pct,
        'max_pct_change': q3_pct + 1.5 * iqr_pct,
        'median_pct_change': median_pct,
        'total_terms': len(df)
    }

    return thresholds


def filter_go_terms_iqr(df, thresholds, min_genes=10,
                       gene_col='no_of_genes_in_hetio_GO_2016',
                       pct_col='pct_change_genes'):
    """
    Filter GO terms using IQR thresholds and minimum gene count.

    Parameters
    ----------
    df : pd.DataFrame
        GO term comparison data
    thresholds : dict
        Dictionary with min/max values from calculate_iqr_thresholds
    min_genes : int
        Minimum number of genes required (default 10)
    gene_col : str
        Column name for gene counts
    pct_col : str
        Column name for percent change values

    Returns
    -------
    pd.DataFrame
        Filtered dataframe containing only GO terms within IQR bounds
        and above minimum gene count
    """
    mask = (
        (df[gene_col] >= min_genes) &
        (df[gene_col] >= thresholds['min_genes']) &
        (df[gene_col] <= thresholds['max_genes']) &
        (df[pct_col] >= thresholds['min_pct_change']) &
        (df[pct_col] <= thresholds['max_pct_change'])
    )
    return df.loc[mask].reset_index(drop=True)


def filter_overlapping_go_terms(jaccard_matrix, dataset, threshold=0.9):
    """
    Remove redundant GO terms with Jaccard similarity > threshold.

    Uses graph-based clustering to handle transitive relationships.
    Selects representative with greatest absolute percent change in
    gene count.

    Parameters
    ----------
    jaccard_matrix : pd.DataFrame
        Symmetric matrix of Jaccard similarities (GO terms x GO terms)
    dataset : pd.DataFrame
        Dataset with columns: go_id, no_of_genes_in_hetio_GO_2016,
        no_of_genes_in_GO_2024
    threshold : float
        Jaccard similarity threshold for redundancy (default 0.9)

    Returns
    -------
    filtered_dataset : pd.DataFrame
        Dataset containing only representative GO terms
    removed_terms : dict
        Mapping {removed_go_id: representative_go_id}
    cluster_stats : pd.DataFrame
        Statistics about identified clusters
    removal_df : pd.DataFrame
        Detailed information about removed terms
    """
    adjacency = (jaccard_matrix > threshold).astype(int)
    np.fill_diagonal(adjacency.values, 0)

    n_components, labels = connected_components(
        csgraph=csr_matrix(adjacency.values),
        directed=False
    )

    clusters = {}
    for go_id, label in zip(jaccard_matrix.index, labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(go_id)

    representatives = []
    removed_terms = {}
    removal_details = []

    for cluster_id, cluster_members in clusters.items():
        if len(cluster_members) == 1:
            representatives.append(cluster_members[0])
        else:
            cluster_data = dataset[
                dataset['go_id'].isin(cluster_members)
            ].copy()

            cluster_data['pct_change'] = abs(
                (cluster_data['no_of_genes_in_GO_2024'] -
                 cluster_data['no_of_genes_in_hetio_GO_2016']) /
                cluster_data['no_of_genes_in_hetio_GO_2016'] * 100
            )

            max_pct_change = cluster_data['pct_change'].max()
            candidates = cluster_data[
                cluster_data['pct_change'] == max_pct_change
            ]
            representative = sorted(candidates['go_id'].tolist())[0]

            representatives.append(representative)

            non_reps = cluster_data[cluster_data['go_id'] != representative]
            for go_id in non_reps['go_id']:
                removed_terms[go_id] = representative

            if len(non_reps) > 0:
                rep_row = cluster_data[
                    cluster_data['go_id'] == representative
                ].iloc[0]

                non_reps = non_reps.copy()
                non_reps['representative_go_id'] = representative
                non_reps['jaccard_similarity'] = non_reps['go_id'].apply(
                    lambda x: jaccard_matrix.loc[x, representative]
                )
                non_reps['kept_genes_2016'] = rep_row['no_of_genes_in_hetio_GO_2016']
                non_reps['kept_genes_2024'] = rep_row['no_of_genes_in_GO_2024']
                non_reps['kept_pct_change'] = rep_row['pct_change']

                removal_details.extend(non_reps.rename(columns={
                    'go_id': 'removed_go_id',
                    'no_of_genes_in_hetio_GO_2016': 'removed_genes_2016',
                    'no_of_genes_in_GO_2024': 'removed_genes_2024',
                    'pct_change': 'removed_pct_change'
                })[[
                    'removed_go_id', 'representative_go_id', 'jaccard_similarity',
                    'removed_genes_2016', 'removed_genes_2024', 'removed_pct_change',
                    'kept_genes_2016', 'kept_genes_2024', 'kept_pct_change'
                ]].to_dict('records'))

    filtered_dataset = dataset[
        dataset['go_id'].isin(representatives)
    ].copy()

    cluster_sizes = [len(members) for members in clusters.values()]
    max_size = max(cluster_sizes)
    cluster_stats = pd.DataFrame({
        'cluster_size': range(1, max_size + 1),
        'count': [
            cluster_sizes.count(size) for size in range(1, max_size + 1)
        ]
    })
    cluster_stats = cluster_stats[cluster_stats['count'] > 0]

    removal_df = (
        pd.DataFrame(removal_details)
        if removal_details
        else pd.DataFrame()
    )

    return filtered_dataset, removed_terms, cluster_stats, removal_df
