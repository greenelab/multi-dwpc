"""
Gene Ontology processing and hierarchy analysis.

This module provides functions for parsing GO OBO files, calculating hierarchy
metrics, and computing percent changes in gene annotations over time.

"""

import re
from collections import defaultdict, deque
from pathlib import Path
import pandas as pd
import numpy as np


def _parse_obo_file(obo_path):
    """
    Parse GO OBO file and extract term information.

    Parameters
    ----------
    obo_path : str or Path
        Path to the OBO file

    Returns
    -------
    dict
        Dictionary mapping GO ID to term information including
        parents (is_a relationships) and namespace
    """
    terms = {}
    current_term = None

    with open(obo_path, 'r') as f:
        for line in f:
            line = line.strip()

            if line == '[Term]':
                if current_term:
                    terms[current_term['id']] = current_term
                current_term = {
                    'id': None,
                    'parents': [],
                    'namespace': None
                }

            elif line.startswith('id:') and current_term is not None:
                current_term['id'] = line.split('id:')[1].strip()

            elif line.startswith('is_a:') and current_term is not None:
                parent_id = line.split('is_a:')[1].split('!')[0].strip()
                current_term['parents'].append(parent_id)

            elif line.startswith('namespace:') and current_term is not None:
                namespace = line.split('namespace:')[1].strip()
                current_term['namespace'] = namespace

            elif line.startswith('is_obsolete:') and current_term is not None:
                current_term['obsolete'] = True

        if current_term and current_term['id']:
            terms[current_term['id']] = current_term

    return terms


def _calculate_depths(terms, namespace_filter='biological_process'):
    """
    Calculate depth from root for each GO term using BFS.

    Parameters
    ----------
    terms : dict
        Dictionary of term information from _parse_obo_file
    namespace_filter : str
        Only calculate depths for terms in this namespace

    Returns
    -------
    dict
        Dictionary mapping GO ID to depth (distance from root)
    """
    filtered_terms = {
        go_id: term for go_id, term in terms.items()
        if term['namespace'] == namespace_filter
        and not term.get('obsolete', False)
    }

    roots = [
        go_id for go_id, term in filtered_terms.items()
        if not term['parents'] or all(
            p not in filtered_terms for p in term['parents']
        )
    ]

    children = defaultdict(list)
    for go_id, term in filtered_terms.items():
        for parent in term['parents']:
            if parent in filtered_terms:
                children[parent].append(go_id)

    depths = {}
    queue = deque([(root, 0) for root in roots])

    while queue:
        go_id, depth = queue.popleft()

        if go_id in depths and depths[go_id] <= depth:
            continue

        depths[go_id] = depth

        for child_id in children[go_id]:
            queue.append((child_id, depth + 1))

    return depths


def _get_normalized_depths(obo_path, namespace='biological_process'):
    """
    Calculate normalized depth scores for GO terms.

    Normalized depth = depth / max_depth, giving values in [0, 1]

    Parameters
    ----------
    obo_path : str or Path
        Path to the OBO file
    namespace : str
        GO namespace to analyze

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: go_id, depth, max_depth, normalized_depth
    """
    terms = _parse_obo_file(obo_path)
    depths = _calculate_depths(terms, namespace)

    if not depths:
        return pd.DataFrame(
            columns=['go_id', 'depth', 'max_depth', 'normalized_depth']
        )

    max_depth = max(depths.values())

    df = pd.DataFrame([
        {
            'go_id': go_id,
            'depth': depth,
            'max_depth': max_depth,
            'normalized_depth': depth / max_depth if max_depth > 0 else 0
        }
        for go_id, depth in depths.items()
    ])

    return df


def identify_leaf_terms(obo_path, namespace='biological_process'):
    """
    Identify leaf GO terms (terms with no children).

    Parameters
    ----------
    obo_path : str or Path
        Path to the OBO file
    namespace : str
        GO namespace to analyze

    Returns
    -------
    set
        Set of GO IDs that are leaf terms
    """
    terms = _parse_obo_file(obo_path)

    namespace_terms = {
        go_id for go_id, term in terms.items()
        if term['namespace'] == namespace
        and not term.get('obsolete', False)
    }

    parent_terms = {
        parent
        for go_id, term in terms.items()
        if term['namespace'] == namespace and not term.get('obsolete', False)
        for parent in term['parents']
        if parent in namespace_terms
    }

    leaf_terms = namespace_terms - parent_terms

    return leaf_terms


def identify_parent_terms(leaf_terms, obo_path,
                          namespace='biological_process'):
    """
    Identify direct parent terms of leaf nodes (1 level up).

    Parameters
    ----------
    leaf_terms : set
        Set of GO IDs that are leaf terms
    obo_path : str or Path
        Path to the OBO file
    namespace : str
        GO namespace to analyze

    Returns
    -------
    set
        Set of GO IDs that are direct parents of leaf terms
    """
    terms = _parse_obo_file(obo_path)

    parent_set = {
        parent
        for leaf_id in leaf_terms
        if leaf_id in terms
        for parent in terms[leaf_id].get('parents', [])
        if parent in terms
        and terms[parent]['namespace'] == namespace
        and not terms[parent].get('obsolete', False)
    }

    return parent_set


def add_hierarchy_metrics(gene_count_df, obo_path_2016, obo_path_2024,
                          go_id_column='go_id'):
    """
    Add hierarchy metrics to gene count comparison dataframe.

    Parameters
    ----------
    gene_count_df : pd.DataFrame
        DataFrame with GO terms and gene counts (from merge of 2016/2024)
    obo_path_2016 : str or Path
        Path to 2016 GO ontology OBO file
    obo_path_2024 : str or Path
        Path to 2024 GO ontology OBO file
    go_id_column : str
        Name of GO ID column

    Returns
    -------
    pd.DataFrame
        Input dataframe with added columns:
        - depth_2016: absolute depth in 2016 hierarchy
        - normalized_depth_2016: depth/max_depth in 2016
        - depth_2024: absolute depth in 2024 hierarchy
        - normalized_depth_2024: depth/max_depth in 2024
        - depth_change: change in absolute depth
        - norm_depth_change: change in normalized depth
    """
    depths_2016 = _get_normalized_depths(obo_path_2016)
    depths_2024 = _get_normalized_depths(obo_path_2024)

    depths_2016 = depths_2016.rename(columns={
        'depth': 'depth_2016',
        'normalized_depth': 'normalized_depth_2016'
    })
    depths_2024 = depths_2024.rename(columns={
        'depth': 'depth_2024',
        'normalized_depth': 'normalized_depth_2024'
    })

    depths_2016 = depths_2016.drop(columns=['max_depth'])
    depths_2024 = depths_2024.drop(columns=['max_depth'])

    result = gene_count_df.copy()
    result = result.merge(depths_2016, on=go_id_column, how='left')
    result = result.merge(depths_2024, on=go_id_column, how='left')

    result['depth_change'] = result['depth_2024'] - result['depth_2016']
    result['norm_depth_change'] = (
        result['normalized_depth_2024'] - result['normalized_depth_2016']
    )

    return result


def classify_genes_stable_added(genes_2016_df, genes_2024_df):
    """
    Classify GO-gene pairs as stable or added between 2016 and 2024.

    Classification is per GO-term-gene pair (a gene may be stable for one
    GO term but added for another):
    - Stable: GO-gene pair present in both 2016 AND 2024 annotations
    - Added: GO-gene pair only in 2024 (gene existed in Hetionet 2016,
             but was not annotated to this GO term)

    Parameters
    ----------
    genes_2016_df : pd.DataFrame
        2016 GO-gene pairs with columns: go_id, entrez_gene_id
    genes_2024_df : pd.DataFrame
        2024 GO-gene pairs with columns: go_id, entrez_gene_id

    Returns
    -------
    genes_2016_stable : pd.DataFrame
        2016 GO-gene pairs for stable genes only
    genes_2024_added : pd.DataFrame
        2024 GO-gene pairs for added genes only
    classification_summary : pd.DataFrame
        Per-GO term statistics with columns:
        - go_id: GO term identifier
        - n_genes_2016: total genes in 2016
        - n_genes_2024: total genes in 2024
        - n_stable: genes present in both years
        - n_added: genes only in 2024
        - n_removed: genes only in 2016
        - pct_stable: percent of 2024 genes that are stable
        - pct_added: percent of 2024 genes that are added
    """
    genes_2016_df = genes_2016_df.copy()
    genes_2024_df = genes_2024_df.copy()

    # Create unique pair identifiers
    genes_2016_df['pair_id'] = (
        genes_2016_df['go_id'] + '|' +
        genes_2016_df['entrez_gene_id'].astype(str)
    )
    genes_2024_df['pair_id'] = (
        genes_2024_df['go_id'] + '|' +
        genes_2024_df['entrez_gene_id'].astype(str)
    )

    # Identify stable and added pairs using set operations
    pairs_2016_set = set(genes_2016_df['pair_id'])
    pairs_2024_set = set(genes_2024_df['pair_id'])

    stable_pairs = pairs_2016_set & pairs_2024_set
    added_pairs = pairs_2024_set - pairs_2016_set

    # Filter datasets
    genes_2016_stable = genes_2016_df[
        genes_2016_df['pair_id'].isin(stable_pairs)
    ].drop(columns=['pair_id']).reset_index(drop=True)

    genes_2024_added = genes_2024_df[
        genes_2024_df['pair_id'].isin(added_pairs)
    ].drop(columns=['pair_id']).reset_index(drop=True)

    # Compute per-GO term statistics
    summary_2016 = genes_2016_df.groupby('go_id')['entrez_gene_id'].nunique()
    summary_2024 = genes_2024_df.groupby('go_id')['entrez_gene_id'].nunique()
    summary_stable = genes_2016_stable.groupby('go_id')['entrez_gene_id'].nunique()
    summary_added = genes_2024_added.groupby('go_id')['entrez_gene_id'].nunique()

    classification_summary = pd.DataFrame({
        'go_id': summary_2016.index.union(summary_2024.index),
    })
    classification_summary['n_genes_2016'] = classification_summary['go_id'].map(
        summary_2016
    ).fillna(0).astype(int)
    classification_summary['n_genes_2024'] = classification_summary['go_id'].map(
        summary_2024
    ).fillna(0).astype(int)
    classification_summary['n_stable'] = classification_summary['go_id'].map(
        summary_stable
    ).fillna(0).astype(int)
    classification_summary['n_added'] = classification_summary['go_id'].map(
        summary_added
    ).fillna(0).astype(int)
    classification_summary['n_removed'] = (
        classification_summary['n_genes_2016'] -
        classification_summary['n_stable']
    )
    classification_summary['pct_stable'] = (
        100 * classification_summary['n_stable'] /
        classification_summary['n_genes_2024'].replace(0, 1)
    )
    classification_summary['pct_added'] = (
        100 * classification_summary['n_added'] /
        classification_summary['n_genes_2024'].replace(0, 1)
    )

    return genes_2016_stable, genes_2024_added, classification_summary


def calculate_percent_change(gene_count_df, count_col_2016,
                             count_col_2024, go_id_column='go_id'):
    """
    Calculate percent change in gene counts between 2016 and 2024.

    Parameters
    ----------
    gene_count_df : pd.DataFrame
        DataFrame with GO terms and gene counts for both years
    count_col_2016 : str
        Column name for 2016 gene counts
    count_col_2024 : str
        Column name for 2024 gene counts
    go_id_column : str
        Name of GO ID column

    Returns
    -------
    pd.DataFrame
        Input dataframe with added column:
        - pct_change_genes: percent change from 2016 to 2024
    """
    result = gene_count_df.copy()

    result['pct_change_genes'] = (
        (result[count_col_2024] - result[count_col_2016]) /
        result[count_col_2016]
    )

    return result
