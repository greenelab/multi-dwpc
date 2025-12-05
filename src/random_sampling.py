"""
Random gene sampling for null distribution generation.

This module provides functions for generating random gene samples that match
the size distribution of real GO-Gene associations, for use in null
distribution analysis.

Major update 2025-01-11:
- Added hash-based deterministic seeding (fixes reproducibility issue)
- Added option to exclude genes with ANY GO annotation (true null)
- Added degree-matched sampling option
- Added permutation-based null distribution functions
"""

import pandas as pd
import numpy as np
import hashlib


def generate_random_gene_samples(go_gene_df, all_genes, go_id_col='go_id',
                                 gene_id_col='neo4j_target_id',
                                 random_state=42):
    """
    Generate random gene samples matching GO term sizes.

    DEPRECATED: Use generate_random_gene_samples_v2() for improved
    reproducibility and statistical validity.

    For each GO term, samples the same number of random genes as real
    associations.

    Parameters
    ----------
    go_gene_df : pd.DataFrame
        DataFrame with GO-Gene associations
    all_genes : list or array
        Complete list of all genes available for sampling
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')
    random_state : int
        Random seed for reproducibility (default 42)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: go_id, neo4j_pseudo_target_id
        containing random gene samples
    """
    np.random.seed(random_state)

    random_samples = []

    for go_id in go_gene_df[go_id_col].unique():
        go_genes = go_gene_df[
            go_gene_df[go_id_col] == go_id
        ][gene_id_col].values

        n_genes = len(go_genes)

        random_genes = np.random.choice(all_genes, size=n_genes,
                                       replace=False)

        for gene in random_genes:
            random_samples.append({
                go_id_col: go_id,
                'neo4j_pseudo_target_id': gene
            })

    return pd.DataFrame(random_samples)


def generate_random_gene_samples_v2(go_gene_df, all_genes,
                                    all_go_annotations=None,
                                    gene_degrees=None,
                                    match_degree=False,
                                    exclude_annotated_genes=True,
                                    go_id_col='go_id',
                                    gene_id_col='neo4j_target_id',
                                    source_id_col='neo4j_source_id',
                                    random_state=42):
    """
    Generate random gene samples with improved statistical controls.

    This function fixes critical reproducibility and validity issues:
    1. Uses hash-based deterministic seeds (reproducible across runs)
    2. Option to exclude genes with ANY GO annotation (true null)
    3. Option to match degree distribution (control for connectivity)

    Parameters
    ----------
    go_gene_df : pd.DataFrame
        DataFrame with GO-Gene associations
        Required columns: go_id_col, gene_id_col, source_id_col
    all_genes : list or array
        Complete gene universe for sampling
        IMPORTANT: Use consistent universe across all datasets for fair
        temporal comparisons (recommend 2016 gene set for all years)
    all_go_annotations : pd.DataFrame, optional
        All GO-gene associations across entire dataset
        Required columns: go_id_col, gene_id_col
        If exclude_annotated_genes=True, genes in this dataframe are excluded
    gene_degrees : pd.DataFrame, optional
        Gene connectivity degrees
        Required columns: gene_id_col, 'degree'
        Used if match_degree=True
    match_degree : bool
        If True, sample random genes matching degree distribution
        Requires gene_degrees to be provided (default False)
    exclude_annotated_genes : bool
        If True, exclude genes with ANY GO annotation from sampling
        Creates true null distribution (default True)
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')
    source_id_col : str
        Column name for GO source IDs (default 'neo4j_source_id')
    random_state : int
        Base random seed for reproducibility (default 42)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - go_id_col: GO term ID
        - source_id_col: Neo4j source node ID
        - 'neo4j_pseudo_target_id': Random gene ID

    Notes
    -----
    Deterministic seeding: Each GO term gets a unique deterministic seed
    based on hash(random_state + go_id). This ensures reproducibility
    across runs regardless of iteration order.

    Examples
    --------
    # Basic usage (improved null distribution)
    random_df = generate_random_gene_samples_v2(
        go_gene_df=hetio_BPpG_filtered,
        all_genes=all_available_genes,
        all_go_annotations=all_go_gene_pairs,
        exclude_annotated_genes=True
    )

    # Degree-matched sampling
    random_df = generate_random_gene_samples_v2(
        go_gene_df=hetio_BPpG_filtered,
        all_genes=all_available_genes,
        gene_degrees=gene_degree_df,
        match_degree=True
    )
    """
    if match_degree and gene_degrees is None:
        raise ValueError("gene_degrees required when match_degree=True")

    all_genes = np.array(all_genes)

    if exclude_annotated_genes and all_go_annotations is not None:
        annotated_genes = set(all_go_annotations[gene_id_col].unique())
        all_genes = np.array([g for g in all_genes
                             if g not in annotated_genes])
        print(f"Excluded {len(annotated_genes)} annotated genes from "
              f"sampling pool")
        print(f"Sampling universe: {len(all_genes)} genes")

    if match_degree and gene_degrees is not None:
        gene_degree_dict = dict(zip(gene_degrees[gene_id_col],
                                   gene_degrees['degree']))
        all_genes_with_degree = [g for g in all_genes
                                 if g in gene_degree_dict]
        if len(all_genes_with_degree) < len(all_genes):
            print(f"Warning: {len(all_genes) - len(all_genes_with_degree)} "
                  f"genes missing degree information")
        all_genes = np.array(all_genes_with_degree)

    random_samples = []

    for go_id in go_gene_df[go_id_col].unique():
        go_subset = go_gene_df[go_gene_df[go_id_col] == go_id]
        go_genes = go_subset[gene_id_col].values
        n_genes = len(go_genes)

        if n_genes > len(all_genes):
            raise ValueError(
                f"GO term {go_id} has {n_genes} genes but sampling "
                f"universe only has {len(all_genes)} genes. "
                f"Cannot sample without replacement."
            )

        seed_string = f"{random_state}_{go_id}"
        seed_hash = hashlib.md5(seed_string.encode()).hexdigest()
        seed = int(seed_hash, 16) % (2**32)
        rng = np.random.RandomState(seed)

        if match_degree:
            real_degrees = [gene_degree_dict.get(g, 0) for g in go_genes]
            degree_quartiles = np.percentile(real_degrees, [25, 50, 75])

            sampled_genes = []
            for degree in real_degrees:
                if degree <= degree_quartiles[0]:
                    bin_genes = [g for g in all_genes
                                if gene_degree_dict.get(g, 0) <=
                                degree_quartiles[0]]
                elif degree <= degree_quartiles[1]:
                    bin_genes = [g for g in all_genes
                                if degree_quartiles[0] <
                                gene_degree_dict.get(g, 0) <=
                                degree_quartiles[1]]
                elif degree <= degree_quartiles[2]:
                    bin_genes = [g for g in all_genes
                                if degree_quartiles[1] <
                                gene_degree_dict.get(g, 0) <=
                                degree_quartiles[2]]
                else:
                    bin_genes = [g for g in all_genes
                                if gene_degree_dict.get(g, 0) >
                                degree_quartiles[2]]

                if len(bin_genes) > 0:
                    sampled_gene = rng.choice(bin_genes, size=1)[0]
                    sampled_genes.append(sampled_gene)
                else:
                    sampled_gene = rng.choice(all_genes, size=1)[0]
                    sampled_genes.append(sampled_gene)

            random_genes = np.array(sampled_genes)
        else:
            random_genes = rng.choice(all_genes, size=n_genes, replace=False)

        source_id = go_subset[source_id_col].iloc[0]

        for gene in random_genes:
            random_samples.append({
                go_id_col: go_id,
                source_id_col: source_id,
                'neo4j_pseudo_target_id': gene
            })

    return pd.DataFrame(random_samples)


def permute_go_labels(go_gene_df, n_permutations=100,
                     go_id_col='go_id',
                     gene_id_col='neo4j_target_id',
                     source_id_col='neo4j_source_id',
                     random_state=42):
    """
    Generate permuted GO-gene associations for null distribution.

    This approach permutes gene assignments among GO terms while preserving:
    - GO term metadata (neo4j_source_id, go_name, etc.) - CRITICAL for API
    - Gene degree distribution (same genes, shuffled labels)
    - Number of genes per GO term
    - Total number of associations

    This is the RECOMMENDED approach for mixed-effects models as it
    preserves biological structure while breaking GO-gene associations.

    IMPORTANT: This function preserves the GO term -> neo4j_source_id mapping
    so that permuted datasets can be used with the Hetionet API. Only gene
    assignments are shuffled, not GO term identities.

    Parameters
    ----------
    go_gene_df : pd.DataFrame
        Real GO-gene associations
        Required columns: go_id_col, gene_id_col, source_id_col
    n_permutations : int
        Number of permutations to generate (default 100)
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')
    source_id_col : str
        Column name for GO source IDs (default 'neo4j_source_id')
    random_state : int
        Random seed for reproducibility (default 42)

    Returns
    -------
    list of pd.DataFrame
        List of permuted dataframes, one per permutation
        Each has same structure as input go_gene_df

    Notes
    -----
    Algorithm:
    1. Extract GO term metadata (go_id, neo4j_source_id, go_name, etc.)
    2. For each permutation:
       a. Shuffle gene IDs among all GO terms
       b. Preserve GO term sizes
       c. Merge with GO metadata to restore neo4j_source_id

    This preserves gene-level characteristics (degree, other annotations)
    while testing whether specific GO-gene associations matter.

    BUG FIX (2025-12-02): Previous version shuffled go_id column but left
    neo4j_source_id unchanged, creating invalid GO-gene pairs that couldn't
    be used with Hetionet API. New version preserves GO term identities and
    shuffles only gene assignments.

    Examples
    --------
    # Generate 100 permutations
    permuted_dfs = permute_go_labels(
        hetio_BPpG_filtered,
        n_permutations=100,
        random_state=42
    )

    # Use in DWPC pipeline
    for i, perm_df in enumerate(permuted_dfs):
        # Save for DWPC computation
        perm_df.to_csv(f'permutation_{i:03d}.csv', index=False)
    """
    go_term_sizes = go_gene_df.groupby(go_id_col).size().to_dict()
    unique_go_ids = list(go_term_sizes.keys())

    # Extract GO term metadata to preserve neo4j_source_id and other columns
    # Each GO term must keep its correct neo4j_source_id for API compatibility
    metadata_cols = [col for col in go_gene_df.columns if col != gene_id_col]
    go_metadata = go_gene_df[metadata_cols].drop_duplicates(subset=[go_id_col])

    # Verify each GO term has unique metadata
    if len(go_metadata) != len(unique_go_ids):
        raise ValueError(
            f"GO term metadata is not unique: expected {len(unique_go_ids)} "
            f"GO terms but found {len(go_metadata)} unique metadata rows. "
            f"Each GO term must have consistent neo4j_source_id and other metadata."
        )

    # Extract all genes
    all_genes = go_gene_df[gene_id_col].values

    permuted_dataframes = []

    for perm_idx in range(n_permutations):
        perm_seed = random_state + perm_idx
        rng = np.random.RandomState(perm_seed)

        # Shuffle genes
        shuffled_genes = rng.permutation(all_genes)

        # Assign shuffled genes to GO terms while preserving sizes
        assignments = []
        gene_idx = 0

        for go_id in unique_go_ids:
            n_genes = go_term_sizes[go_id]
            assigned_genes = shuffled_genes[gene_idx:gene_idx + n_genes]

            for gene in assigned_genes:
                assignments.append({
                    go_id_col: go_id,
                    gene_id_col: gene
                })

            gene_idx += n_genes

        # Create dataframe from gene assignments
        perm_df = pd.DataFrame(assignments)

        # Merge with GO metadata to restore neo4j_source_id and other columns
        perm_df = perm_df.merge(go_metadata, on=go_id_col, how='left')

        # Reorder columns to match original
        perm_df = perm_df[go_gene_df.columns]

        permuted_dataframes.append(perm_df)

    return permuted_dataframes


def calculate_gene_promiscuity(go_gene_df, go_id_col='go_id',
                               gene_id_col='neo4j_target_id'):
    """
    Calculate gene promiscuity (number of GO terms per gene).

    Parameters
    ----------
    go_gene_df : pd.DataFrame
        GO-gene associations
        Required columns: go_id_col, gene_id_col
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - gene_id_col: Gene ID
        - 'promiscuity': Number of GO terms the gene belongs to
    """
    gene_promiscuity = go_gene_df.groupby(
        gene_id_col
    )[go_id_col].nunique().reset_index()

    gene_promiscuity.columns = [gene_id_col, 'promiscuity']

    return gene_promiscuity


def _sample_single_gene(go_id, real_gene, real_prom, candidate_pool,
                       promiscuity_tolerance, random_state, gene_id_col):
    """
    Helper function to sample one gene matching promiscuity.

    Separated for clarity and to avoid nested function definitions.

    Parameters
    ----------
    go_id : str
        GO term ID for deterministic seeding
    real_gene : int
        Real gene ID for deterministic seeding
    real_prom : int
        Promiscuity of the real gene (number of GO terms)
    candidate_pool : pd.DataFrame
        DataFrame of candidate genes with promiscuity values
        Required columns: gene_id_col, 'promiscuity'
    promiscuity_tolerance : int
        Allowed difference in promiscuity count
    random_state : int
        Base random seed
    gene_id_col : str
        Column name for gene IDs in candidate_pool

    Returns
    -------
    dict
        Dictionary with 'gene_id' and 'promiscuity' keys
    """
    prom_diff = np.abs(candidate_pool['promiscuity'] - real_prom)
    candidates = candidate_pool[prom_diff <= promiscuity_tolerance]

    tol = promiscuity_tolerance
    while len(candidates) == 0 and tol <= 20:
        tol += 2
        candidates = candidate_pool[prom_diff <= tol]

    if len(candidates) == 0:
        candidates = candidate_pool

    seed_string = f"{random_state}_{go_id}_{real_gene}"
    seed_hash = hashlib.md5(seed_string.encode()).hexdigest()
    seed = int(seed_hash, 16) % (2**32)
    rng = np.random.RandomState(seed)

    idx = rng.choice(len(candidates))
    sampled_row = candidates.iloc[idx]

    return {
        'gene_id': sampled_row[gene_id_col],
        'promiscuity': sampled_row['promiscuity']
    }


def generate_promiscuity_controlled_samples(go_gene_df, all_go_annotations,
                                            go_id_col='go_id',
                                            gene_id_col='neo4j_target_id',
                                            source_id_col='neo4j_source_id',
                                            promiscuity_tolerance=2,
                                            random_state=42):
    """
    Generate random gene samples controlling for gene promiscuity.

    For each GO-gene pair, samples a random gene from other GO terms with
    similar promiscuity (number of GO term annotations).

    Parameters
    ----------
    go_gene_df : pd.DataFrame
        Target GO-gene associations
        Required columns: go_id_col, gene_id_col, source_id_col
    all_go_annotations : pd.DataFrame
        All GO-gene associations for promiscuity calculation
        Required columns: go_id_col, gene_id_col
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')
    source_id_col : str
        Column name for GO source IDs (default 'neo4j_source_id')
    promiscuity_tolerance : int
        Allowed promiscuity difference (default 2)
    random_state : int
        Base random seed (default 42)

    Returns
    -------
    pd.DataFrame
        Random samples with columns:
        - go_id_col, source_id_col
        - 'neo4j_pseudo_target_id': Sampled gene
        - 'real_promiscuity', 'sampled_promiscuity'

    Notes
    -----
    Uses hash-based seeding for reproducibility independent of iteration
    order (see REPOSITORY_REVIEW.md Section 6.1).
    """
    promiscuity_df = calculate_gene_promiscuity(
        all_go_annotations, go_id_col, gene_id_col
    )

    df_with_prom = go_gene_df.merge(
        promiscuity_df, on=gene_id_col, how='left'
    )
    df_with_prom['promiscuity'] = df_with_prom['promiscuity'].fillna(1)

    go_gene_sets = go_gene_df.groupby(go_id_col)[gene_id_col].apply(
        set
    ).to_dict()

    all_results = []

    grouped = df_with_prom.groupby(go_id_col)
    for go_id, go_group in grouped:
        genes_in_term = go_gene_sets[go_id]
        candidate_pool = promiscuity_df[
            ~promiscuity_df[gene_id_col].isin(genes_in_term)
        ].copy()

        source_id = go_group[source_id_col].iloc[0]

        sampled = go_group.apply(
            lambda row: _sample_single_gene(
                go_id, row[gene_id_col], row['promiscuity'],
                candidate_pool, promiscuity_tolerance, random_state,
                gene_id_col
            ),
            axis=1
        )

        go_results = pd.DataFrame({
            go_id_col: go_id,
            source_id_col: source_id,
            'neo4j_pseudo_target_id': sampled.apply(lambda x: x['gene_id']),
            'real_promiscuity': go_group['promiscuity'].values,
            'sampled_promiscuity': sampled.apply(lambda x: x['promiscuity'])
        })

        all_results.append(go_results)

    return pd.concat(all_results, ignore_index=True)


def validate_random_samples(real_df, random_df, gene_degrees=None,
                            all_go_annotations=None,
                            go_id_col='go_id',
                            gene_id_col='neo4j_target_id'):
    """
    Validate random gene samples against expected properties.

    Checks:
    1. Sample sizes match (same N per GO term)
    2. Degree distributions similar (if gene_degrees provided)
    3. GO annotation enrichment (random < real)
    4. No overlap between real and random genes per GO term

    Parameters
    ----------
    real_df : pd.DataFrame
        Real GO-gene associations
    random_df : pd.DataFrame
        Random GO-gene samples
    gene_degrees : pd.DataFrame, optional
        Gene degree information
    all_go_annotations : pd.DataFrame, optional
        All GO annotations for enrichment check
    go_id_col : str
        Column name for GO IDs (default 'go_id')
    gene_id_col : str
        Column name for gene IDs (default 'neo4j_target_id')

    Returns
    -------
    dict
        Validation results with pass/fail for each check
    """
    results = {}

    real_sizes = real_df.groupby(go_id_col).size().sort_index()
    random_sizes = random_df.groupby(go_id_col).size().sort_index()
    sizes_match = (real_sizes == random_sizes).all()
    results['sample_sizes_match'] = sizes_match

    if gene_degrees is not None:
        gene_degree_dict = dict(zip(gene_degrees[gene_id_col],
                                   gene_degrees['degree']))
        real_degrees = [gene_degree_dict.get(g, 0)
                       for g in real_df[gene_id_col]]
        random_degrees = [gene_degree_dict.get(g, 0)
                         for g in random_df[gene_id_col]]

        from scipy import stats
        ks_stat, ks_pval = stats.ks_2samp(real_degrees, random_degrees)
        results['degree_distribution_similar'] = ks_pval > 0.05
        results['degree_ks_statistic'] = ks_stat
        results['degree_ks_pvalue'] = ks_pval

    if all_go_annotations is not None:
        real_genes = set(real_df[gene_id_col])
        random_genes = set(random_df[gene_id_col])
        annotated_genes = set(all_go_annotations[gene_id_col])

        real_annotated_pct = len(real_genes & annotated_genes) / len(
            real_genes)
        random_annotated_pct = len(random_genes & annotated_genes) / len(
            random_genes)

        results['real_annotated_percent'] = real_annotated_pct * 100
        results['random_annotated_percent'] = random_annotated_pct * 100
        results['random_less_annotated'] = (random_annotated_pct <
                                           real_annotated_pct)

    overlap_count = 0
    for go_id in real_df[go_id_col].unique():
        real_genes_go = set(real_df[real_df[go_id_col] == go_id][
            gene_id_col])
        random_genes_go = set(random_df[random_df[go_id_col] == go_id][
            'neo4j_pseudo_target_id'])
        overlap = len(real_genes_go & random_genes_go)
        overlap_count += overlap

    results['total_overlap_genes'] = overlap_count
    results['no_overlap'] = overlap_count == 0

    return results
