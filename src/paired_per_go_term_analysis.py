"""
Paired Per-GO-Term DWPC Analysis

This module implements the corrected statistical analysis approach for
comparing DWPC values across datasets. The key insight is that GO terms are
the experimental units, not individual gene pairs.

Statistical Approach:
- Aggregates DWPC by (GO term, metapath, dataset)
- Performs paired statistical tests (same GO terms across datasets)
- Computes effect sizes for paired data using Cohen's d

Why Paired Analysis:
- Comparing same GO terms across conditions (real, permuted, random)
- Paired tests have higher statistical power than unpaired tests
- Controls for GO-term-specific confounds (size, hierarchy level, etc.)

Usage:
    from src.paired_per_go_term_analysis import run_paired_analysis

    datasets = {
        '2016_real': df1,
        '2016_perm': df2,
        '2016_random': df3,
        '2024_real': df4,
        '2024_random': df5
    }

    agg_df = aggregate_by_go_term_metapath(datasets)

    comparisons = [
        {'year': 2016, 'type': 'real_vs_perm', 'name': '2016_real_vs_perm'},
        {'year': 2016, 'type': 'real_vs_random', 'name': '2016_real_vs_random'},
        {'year': 2024, 'type': 'real_vs_random', 'name': '2024_real_vs_random'}
    ]

    results_df = run_paired_analysis(agg_df, comparisons)
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path


def aggregate_by_go_term_metapath(datasets):
    """
    Aggregate DWPC values by (GO term, metapath, dataset).

    This function groups individual (GO term, Gene) DWPC values by GO term
    and computes summary statistics. This is essential for proper statistical
    analysis because GO terms are the experimental units, not individual genes.

    Parameters
    ----------
    datasets : dict
        Dictionary mapping dataset labels (e.g., '2016_real') to DataFrames.
        Each DataFrame must have columns:
        - neo4j_source_id: GO term identifier
        - metapath_abbreviation: Metapath string
        - dwpc: DWPC value

    Returns
    -------
    pd.DataFrame
        Aggregated data with columns:
        - neo4j_source_id: GO term identifier
        - metapath_abbreviation: Metapath string
        - dataset_label: Dataset name
        - dwpc_mean: Mean DWPC for this (GO term, metapath, dataset)
        - dwpc_std: Standard deviation of DWPC
        - dwpc_count: Number of genes in this GO term
        - dwpc_nonzero_mean: Mean of nonzero DWPC values
        - dwpc_nonzero_count: Number of genes with nonzero DWPC
    """
    agg_results = []

    for label, df in datasets.items():
        grouped = df.groupby(['neo4j_source_id', 'metapath_abbreviation'])

        for (go_term, metapath), group in grouped:
            dwpc_vals = group['dwpc'].dropna()
            dwpc_nonzero = dwpc_vals[dwpc_vals > 0]

            agg_results.append({
                'neo4j_source_id': go_term,
                'metapath_abbreviation': metapath,
                'dataset_label': label,
                'dwpc_mean': dwpc_vals.mean() if len(dwpc_vals) > 0 else np.nan,
                'dwpc_std': dwpc_vals.std() if len(dwpc_vals) > 1 else np.nan,
                'dwpc_count': len(dwpc_vals),
                'dwpc_nonzero_mean': dwpc_nonzero.mean() if len(dwpc_nonzero) > 0 else np.nan,
                'dwpc_nonzero_count': len(dwpc_nonzero)
            })

    agg_df = pd.DataFrame(agg_results)

    print(f"Aggregated {len(agg_df):,} (GO term, metapath, dataset) combinations")
    print(f"  {agg_df['neo4j_source_id'].nunique()} unique GO terms")
    print(f"  {agg_df['metapath_abbreviation'].nunique()} unique metapaths")
    print(f"  {agg_df['dataset_label'].nunique()} datasets")

    return agg_df


def pivot_for_paired_comparison(agg_df, year, comparison_type):
    """
    Pivot aggregated data to wide format for paired comparisons.

    Creates matched pairs of GO terms across two datasets for paired
    statistical testing. Only GO terms present in both datasets are retained.

    Parameters
    ----------
    agg_df : pd.DataFrame
        Aggregated data from aggregate_by_go_term_metapath()
    year : int or None
        Year to analyze (2016 or 2024). Use None for temporal comparisons.
    comparison_type : str
        One of: 'real_vs_perm', 'real_vs_random', or 'real_2016_vs_2024'

    Returns
    -------
    pd.DataFrame
        Wide format with columns:
        - neo4j_source_id: GO term identifier
        - metapath_abbreviation: Metapath string
        - dwpc_real: Mean DWPC in first dataset
        - dwpc_comparison: Mean DWPC in second dataset

    Raises
    ------
    ValueError
        If comparison_type is not recognized
    """
    if comparison_type == 'real_2016_vs_2024':
        dataset_labels = ['2016_real', '2024_real']
        col_mapping = {
            '2016_real': 'dwpc_real',
            '2024_real': 'dwpc_comparison'
        }
        display_name = 'real_2016_vs_2024'
    elif comparison_type == 'real_vs_perm':
        dataset_labels = [f'{year}_real', f'{year}_perm']
        col_mapping = {
            f'{year}_real': 'dwpc_real',
            dataset_labels[1]: 'dwpc_comparison'
        }
        display_name = f'{year} {comparison_type}'
    elif comparison_type == 'real_vs_random':
        dataset_labels = [f'{year}_real', f'{year}_random']
        col_mapping = {
            f'{year}_real': 'dwpc_real',
            dataset_labels[1]: 'dwpc_comparison'
        }
        display_name = f'{year} {comparison_type}'
    else:
        raise ValueError(f"Unknown comparison_type: {comparison_type}")

    subset = agg_df[agg_df['dataset_label'].isin(dataset_labels)].copy()

    pivot_df = subset.pivot_table(
        index=['neo4j_source_id', 'metapath_abbreviation'],
        columns='dataset_label',
        values='dwpc_mean'
    ).reset_index()

    pivot_df.columns.name = None
    pivot_df = pivot_df.rename(columns=col_mapping)

    n_before = len(pivot_df)
    pivot_df = pivot_df.dropna(subset=['dwpc_real', 'dwpc_comparison'])
    n_after = len(pivot_df)
    n_dropped = n_before - n_after

    print(f"\n{display_name}:")
    print(f"  {n_after:,} GO term-metapath pairs with data in both datasets")
    if n_dropped > 0:
        print(f"  {n_dropped:,} pairs dropped (missing in one dataset)")

    return pivot_df


def cohens_d_paired(diff_values):
    """
    Compute Cohen's d effect size for paired samples.

    For paired data, effect size is computed from the differences:
    Cohen's d = mean(differences) / std(differences)

    This differs from unpaired Cohen's d which uses pooled standard deviation.

    Parameters
    ----------
    diff_values : array-like
        Array of paired differences (real - comparison)

    Returns
    -------
    float
        Cohen's d effect size. Returns np.nan if computation is not possible.

    Interpretation:
        Small effect: |d| ~ 0.2
        Medium effect: |d| ~ 0.5
        Large effect: |d| ~ 0.8
    """
    if len(diff_values) == 0:
        return np.nan

    mean_diff = np.mean(diff_values)
    std_diff = np.std(diff_values, ddof=1)

    if std_diff == 0:
        return np.nan

    return mean_diff / std_diff


def perform_paired_comparison_per_metapath(pivot_df, metapath, comparison_name, min_n=10):
    """
    Perform paired t-test for one metapath.

    Uses scipy.stats.ttest_rel for paired samples (same GO terms across
    datasets). This is appropriate because we are comparing the same GO terms
    in different conditions.

    Parameters
    ----------
    pivot_df : pd.DataFrame
        Pivoted data with dwpc_real and dwpc_comparison columns
    metapath : str
        Metapath abbreviation to test
    comparison_name : str
        Name of comparison (e.g., "2016_real_vs_perm")
    min_n : int, default=10
        Minimum number of paired observations required for testing

    Returns
    -------
    dict or None
        Test results dictionary with keys:
        - metapath: Metapath string
        - comparison: Comparison name
        - n_pairs: Number of paired GO terms
        - mean_real: Mean DWPC in real dataset
        - mean_comparison: Mean DWPC in comparison dataset
        - mean_diff: Mean of paired differences
        - std_diff: Standard deviation of paired differences
        - t_statistic: t-statistic from paired t-test
        - p_value: Two-tailed p-value
        - cohens_d: Effect size (paired Cohen's d)

        Returns None if insufficient data (n < min_n)
    """
    mp_data = pivot_df[pivot_df['metapath_abbreviation'] == metapath].copy()

    if len(mp_data) < min_n:
        return None

    real_vals = mp_data['dwpc_real'].values
    comp_vals = mp_data['dwpc_comparison'].values

    t_stat, p_val = stats.ttest_rel(real_vals, comp_vals)

    diff_vals = real_vals - comp_vals
    effect_size = cohens_d_paired(diff_vals)

    mean_real = np.mean(real_vals)
    mean_comp = np.mean(comp_vals)
    mean_diff = np.mean(diff_vals)
    std_diff = np.std(diff_vals, ddof=1)

    return {
        'metapath': metapath,
        'comparison': comparison_name,
        'n_pairs': len(mp_data),
        'mean_real': mean_real,
        'mean_comparison': mean_comp,
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        't_statistic': t_stat,
        'p_value': p_val,
        'cohens_d': effect_size
    }


def run_paired_analysis(agg_df, comparisons_config):
    """
    Run all paired comparisons and apply FDR correction.

    This function orchestrates the complete paired analysis workflow:
    1. Pivot data for each comparison
    2. Run paired t-tests for all metapaths
    3. Apply FDR correction (both per-comparison and global)
    4. Return comprehensive results

    Parameters
    ----------
    agg_df : pd.DataFrame
        Aggregated data from aggregate_by_go_term_metapath()
    comparisons_config : list of dict
        List of comparisons to run. Each dict must have:
        - year: int (2016 or 2024)
        - type: str ('real_vs_perm' or 'real_vs_random')
        - name: str (comparison identifier)

    Returns
    -------
    pd.DataFrame
        All paired test results with columns:
        - All columns from perform_paired_comparison_per_metapath()
        - p_fdr: FDR-corrected p-value (per comparison)
        - p_fdr_global: FDR-corrected p-value (across all tests)

    Notes
    -----
    FDR correction uses Benjamini-Hochberg procedure (statsmodels.stats.multitest)
    Both per-comparison and global FDR are provided for flexibility.
    """
    all_results = []

    print("\n" + "="*80)
    print("PAIRED STATISTICAL TESTS")
    print("="*80)

    for config in comparisons_config:
        year = config['year']
        comp_type = config['type']
        comp_name = config['name']

        print(f"\nComparison: {comp_name}")
        print("-" * 40)

        pivot_df = pivot_for_paired_comparison(agg_df, year, comp_type)

        if len(pivot_df) == 0:
            print(f"  No data available for {comp_name}. Skipping.")
            continue

        metapaths = sorted(pivot_df['metapath_abbreviation'].unique())
        print(f"  Testing {len(metapaths)} metapaths...")

        for metapath in metapaths:
            result = perform_paired_comparison_per_metapath(
                pivot_df, metapath, comp_name
            )
            if result is not None:
                all_results.append(result)

        n_tests = len([r for r in all_results if r['comparison'] == comp_name])
        print(f"  Completed {n_tests} paired tests")

    if len(all_results) == 0:
        print("\nWARNING: No tests completed!")
        return pd.DataFrame()

    results_df = pd.DataFrame(all_results)

    print("\n" + "="*80)
    print("APPLYING FDR CORRECTION")
    print("="*80)

    for comp_name in results_df['comparison'].unique():
        comp_mask = results_df['comparison'] == comp_name
        p_values = results_df.loc[comp_mask, 'p_value'].values

        _, p_fdr, _, _ = multipletests(p_values, method='fdr_bh')
        results_df.loc[comp_mask, 'p_fdr'] = p_fdr

        n_sig = (p_fdr < 0.05).sum()
        print(f"{comp_name}: {n_sig}/{len(p_values)} significant at FDR < 0.05")

    all_p_values = results_df['p_value'].values
    _, p_fdr_global, _, _ = multipletests(all_p_values, method='fdr_bh')
    results_df['p_fdr_global'] = p_fdr_global

    n_sig_global = (p_fdr_global < 0.05).sum()
    print(f"\nGlobal FDR: {n_sig_global}/{len(results_df)} significant at FDR < 0.05")

    return results_df


def summarize_paired_results(results_df, output_dir):
    """
    Create summary statistics and save results to disk.

    Parameters
    ----------
    results_df : pd.DataFrame
        Paired test results from run_paired_analysis()
    output_dir : str or Path
        Output directory for saving results

    Returns
    -------
    str
        Summary text with key statistics

    Side Effects
    ------------
    Creates two files in output_dir:
    - paired_tests_per_metapath.csv: Full results table
    - paired_tests_summary.txt: Text summary
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_df.to_csv(output_dir / 'paired_tests_per_metapath.csv', index=False)
    print(f"\nSaved: {output_dir / 'paired_tests_per_metapath.csv'}")

    summary_lines = []
    summary_lines.append("PAIRED STATISTICAL TESTS SUMMARY")
    summary_lines.append("="*80)
    summary_lines.append(f"Total paired tests: {len(results_df)}")
    summary_lines.append("")

    summary_lines.append("Significant results by comparison (FDR < 0.05):")
    for comp_name in results_df['comparison'].unique():
        comp_data = results_df[results_df['comparison'] == comp_name]
        n_sig = (comp_data['p_fdr'] < 0.05).sum()
        n_total = len(comp_data)
        summary_lines.append(f"  {comp_name}: {n_sig}/{n_total}")

    summary_lines.append("")
    summary_lines.append(f"Global FDR < 0.05: {(results_df['p_fdr_global'] < 0.05).sum()}/{len(results_df)}")

    summary_lines.append("")
    summary_lines.append("Effect size statistics (Cohen's d for paired data):")
    summary_lines.append(f"  Mean: {results_df['cohens_d'].mean():.4f}")
    summary_lines.append(f"  Median: {results_df['cohens_d'].median():.4f}")
    summary_lines.append(f"  Std: {results_df['cohens_d'].std():.4f}")
    summary_lines.append(f"  Range: [{results_df['cohens_d'].min():.4f}, {results_df['cohens_d'].max():.4f}]")

    summary_text = "\n".join(summary_lines)

    with open(output_dir / 'paired_tests_summary.txt', 'w') as f:
        f.write(summary_text)

    print(f"Saved: {output_dir / 'paired_tests_summary.txt'}")
    print("\n" + summary_text)

    return summary_text
