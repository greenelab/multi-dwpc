import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.result_normalization import (  # noqa: E402
    discover_year_result_files,
    load_normalized_year_results,
    summarize_normalized_year_results,
)


class ResultNormalizationTests(unittest.TestCase):
    def test_load_normalized_year_results_direct(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results_dir = Path(tmpdir)
            pd.DataFrame(
                [
                    {
                        "go_id": "GO:1",
                        "entrez_gene_id": 1,
                        "metapath": "BPpG<rG",
                        "dwpc": 0.2,
                    },
                    {
                        "go_id": "GO:1",
                        "entrez_gene_id": 2,
                        "metapath": "BPpG<rG",
                        "dwpc": 0.4,
                    },
                ]
            ).to_csv(results_dir / "dwpc_all_GO_positive_growth_2016_real.csv", index=False)

            norm_df = load_normalized_year_results(results_dir, score_source="direct")
            summary_df = summarize_normalized_year_results(norm_df)

        self.assertEqual(set(norm_df["score_source"].astype(str).unique().tolist()), {"direct"})
        self.assertEqual(set(norm_df["control"].astype(str).unique().tolist()), {"real"})
        self.assertEqual(int(summary_df.loc[0, "year"]), 2016)
        self.assertAlmostEqual(float(summary_df.loc[0, "mean_score"]), 0.3, places=7)
        self.assertEqual(int(summary_df.loc[0, "n_pairs"]), 2)

    def test_load_normalized_year_results_api(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results_dir = root / "results"
            data_dir = root / "data"
            results_dir.mkdir(parents=True, exist_ok=True)
            data_dir.mkdir(parents=True, exist_ok=True)

            pd.DataFrame(
                [
                    {"neo4j_id": 1001, "identifier": 1},
                    {"neo4j_id": 1002, "identifier": 2},
                ]
            ).to_csv(data_dir / "neo4j_gene_mapping.csv", index=False)
            pd.DataFrame([{"neo4j_id": 2001, "identifier": "GO:1"}]).to_csv(
                data_dir / "neo4j_bp_mapping.csv",
                index=False,
            )
            pd.DataFrame(
                [
                    {
                        "neo4j_source_id": 2001,
                        "neo4j_target_id": 1001,
                        "metapath_abbreviation": "BPpG<rG",
                        "dwpc": 0.5,
                    },
                    {
                        "neo4j_source_id": 2001,
                        "neo4j_target_id": 1002,
                        "metapath_abbreviation": "BPpG<rG",
                        "dwpc": 0.7,
                    },
                ]
            ).to_csv(results_dir / "res_all_GO_positive_growth_2024_perm_001.csv", index=False)

            files = discover_year_result_files(results_dir, score_source="api")
            norm_df = load_normalized_year_results(
                results_dir,
                score_source="api",
                data_dir=data_dir,
            )
            summary_df = summarize_normalized_year_results(norm_df)

        self.assertEqual(len(files), 1)
        self.assertEqual(set(norm_df["score_source"].astype(str).unique().tolist()), {"api"})
        self.assertEqual(set(norm_df["control"].astype(str).unique().tolist()), {"permuted"})
        self.assertEqual(int(norm_df.loc[0, "year"]), 2024)
        self.assertEqual(int(summary_df.loc[0, "replicate"]), 1)
        self.assertAlmostEqual(float(summary_df.loc[0, "mean_score"]), 0.6, places=7)


if __name__ == "__main__":
    unittest.main()
