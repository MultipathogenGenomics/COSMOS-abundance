import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import polars as pl


@unittest.skipUnless(os.environ.get("COSMOS_INTEGRATION") == "1", "set COSMOS_INTEGRATION=1 to run MAFFT/KMA/kallisto")
class IntegrationTests(unittest.TestCase):
    def test_prepare_run_resume_and_empty_run(self):
        example = Path(__file__).resolve().parents[1] / "examples/make_example.py"
        spec = importlib.util.spec_from_file_location("example", example)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory(prefix="cosmos integration ") as directory:
            root = Path(directory)
            module.make_example(root / "data")

            def run(*args):
                result = subprocess.run([sys.executable, "-m", "cosmos_abundance", *map(str, args)],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                return result

            run("prepare", "--fasta", root / "data/references.fasta", "--species", root / "data/species.tsv",
                "--out", root / "reference")
            assignments = pl.read_csv(root / "reference/assignments.tsv", separator="\t")
            self.assertEqual(set(assignments["assignment_name"]), {"Species alpha + Species beta", "Species gamma"})
            for extension in ["png", "svg"]:
                self.assertGreater((root / f"reference/mash_tree.{extension}").stat().st_size, 100)
            command = ["run", "--reference", root / "reference", "--samples", root / "data/samples.tsv",
                       "--out", root / "results", "--cores", "2", "--bootstraps", "3"]
            run(*command)
            table_path = root / "results/abundance/abundance_within_filter_pass_reads.tsv"
            table = pl.read_csv(table_path, separator="\t")
            mixture = table.filter(pl.col("sample") == "mixture")
            self.assertAlmostEqual(mixture["Species alpha + Species beta"].item(), 0.75, delta=0.05)
            self.assertEqual(table.filter(pl.col("sample") == "unassigned")["Species gamma"].item(), 0)
            modified = table_path.stat().st_mtime_ns
            run(*command)
            self.assertEqual(table_path.stat().st_mtime_ns, modified)
            run("estimate", "--run-dir", root / "results", "--min-kma-count", "1000000")
            self.assertEqual(pl.read_csv(table_path, separator="\t").columns, ["sample"])
            run(*command)
            self.assertIn("Species gamma", pl.read_csv(table_path, separator="\t").columns)

            # An entirely unassigned batch, using the FASTA bypass interface.
            (root / "data/empty.tsv").write_text("sample\tread1\tread2\n001\tunassigned_R1.fastq.gz\tunassigned_R2.fastq.gz\n")
            run("run", "--fasta", root / "reference/references.fasta", "--samples", root / "data/empty.tsv",
                "--out", root / "empty-results", "--cores", "2", "--bootstraps", "3")
            empty = pl.read_csv(root / "empty-results/abundance/abundance_within_all_assigned_reads.tsv",
                                separator="\t", schema_overrides={"sample": pl.String})
            self.assertEqual(empty["sample"].to_list(), ["001"])
            self.assertEqual(empty["non_confidently_assigned"].item(), 0)


if __name__ == "__main__":
    unittest.main()
