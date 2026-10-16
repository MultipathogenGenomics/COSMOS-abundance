from pathlib import Path
import subprocess
import tempfile
import unittest
import io
from contextlib import redirect_stdout, redirect_stderr
from unittest.mock import patch

from cosmos_abundance.cli import parse_args
from cosmos_abundance.io import clean_run_outputs, write_fasta
from cosmos_abundance.run import run_pipeline, execute_workflow


class OutputTests(unittest.TestCase):
    def test_progress_hides_routine_output_and_preserves_failure_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "workflow.log"
            for status in [0, 1]:
                stdout, stderr = io.StringIO(), io.StringIO()
                with patch("cosmos_abundance.run.subprocess.Popen") as popen:
                    process = popen.return_value.__enter__.return_value
                    process.stdout = iter(["tool diagnostic\n", "1 of 2 steps (50%) done\n", "useful error detail\n"])
                    process.wait.return_value = status
                    with redirect_stdout(stdout), redirect_stderr(stderr):
                        if status:
                            with self.assertRaisesRegex(ValueError, "full log"):
                                execute_workflow(["snakemake"], log)
                        else:
                            execute_workflow(["snakemake"], log)
                self.assertNotIn("tool diagnostic", stdout.getvalue())
                self.assertIn("50%", stdout.getvalue())
                self.assertIn("tool diagnostic", log.read_text())
                if status:
                    self.assertIn("useful error detail", stderr.getvalue())
                else:
                    self.assertIn("100% complete", stdout.getvalue())

    def test_cleanup_removes_only_known_intermediates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in [".snakemake", "read-assignment", "relative_abundance"]:
                (root / name).mkdir()
            for name in ["config.json", "notes.txt", "relative_abundance/parameters.json",
                         "relative_abundance/filter_pass_read_counts.tsv"]:
                (root / name).write_text("test")
            clean_run_outputs(root, keep_intermediate=True)
            self.assertFalse((root / ".snakemake").exists())
            self.assertTrue((root / "read-assignment").is_dir())
            clean_run_outputs(root, keep_intermediate=False)
            self.assertFalse((root / "read-assignment").exists())
            self.assertFalse((root / "config.json").exists())
            self.assertFalse((root / "relative_abundance/filter_pass_read_counts.tsv").exists())
            self.assertTrue((root / "relative_abundance/parameters.json").is_file())
            self.assertTrue((root / "notes.txt").is_file())

    def test_failed_and_dry_runs_preserve_intermediates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fasta([{"id": "a", "sequence": "ACGT" * 150}], root / "ref.fa")
            for name in ["r1.fq", "r2.fq"]:
                (root / name).write_text("")
            (root / "samples.tsv").write_text("sample\tread1\tread2\ns\tr1.fq\tr2.fq\n")
            args = parse_args(["run", "--fasta", str(root / "ref.fa"), "--samples", str(root / "samples.tsv"),
                               "--out", str(root / "results")])
            with patch("cosmos_abundance.run.shutil.which", return_value="tool"), \
                 patch("cosmos_abundance.run.execute_workflow", side_effect=subprocess.CalledProcessError(1, "snakemake")), \
                 patch("cosmos_abundance.run.clean_run_outputs") as cleanup:
                with self.assertRaises(subprocess.CalledProcessError):
                    run_pipeline(args)
                cleanup.assert_not_called()
            self.assertTrue((root / "results/inputs.json").is_file())
            args.dry_run = True
            with patch("cosmos_abundance.run.shutil.which", return_value="tool"), \
                 patch("cosmos_abundance.run.subprocess.run") as command, \
                 patch("cosmos_abundance.run.clean_run_outputs") as cleanup:
                run_pipeline(args)
                self.assertIn("--dry-run", command.call_args.args[0])
                cleanup.assert_not_called()
