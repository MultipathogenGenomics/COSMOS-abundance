from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import shutil
import warnings

import numpy as np

from cosmos_abundance.io import parse_fasta_file, read_species_map, read_species_metadata
from cosmos_abundance.cli import parse_args
from cosmos_abundance.reference.prepare import prepare_reference
from cosmos_abundance.reference.plot import expected_figure_width, plot_mash_tree
from cosmos_abundance.reference.cluster import cluster_sequences_by_distance, identify_best_representatives
from cosmos_abundance.reference.fasta import find_best_gap_free_region, trim_alignment


class ReferenceTests(unittest.TestCase):
    def test_n_filter_runs_after_length_filter_before_alignment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input.fa").write_text(
                f">short\n{'N' * 10}\n>high\n{'N' * 31 + 'A' * 569}\n"
                f">boundary\n{'n' * 30 + 'a' * 570}\n>clear\n{'ACGT' * 150}\n")
            args = parse_args(["prepare", "--fasta", str(root / "input.fa"), "--out", str(root / "output")])
            self.assertEqual(args.max_n_fraction, 0.05)
            with patch("cosmos_abundance.reference.prepare.run_mafft", side_effect=RuntimeError("stop before alignment")) as align:
                with self.assertRaisesRegex(RuntimeError, "stop before alignment"):
                    prepare_reference(args)
                align.assert_called_once()
            kept = parse_fasta_file(root / "output/input.fasta")
            self.assertEqual([record["id"] for record in kept], ["boundary", "clear"])
            report = (root / "output/removed_sequences.tsv").read_text()
            self.assertIn("high\thigh_n_content\t", report)
            self.assertIn("short\tshort_input\t", report)

    def test_strict_n_filter_writes_report_when_all_inputs_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input.fa").write_text(f">one\n{'A' * 599}N\n")
            args = parse_args(["prepare", "--fasta", str(root / "input.fa"), "--out", str(root / "output"),
                               "--max-n-fraction", "0"])
            with patch("cosmos_abundance.reference.prepare.run_mafft") as align:
                with self.assertRaisesRegex(ValueError, "N-fraction"):
                    prepare_reference(args)
                align.assert_not_called()
            self.assertIn("one\thigh_n_content", (root / "output/removed_sequences.tsv").read_text())

    def test_n_fraction_argument_bounds(self):
        from cosmos_abundance.cli import proportion
        import argparse
        self.assertEqual(proportion("0"), 0)
        self.assertEqual(proportion("1"), 1)
        for value in ["-0.1", "1.1", "nan", "inf"]:
            with self.assertRaises(argparse.ArgumentTypeError):
                proportion(value)

    def test_optional_tree_labels_and_legacy_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "species.tsv"
            path.write_text("id\tspecies\ttree_label\na\tSpecies alpha\tShort alpha\nb\tSpecies beta\t\n")
            species, labels = read_species_metadata(path, ["a", "b"])
            self.assertEqual(labels, {"a": "Short alpha", "b": "Species beta | b"})
            self.assertEqual(read_species_map(path, ["a", "b"]), species)
            path.write_text("a\tSpecies alpha\nb\tSpecies beta\n")
            self.assertEqual(read_species_map(path, ["a", "b"]), species)
            path.write_text("tree_label\tspecies\tid\nShort alpha\tSpecies alpha\ta\n")
            self.assertEqual(read_species_metadata(path, ["a"])[1], {"a": "Short alpha"})

    def test_long_tree_labels_fit_without_layout_warnings(self):
        labels = {"a": "Wide W label " * 25, "b": "Another label"}
        self.assertGreater(expected_figure_width(list(labels.values())), expected_figure_width(["short"]))
        with tempfile.TemporaryDirectory() as directory, warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            for ids, distances in [(["a", "b"], np.array([[0, 0.02], [0.02, 0]])),
                                    (["a"], np.zeros((1, 1)))]:
                path = Path(directory) / f"tree_{len(ids)}.svg"
                plot_mash_tree(ids, distances, ["a"], 0.01, path, labels)
                self.assertIn(labels["a"], path.read_text())
            self.assertFalse(any("layout" in str(item.message).lower() for item in caught), caught)

    def test_prepare_reports_filtered_and_clustered_sequences(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input.fa").write_text(f">a\n{'ACGT' * 150}\n>b\n{'ACGT' * 150}\n>short\nACGT\n")
            (root / "species.tsv").write_text("id\tspecies\ttree_label\na\tAlpha\tCustom A\nb\tBeta\tCustom B\nshort\tOther\t\n")
            args = parse_args(["prepare", "--fasta", str(root / "input.fa"), "--species",
                               str(root / "species.tsv"), "--out", str(root / "output"), "--keep-intermediate"])
            def align(source, output, threads):
                shutil.copyfile(source, output)
            with patch("cosmos_abundance.reference.prepare.run_mafft", side_effect=align), \
                 patch("cosmos_abundance.reference.prepare.subprocess.run", return_value=SimpleNamespace(stdout="test", stderr="")), \
                 patch("cosmos_abundance.reference.prepare.plot_mash_tree") as plot:
                prepare_reference(args)
            import csv
            with open(root / "output/removed_sequences.tsv") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows, [{"id": "short", "reason": "short_input", "representative_id": ""},
                                    {"id": "b", "reason": "clustered", "representative_id": "a"}])
            self.assertIn("a\tAlpha + Beta", (root / "output/assignment_mapping.tsv").read_text())
            self.assertEqual(plot.call_args.args[-1]["a"], "Custom A")

    def test_small_alignment_retains_sequence_and_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "msa.fa").write_text(">one\n--ACGTACGT--\n>two\n--ACGTACGT--\n")
            records, removed, bounds = trim_alignment(root / "msa.fa", root, 8, 0, 0.5)
            self.assertEqual([record["sequence"] for record in records], ["ACGTACGT"] * 2)
            self.assertEqual([record["id"] for record in records], ["one", "two"])
            self.assertEqual(removed, [])
            self.assertEqual(bounds["alignment_start_1based"], 3)

    def test_gap_parameter_counts_intervening_columns(self):
        self.assertEqual(find_best_gap_free_region(np.array([1, 1, 0, 1, 1], dtype=bool), 1), (0, 4, 5))
        self.assertEqual(find_best_gap_free_region(np.array([1, 1, 0, 1, 1], dtype=bool), 0), (0, 1, 2))

    def test_minimum_length_is_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "msa.fa").write_text(">one\nACGT\n")
            with self.assertRaisesRegex(ValueError, "No shared region"):
                trim_alignment(root / "msa.fa", root, 5, 0, 0.5)

    def test_duplicate_ids_and_incomplete_mapping_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input.fa").write_text(">same\nACGT\n>same\nACGT\n")
            with self.assertRaisesRegex(ValueError, "Duplicate FASTA"):
                parse_fasta_file(root / "input.fa")
            (root / "species.tsv").write_text("id\tspecies\na\tSpecies alpha\n")
            with self.assertRaisesRegex(ValueError, "Missing species"):
                read_species_map(root / "species.tsv", ["a", "b"])

    def test_complete_linkage_and_deterministic_representative(self):
        distances = np.array([[0, 0.001, 0.2], [0.001, 0, 0.2], [0.2, 0.2, 0]])
        clusters = cluster_sequences_by_distance(["b", "a", "c"], distances, 0.01)
        representatives = identify_best_representatives(clusters, distances)
        self.assertEqual(set(representatives.filter(representatives["is_best_representative"])["seq_id"]), {"a", "c"})


if __name__ == "__main__":
    unittest.main()
