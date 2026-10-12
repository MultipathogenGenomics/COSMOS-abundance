from pathlib import Path
import tempfile
import unittest

import numpy as np

from cosmos_abundance.io import parse_fasta_file, read_species_map
from cosmos_abundance.reference.cluster import cluster_sequences_by_distance, identify_best_representatives
from cosmos_abundance.reference.fasta import find_best_gap_free_region, trim_alignment


class ReferenceTests(unittest.TestCase):
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
