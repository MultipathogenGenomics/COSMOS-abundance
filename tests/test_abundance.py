import unittest

import polars as pl

from cosmos_abundance.abundance import convert_counts_to_abundance, make_predictor_table


class AbundanceTests(unittest.TestCase):
    def test_fused_names_are_preserved_and_counts_aggregate(self):
        data = pl.DataFrame({"sample": ["s", "s", "s"], "target_id": ["a", "b", "c"],
                             "kma:readCount": [10, 20, 30]})
        names = {"a": "Alpha + Beta", "b": "Alpha + Beta", "c": "Gamma"}
        counts = make_predictor_table(data, "kma:readCount", names)
        self.assertEqual(counts["Alpha + Beta"].item(), 30)
        self.assertEqual(convert_counts_to_abundance(counts)["Alpha + Beta"].item(), 0.5)

    def test_no_passing_targets_returns_valid_empty_table(self):
        data = pl.DataFrame({"sample": ["s"], "target_id": ["a"], "kma:readCount": [0]})
        counts = make_predictor_table(data, "kma:readCount", {"a": "Alpha"},
                                      filter_threshold_dict={"kma:readCount": 5})
        self.assertEqual(counts.columns, ["sample"])
        self.assertEqual(counts.height, 0)

    def test_zero_counts_do_not_create_nan(self):
        result = convert_counts_to_abundance(pl.DataFrame({"sample": ["s"], "Alpha": [0]}))
        self.assertEqual(result["Alpha"].item(), 0)


if __name__ == "__main__":
    unittest.main()
