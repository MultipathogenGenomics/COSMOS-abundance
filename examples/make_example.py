"""Write small synthetic references and paired reads; no biological data required."""

import gzip
from pathlib import Path
import random


def make_example(output_dir):
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    sequence = "".join(rng.choices("ACGT", k=650))
    variant = list(sequence)
    for index in rng.sample(range(650), 40):
        variant[index] = rng.choice([base for base in "ACGT" if base != variant[index]])
    variant = "".join(variant)
    (root / "references.fasta").write_text(f">ref_a\n{sequence}\n>ref_b\n{sequence}\n>ref_c\n{variant}\n>short\nACGT\n")
    (root / "species.tsv").write_text("id\tspecies\nref_a\tSpecies alpha\nref_b\tSpecies beta\nref_c\tSpecies gamma\nshort\tSpecies delta\n")
    complement = str.maketrans("ACGT", "TGCA")
    for sample in ["mixture", "unassigned"]:
        with gzip.open(root / f"{sample}_R1.fastq.gz", "wt") as read1, gzip.open(root / f"{sample}_R2.fastq.gz", "wt") as read2:
            for index in range(200):
                target = sequence if index < 150 else variant
                if sample == "unassigned":
                    target = "".join(rng.choices("ACGT", k=650))
                start = rng.randrange(351)
                fragment = target[start:start + 300]
                for handle, read, mate in [(read1, fragment[:100], 1),
                                            (read2, fragment[-100:].translate(complement)[::-1], 2)]:
                    handle.write(f"@{sample}_{index}/{mate}\n{read}\n+\n{'I' * len(read)}\n")
    (root / "samples.tsv").write_text("sample\tread1\tread2\nmixture\tmixture_R1.fastq.gz\tmixture_R2.fastq.gz\nunassigned\tunassigned_R1.fastq.gz\tunassigned_R2.fastq.gz\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="example-data")
    make_example(parser.parse_args().out)
