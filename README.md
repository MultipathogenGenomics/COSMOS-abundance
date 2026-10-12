# COSMOS-abundance

Prepare a representative marker reference set and estimate read-based relative
abundance using KMA and kallisto. This package extracts the reference-set,
read-assignment and estimate-abundance methods from COSMOS-analysis.

## Installation

From this repository:

```bash
conda env create -f environment.yml
conda activate cosmos-abundance
python -m pip install --no-deps --no-build-isolation .
cosmos-abundance --help
```

MAFFT, KMA (including `kma_index`), and kallisto must be on PATH. The Conda
environment includes these programs. The Python package alone does not install
them. This is not yet a published Bioconda package. Linux and macOS are the
intended platforms; the first local smoke test used macOS, MAFFT 7.526,
KMA 1.6.8 and kallisto 0.51.1.

## Prepare references

Supply one ungapped DNA FASTA containing homologous sequences from the same
marker/locus, in the same orientation. FASTA IDs are the first whitespace-delimited
token in each header and must be unique.

An optional tab-separated species file maps those IDs to names:

```text
id	species	tree_label
reference_1	Aedes aegypti	Ae. aegypti (reference 1)
reference_2	Aedes albopictus	Ae. albopictus (reference 2)
reference_3	Culex quinquefasciatus	Cx. quinquefasciatus
```

The `tree_label` column is optional and changes only the displayed tree label,
not species names, clustering or abundance assignments. Blank labels fall back
to `species | FASTA ID`. Include a header when using `tree_label`; the original
two-column `id` / `species` format still accepts files with or without a header.
Every FASTA ID must have a mapping when a species file is supplied. Without a
species file, sequence IDs are used as names.

```bash
cosmos-abundance prepare \
    --fasta references.fasta \
    --species species.tsv \
    --out reference-set/ \
    --mash-distance-threshold 0.01 \
    --min-reference-length 500 \
    --threads 4
```

Preparation filters short sequences, aligns with MAFFT, trims alignment ends,
selects a shared region, and clusters by complete linkage. Each cluster retains
the sequence with the smallest mean distance to its members (the medoid), with
ties broken by ID. Names are the sorted unique species joined with ` + `:
`Aedes aegypti + Aedes albopictus`. This means the cluster cannot distinguish
those species; it does not split abundance between them. The separator ` + ` is
reserved and cannot appear inside a supplied species name.

The representative keeps its original ID. `assignments.tsv` connects that ID to
the combined name, and `clusters.tsv` records all retained cluster members.
Species removed during filtering/trimming do not contribute to cluster names.

| Option | Default | Meaning |
| --- | --- | --- |
| `--mash-distance-threshold` | 0.01 | Complete-linkage clustering distance |
| `--min-reference-length` | 500 | Minimum ungapped bases before and after trimming |
| `--kmer-size` | 31 | K-mer size for distance calculation |
| `--allowed-gap-size` | 3 | Maximum intervening gappy columns joining gap-free stretches |
| `--min-occupancy` | 0.5 | Minimum occupied fraction defining alignment ends |
| `--threads` | 1 | MAFFT threads |

Sequences with no unambiguous A/C/G/T k-mers are excluded. If the shared region
is too short, the most gappy sequence is removed and the bounds are recomputed,
following the original approach. Preparation fails if no suitable reference
remains. On successful runs, `removed_sequences.tsv` lists every input sequence absent
from the final representative FASTA. Filtered sequences have their rejection
reason; non-representative cluster members have reason `clustered` and a
`representative_id` linking them to the retained sequence. A header-only report
means every input sequence was retained as a representative. The manifest and
console summary distinguish filtering from replacement by a representative.

Outputs include:

- `references.fasta`: representative sequences for read assignment.
- `assignments.tsv`, `clusters.tsv`, `removed_sequences.tsv`: names and provenance.
- `mash_tree.png` and `mash_tree.svg`: labelled dendrograms with representative
  labels in bold and the clustering threshold marked for multi-sequence inputs.
  Figure width is estimated from the longest displayed label, leaving room for
  the tree before applying tight layout.
- `alignment.fasta`, `alignment.trimmed.fasta`, `references.trimmed.fasta`:
  intermediate alignments/sequences for inspection.
- `mash_distances.tsv`, `manifest.json`: distances, parameters, checksums and summary.

The tree is a complete-linkage distance dendrogram, not an inferred phylogeny.
Distances use full, non-canonical k-mer sets and the Mash transformation of
Jaccard similarity, preserving the COSMOS-analysis calculation. No Mash executable
or sketches are used. All-pairs distances require quadratic memory, so this is
intended for marker reference panels rather than enormous reference databases.

Use a new or empty output directory for preparation.

## Assign reads and estimate abundance

Provide paired-end FASTQ files (plain or gzip-compressed) in a tab-separated sheet:

```text
sample	read1	read2
sample_A	reads/A_R1.fastq.gz	reads/A_R2.fastq.gz
sample_B	reads/B_R1.fastq.gz	reads/B_R2.fastq.gz
```

Relative FASTQ paths are resolved against the sample sheet directory. Sample
names must start with a letter or digit and contain only letters, digits, `.`,
`_` or `-`.

```bash
cosmos-abundance run \
    --reference reference-set/ \
    --samples samples.tsv \
    --out results/ \
    --cores 8 --threads 2
```

Alternatively, bypass reference preparation entirely:

```bash
cosmos-abundance run \
    --fasta curated-references.fasta \
    --species species.tsv \
    --samples samples.tsv \
    --out results/ --cores 8
```

This validates but does not trim or cluster the supplied FASTA. The species file
is optional; without it, results use reference IDs. Assignment requires sequences
of at least 31 bp. Prepared references already carry their species mapping and
do not accept a replacement `--species` file.

Snakemake runs indexing, KMA, kallisto, collation and abundance estimation.
`--cores` sets the total CPU budget; `--threads` caps threads per assignment job.
Kallisto uses 100 bootstraps and seed 42 by default (`--bootstraps`, `--seed`).
`--dry-run` writes the run configuration and shows the pending jobs.

Repeat the same command to resume. Use a new directory if references, assignment
names, FASTQ paths/sizes/modification times, bootstrap count, seed or package
version change. FASTQ contents are not checksummed. Do not modify run inputs in
place. Changing filters through `run` reruns the abundance stage without remapping.

## Abundance settings and results

Defaults match COSMOS-analysis:

| Option | Default |
| --- | --- |
| `--min-coverage` | 90% KMA template coverage |
| `--min-identity` | 97.5% identity across covered length |
| `--min-kma-count` | 5 reads |
| `--min-kallisto-count` | 5 mean estimated counts across bootstraps |

Both tools must identify a target by default. `--allow-single-tool` disables
that requirement but leaves all numerical filters active. KMA quality filters
and counts still determine the reported abundance.

Counts for references with the same assignment name are summed **before**
filtering. Coverage and identity use the maximum across those references,
preserving the original method. Different references may supply those maxima.
Mixed-species names remain distinct from individual species names.

Final tables are in `results/abundance/`:

- `filter_pass_read_counts.tsv`: passing KMA counts by sample and assignment name.
- `abundance_within_filter_pass_reads.tsv`: counts divided by all passing counts.
- `abundance_within_all_assigned_reads.tsv`: includes `non_confidently_assigned`;
  the denominator follows the original paired-end convention of KMA assigned
  fragments multiplied by two. Inconsistent negative residuals raise an error.
- `sample_read_assignment_summary.tsv`: mapping totals and zero-passing-read flags.
- `kallisto_kma_outputs_joined.tsv`: reference-level metrics and assignment names.
- `parameters.json`: filters and count conventions.

All input samples are retained. Zero-assignment samples have zeros rather than
NaNs; if no target passes anywhere, the passing-count table has only its sample
column. These are relative read abundances, not calibrated organism counts.

Recalculate with different filters from existing assignment results:

```bash
cosmos-abundance estimate --run-dir results/ --min-coverage 95
```

This replaces the abundance tables and their parameters. Subsequent `run` calls
should use the same filter options if those results are to be retained.

## Synthetic example and tests

```bash
python examples/make_example.py --out example-data
cosmos-abundance prepare --fasta example-data/references.fasta \
    --species example-data/species.tsv --out example-reference
cosmos-abundance run --reference example-reference \
    --samples example-data/samples.tsv --out example-results --cores 2 --bootstraps 3
python -m unittest discover -s tests -v
COSMOS_INTEGRATION=1 python -m unittest discover -s tests -v
```

The example includes an identical pair labelled as different species, a divergent
reference, a short reference to remove, a mixed sample and an unassigned sample.
The three bootstraps are for this quick example only.
The integration test additionally checks restart behaviour, changed filters,
the direct-FASTA interface, numeric sample names and an entirely unassigned batch.

