# COSMOS-abundance

COSMOS-abundance is a command-line workflow for estimating mosquito species composition and relative read abundance from sequencing data containing the mitochondrial COX1 marker.

The workflow was developed for analysis of mixed-species mosquito pools generated using the COSMOS targeted-metagenomics panel described in:

[COSMOS: integrated vector-pathogen metagenomics for mosquito-borne disease surveillance](LINK)

Although developed and benchmarked for mosquitoes, the workflow could potentially be applied to other taxa where a suitable homologous marker and reference set are available. The workflow may also be applied to shotgun metagenomics data, given sufficient marker sequence recovery, but this approach should be carefully tested.

N.B. The filtering and assignment parameters reproduce the analysis described in the manuscript above, and should not be assumed to be optimal for every dataaset. The differential effects of sample-reference divergence and similarity among competing references on KMA and Kallisto read assignment are described in [this preprint](LINK).

## Workflow overview

The workflow has two main stages:

1. Prepare a marker reference set
    - filter sequences based on length and missing sequence
    - trim markers to a shared sequence region
    - cluster references that are too similar to resolve and select a representative
    - draw a distance tree for visualisation

2. Assign reads and estimate relative abundance
    - index the reference set for KMA and Kallisto
    - assign reads using both tools
    - combine results belonging to the same speces/analytical group
    - filter based on coverage, identity and read counts
    - retain hits supported by both tools
    - estimate relative abundance from KMA-assigned reads among retained targets
    - re-estimate using different filters

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
marker/locus, in the same orientation. 

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

Preparation filters short sequences and sequences with too many Ns, aligns with MAFFT, trims alignment ends,
selects a shared region, and clusters by complete linkage. Each cluster retains
the sequence with the smallest mean distance to its members (the medoid), with
ties broken by ID. Names are the sorted unique species joined with ` + `:
`Aedes aegypti + Aedes albopictus`. This means the cluster cannot distinguish
those species; it does not split abundance between them. The separator ` + ` is
reserved and cannot appear inside a supplied species name.

The representative keeps its original ID. `assignment_mapping.tsv` has columns
`target_id` and `assignment_name`, connecting each ID to its combined species name.
Assignment names are not added to FASTA headers. With `--keep-intermediate`,
`clusters.tsv` also records all retained cluster members.
Species removed during filtering/trimming do not contribute to cluster names.

| Option | Default | Meaning |
| --- | --- | --- |
| `--mash-distance-threshold` | 0.01 | Complete-linkage clustering distance |
| `--min-reference-length` | 500 | Minimum ungapped bases before and after trimming |
| `--max-n-fraction` | 0.05 | Maximum proportion of N bases before alignment |
| `--kmer-size` | 31 | K-mer size for distance calculation |
| `--allowed-gap-size` | 3 | Maximum intervening gappy columns joining gap-free stretches |
| `--min-occupancy` | 0.5 | Minimum occupied fraction defining alignment ends |
| `--threads` | 1 | MAFFT threads |

Before alignment, sequences shorter than `--min-reference-length` are removed
first, then sequences with `N count / sequence length > --max-n-fraction`.
The limit is inclusive: exactly 5% Ns passes the default 0.05 threshold.
Lowercase `n` is counted after FASTA normalisation to uppercase. Other IUPAC
ambiguity codes do not count as Ns. Use `0` to reject any Ns or `1` to disable
this N-content filter. The default is 5% and applies to all input references.
Rejections appear as `short_input` or `high_n_content` in `removed_sequences.tsv`;
this initial report is saved even if all inputs are rejected. The threshold is
also recorded in the manifest when `--keep-intermediate` is used. Length includes N bases.

Sequences with no unambiguous A/C/G/T k-mers are excluded. If the shared region
is too short, the most gappy sequence is removed and the bounds are recomputed,
following the original approach. Preparation fails if no suitable reference
remains. On successful runs, `removed_sequences.tsv` lists every input sequence absent
from the final representative FASTA. Filtered sequences have their rejection
reason; non-representative cluster members have reason `clustered` and a
`representative_id` linking them to the retained sequence. A header-only report
means every input sequence was retained as a representative. The manifest and
console summary distinguish filtering from replacement by a representative.

By default, preparation leaves only:

```text
reference-set/
├── prepared_references.fasta
├── assignment_mapping.tsv
├── removed_sequences.tsv
└── mash_tree.png
```

The PNG labels representatives in bold and marks the clustering threshold for
multi-sequence inputs. Figure width is estimated from the longest displayed label.
SVG output is no longer generated.

Add `--keep-intermediate` to retain input/alignment FASTAs, the MAFFT log,
`clusters.tsv`, `mash_distances.tsv` and `manifest.json` as well. Intermediate
files are removed only after preparation completes successfully.

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

Normal execution shows a progress bar based on completed workflow jobs; it is
not an estimate of elapsed time. Snakemake and tool output go to `workflow.log`.
On failure, error details are printed and the log and intermediates are retained.
On success, the log is retained only with `--keep-intermediate`.
`--dry-run` still prints the proposed jobs and retains its configuration.

Use `--keep-intermediate` when you want to resume a completed run or change
filters without remapping. Without it, successful runs retain only the three
final files below; subsequent runs need a new output directory. Failed runs can
be resumed by repeating the same command. `.snakemake` is removed after success,
even with `--keep-intermediate`. This discards Snakemake's scheduling metadata;
the saved inputs and parameters still allow unchanged runs and filter changes.

Use a new directory if references, assignment names, FASTQ paths/sizes/modification
times, bootstrap count, seed or package version change. FASTQ contents are not
checksummed. Do not modify run inputs in place.

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

Successful runs leave only:

```text
results/
└── relative_abundance/
    ├── parameters.json
    ├── sample_read_assignment_summary.tsv
    └── abundance_within_filter_pass_reads.tsv
```

`abundance_within_filter_pass_reads.tsv` contains counts divided by all passing
counts. The summary records mapping totals and zero-passing-read flags.
`parameters.json` records filters, count conventions, input metadata and workflow
settings so that these remain available after cleanup.

With `--keep-intermediate`, the run also retains read-assignment outputs and
indexes, configuration and input files, `workflow.log`, and these additional tables:

- `filter_pass_read_counts.tsv`: passing KMA counts by sample and assignment name.
- `abundance_within_all_assigned_reads.tsv`: includes `non_confidently_assigned`;
  the denominator follows the original paired-end convention of KMA assigned
  fragments multiplied by two. Inconsistent negative residuals raise an error.
- `kallisto_kma_outputs_joined.tsv`: reference-level metrics and assignment names.

All input samples are retained. Zero-assignment samples have zeros rather than
NaNs; if no target passes anywhere, the passing-count table has only its sample
column. These are relative read abundances, not calibrated organism counts.

For a run created with `--keep-intermediate`, recalculate with different filters:

```bash
cosmos-abundance estimate --run-dir results/ --min-coverage 95 --keep-intermediate
```

This replaces the relative abundance tables and their parameters. Omit
`--keep-intermediate` on `estimate` to discard intermediates after recalculation.
Subsequent `run` calls
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
the direct-FASTA interface, numeric sample names, an entirely unassigned batch,
compact outputs, retained intermediates and quiet progress reporting.

