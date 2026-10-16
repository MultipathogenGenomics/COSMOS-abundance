import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections import deque

import polars as pl

from . import __version__
from .io import file_sha256, parse_fasta_file, read_species_map, write_json, write_tsv, clean_run_outputs


def execute_workflow(command: list[str], log_path: Path) -> None:
    print("Workflow [--------------------] starting", end="", flush=True)
    tail = deque(maxlen=25)
    with open(log_path, "w") as log, subprocess.Popen(
        command + ["--nocolor"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    ) as process:
        for line in process.stdout:
            log.write(line)
            log.flush()
            tail.append(line)
            match = re.search(r"(\d+) of (\d+) steps", line)
            if match:
                done, total = map(int, match.groups())
                if total:
                    filled = min(20, done * 20 // total)
                    print(f"\rWorkflow [{'#' * filled}{'-' * (20 - filled)}] {done}/{total} ({done / total:.0%})", end="", flush=True)
        status = process.wait()
    if status:
        print("\nWorkflow failed.", flush=True)
        print("".join(tail), file=sys.stderr)
        raise ValueError(f"Snakemake failed (exit {status}); full log: {log_path}")
    print("\rWorkflow [####################] 100% complete", flush=True)


def read_samples(path: str) -> dict:
    table = pl.read_csv(path, separator="\t", infer_schema=False)
    if not {"sample", "read1", "read2"}.issubset(table.columns) or table.is_empty():
        raise ValueError("Sample sheet must have sample, read1 and read2 columns and at least one sample")
    samples = {}
    for row in table.to_dicts():
        sample = row["sample"]
        if not sample or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"Invalid sample name: {sample!r}; use letters, digits, '.', '_' or '-'")
        if sample in samples:
            raise ValueError(f"Duplicate sample name: {sample}")
        paths = {}
        for key in ["read1", "read2"]:
            if not row[key]:
                raise ValueError(f"Missing {key} for {sample}")
            read = Path(row[key])
            if not read.is_absolute():
                read = Path(path).resolve().parent / read
            read = read.resolve()
            if not read.is_file():
                raise ValueError(f"FASTQ does not exist: {read}")
            paths[key] = str(read)
        if paths["read1"] == paths["read2"]:
            raise ValueError(f"R1 and R2 must be different files for {sample}")
        samples[sample] = paths
    return samples


def run_pipeline(args) -> None:
    samples = read_samples(args.samples)
    if args.reference:
        if args.species:
            raise ValueError("--species is only used with --fasta; prepared references already contain names")
        reference = Path(args.reference).resolve()
        fasta = reference / "prepared_references.fasta" if reference.is_dir() else reference
        if fasta.is_file():
            mapping_file = fasta.parent / "assignment_mapping.tsv"
        else:
            # Support reference directories produced by earlier versions.
            fasta = reference / "references.fasta"
            mapping_file = reference / "assignments.tsv"
        manifest_file = fasta.parent / "manifest.json"
        if manifest_file.exists():
            manifest = json.loads(manifest_file.read_text())
            if file_sha256(fasta) != manifest["reference_sha256"] or file_sha256(mapping_file) != manifest["assignments_sha256"]:
                raise ValueError("Prepared reference files have changed since preparation")
        names = pl.read_csv(mapping_file, separator="\t", infer_schema=False)
    else:
        fasta = Path(args.fasta).resolve()
        records = parse_fasta_file(fasta)
        mapping = read_species_map(args.species, [record["id"] for record in records])
        names = pl.DataFrame({"target_id": list(mapping), "assignment_name": list(mapping.values())})
    records = parse_fasta_file(fasta)
    if set(names["target_id"]) != {record["id"] for record in records} or names["target_id"].n_unique() != names.height:
        raise ValueError("Reference IDs and assignment names must match one-to-one")
    if any(name in {"sample", "non_confidently_assigned", "total_reads", "total_filter_pass_reads", "kma:n_assigned"}
           for name in names["assignment_name"]):
        raise ValueError("Assignment names conflict with reserved output column names")
    if any(len(record["sequence"]) < 31 for record in records):
        raise ValueError("Read assignment requires reference sequences of at least 31 bases")
    for program in ["kma", "kma_index", "kallisto"]:
        if not shutil.which(program):
            raise ValueError(f"Required executable is not on PATH: {program}")
    root = Path(args.out).resolve()
    if root == fasta.parent:
        raise ValueError("Use separate reference and run output directories")
    inputs = {"reference_sha256": file_sha256(fasta), "assignments": names.to_dicts(), "samples": samples,
              "reads": {path: {"size": Path(path).stat().st_size, "mtime_ns": Path(path).stat().st_mtime_ns}
                        for sample in samples.values() for path in sample.values()},
              "bootstraps": args.bootstraps, "seed": args.seed, "version": __version__}
    if (root / "inputs.json").exists():
        if json.loads((root / "inputs.json").read_text()) != inputs:
            raise ValueError("Run inputs have changed; use a new output directory")
    elif root.exists() and any(root.iterdir()):
        raise ValueError("Output directory has no retained run inputs; use a new directory. Use --keep-intermediate to allow future resumes")
    else:
        root.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(fasta, root / "references.fasta")
        names.write_csv(root / "assignments.tsv", separator="\t")
        write_tsv(root / "samples.tsv", [{"sample": name, **paths} for name, paths in samples.items()],
                  ["sample", "read1", "read2"])
        write_json(root / "inputs.json", inputs)
    if file_sha256(root / "references.fasta") != inputs["reference_sha256"]:
        raise ValueError("The run's reference FASTA has been modified")
    if pl.read_csv(root / "assignments.tsv", separator="\t", infer_schema=False).to_dicts() != names.to_dicts():
        raise ValueError("The run's assignment names have been modified")
    config = {"fasta": str(root / "references.fasta"), "samples": samples,
              **{key: getattr(args, key) for key in ["threads", "bootstraps", "seed", "min_coverage",
                  "min_identity", "min_kma_count", "min_kallisto_count", "allow_single_tool"]}}
    write_json(root / "config.json", config)
    workflow = Path(__file__).parent / "workflow/Snakefile"
    command = [sys.executable, "-m", "snakemake", "--snakefile", str(workflow),
               "--directory", str(root), "--configfile", str(root / "config.json"),
               "--cores", str(args.cores), "--rerun-incomplete"]
    if args.keep_intermediate:
        command.append("--notemp")
    if args.dry_run:
        command.append("--dry-run")
    parameters_file = root / "relative_abundance/parameters.json"
    if parameters_file.exists():
        expected_filters = {"kma:template_coverage": args.min_coverage,
                            "kma:identity_across_covered_length": args.min_identity,
                            "kma:readCount": args.min_kma_count,
                            "kallisto:est_counts_mean": args.min_kallisto_count}
        if not args.allow_single_tool:
            expected_filters["identified_by_both_tools"] = True
        if json.loads(parameters_file.read_text())["filters"] != expected_filters:
            command.extend(["--forcerun", "estimate_abundance"])
    if args.dry_run:
        subprocess.run(command, check=True)
    else:
        execute_workflow(command, root / "workflow.log")
    if not args.dry_run:
        clean_run_outputs(root, args.keep_intermediate)
