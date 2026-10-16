import csv
import hashlib
import json
import shutil
from pathlib import Path

from Bio import SeqIO


def parse_fasta_file(fasta_path: str) -> list[dict]:
    records = []
    seen = set()
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in seen:
                raise ValueError(f"Duplicate FASTA ID: {record.id}")
            sequence = str(record.seq).upper()
            if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
                raise ValueError(f"Expected ungapped DNA sequence for {record.id}")
            seen.add(record.id)
            records.append({"id": record.id, "sequence": sequence})
    if not records:
        raise ValueError(f"No sequences found in {fasta_path}")
    return records


def write_fasta(records: list[dict], output_fasta: str) -> None:
    with open(output_fasta, "w") as handle:
        for record in records:
            handle.write(f">{record['id']}\n{record['sequence']}\n")


def read_species_metadata(path: str | None, ids: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    if path is None:
        return dict(zip(ids, ids)), {seq_id: f"{seq_id} | {seq_id}" for seq_id in ids}
    species_map, tree_labels = {}, {}
    with open(path) as handle:
        rows = [[value.strip() for value in row] for row in csv.reader(handle, delimiter="\t") if row]
    columns = ["id", "species"]
    if rows and {"id", "species"}.issubset(rows[0]):
        columns = rows.pop(0)
        if len(columns) != len(set(columns)) or set(columns) - {"id", "species", "tree_label"}:
            raise ValueError("Species file header must contain id, species and optionally tree_label")
    for row in rows:
        if len(row) != len(columns):
            raise ValueError("Species file rows must match the tab-separated id, species and optional tree_label columns")
        values = dict(zip(columns, row))
        seq_id, species = values["id"], values["species"]
        if not seq_id or not species:
            raise ValueError("Species file requires non-empty id and species values")
        if seq_id in species_map:
            raise ValueError(f"Duplicate ID in species file: {seq_id}")
        if " + " in species:
            raise ValueError("Species names cannot contain the reserved cluster separator ' + '")
        species_map[seq_id] = species
        tree_labels[seq_id] = values.get("tree_label") or f"{species} | {seq_id}"
    missing = set(ids) - set(species_map)
    if missing:
        raise ValueError(f"Missing species for FASTA IDs: {', '.join(sorted(missing))}")
    return ({seq_id: species_map[seq_id] for seq_id in ids},
            {seq_id: tree_labels[seq_id] for seq_id in ids})


def read_species_map(path: str | None, ids: list[str]) -> dict[str, str]:
    species_map, _ = read_species_metadata(path, ids)
    return species_map


def write_tsv(path: str, rows: list[dict], columns: list[str]) -> None:
    with open(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: str, data: dict) -> None:
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def remove_intermediate_files(root: Path, names: list[str]) -> None:
    # Remove only known workflow outputs, never unrelated files in the directory.
    for name in names:
        path = root / name
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def clean_run_outputs(root: Path, keep_intermediate: bool) -> None:
    remove_intermediate_files(root, [".snakemake"])
    if keep_intermediate:
        return
    remove_intermediate_files(root, ["read-assignment", "references.fasta", "assignments.tsv",
                                    "samples.tsv", "inputs.json", "config.json", "workflow.log"])
    remove_intermediate_files(root / "relative_abundance", ["filter_pass_read_counts.tsv",
        "abundance_within_all_assigned_reads.tsv", "kallisto_kma_outputs_joined.tsv"])
