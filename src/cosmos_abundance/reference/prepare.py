from pathlib import Path
import subprocess

from .. import __version__
from ..io import file_sha256, parse_fasta_file, read_species_metadata, write_fasta, write_json, write_tsv
from .cluster import cluster_sequences_by_distance, identify_best_representatives, select_representative_records
from .fasta import run_mafft, trim_alignment
from .mash import compute_mash_distance_matrix, sequence_to_kmers, write_distance_matrix
from .plot import plot_mash_tree


def prepare_reference(args) -> None:
    records = sorted(parse_fasta_file(args.fasta), key=lambda record: record["id"])
    species_map, labels = read_species_metadata(args.species, [record["id"] for record in records])
    output_dir = Path(args.out)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Preparation output directory must be empty; use a new directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    kept, removed = [], []
    for record in records:
        if len(record["sequence"]) < args.min_reference_length:
            removed.append({"id": record["id"], "reason": "short_input"})
        elif not sequence_to_kmers(record["sequence"], args.kmer_size):
            removed.append({"id": record["id"], "reason": "no_valid_kmers"})
        else:
            kept.append(record)
    if not kept:
        raise ValueError("No input references meet the length and k-mer requirements")
    write_fasta(kept, output_dir / "input.fasta")
    print(f"Aligning {len(kept)} references with MAFFT...", flush=True)
    run_mafft(output_dir / "input.fasta", output_dir / "alignment.fasta", args.threads)
    trimmed, trim_removed, trim_bounds = trim_alignment(
        output_dir / "alignment.fasta", output_dir, args.min_reference_length,
        args.allowed_gap_size, args.min_occupancy)
    removed.extend(trim_removed)
    kept = []
    for record in trimmed:
        if sequence_to_kmers(record["sequence"], args.kmer_size):
            kept.append(record)
        else:
            removed.append({"id": record["id"], "reason": "no_valid_kmers_after_trimming"})
    if not kept:
        raise ValueError("No trimmed references contain valid k-mers")
    ids, distances = compute_mash_distance_matrix(kept, args.kmer_size)
    write_distance_matrix(ids, distances, output_dir / "mash_distances.tsv")
    clusters = cluster_sequences_by_distance(ids, distances, args.mash_distance_threshold)
    clusters = identify_best_representatives(clusters, distances)
    representatives = select_representative_records(kept, clusters)
    filtered_count = len(removed)
    cluster_rows, assignment_rows = [], []
    for cluster_id in sorted(clusters["cluster_id"].unique().to_list()):
        rows = [row for row in clusters.to_dicts() if row["cluster_id"] == cluster_id]
        representative_id = next(row["seq_id"] for row in rows if row["is_best_representative"])
        name = " + ".join(sorted({species_map[row["seq_id"]] for row in rows}))
        assignment_rows.append({"target_id": representative_id, "assignment_name": name})
        for row in rows:
            if not row["is_best_representative"]:
                removed.append({"id": row["seq_id"], "reason": "clustered", "representative_id": representative_id})
            cluster_rows.append({"id": row["seq_id"], "species": species_map[row["seq_id"]],
                                 "cluster_id": cluster_id, "representative_id": representative_id,
                                 "assignment_name": name, "is_representative": row["is_best_representative"]})
    write_fasta(representatives, output_dir / "references.fasta")
    write_tsv(output_dir / "assignments.tsv", assignment_rows, ["target_id", "assignment_name"])
    write_tsv(output_dir / "clusters.tsv", cluster_rows,
              ["id", "species", "cluster_id", "representative_id", "assignment_name", "is_representative"])
    write_tsv(output_dir / "removed_sequences.tsv", removed, ["id", "reason", "representative_id"])
    for extension in ["png", "svg"]:
        plot_mash_tree(ids, distances, [record["id"] for record in representatives],
                       args.mash_distance_threshold, output_dir / f"mash_tree.{extension}", labels)
    mafft_version = subprocess.run(["mafft", "--version"], capture_output=True, text=True, check=True)
    write_json(output_dir / "manifest.json", {
        "version": __version__, "mafft_version": (mafft_version.stdout + mafft_version.stderr).strip(),
        "input_fasta": str(Path(args.fasta).resolve()), "input_sha256": file_sha256(args.fasta),
        "species_file": str(Path(args.species).resolve()) if args.species else None,
        "species_sha256": file_sha256(args.species) if args.species else None,
        "reference_sha256": file_sha256(output_dir / "references.fasta"),
        "assignments_sha256": file_sha256(output_dir / "assignments.tsv"),
        "parameters": {key: value for key, value in vars(args).items() if key != "command"},
        "input_sequences": len(records), "retained_sequences": len(kept),
        "representatives": len(representatives), "removed_sequences": len(removed),
        "filtered_sequences": filtered_count, "clustered_sequences": len(removed) - filtered_count, **trim_bounds,
    })
    print(f"Excluded {filtered_count} sequences by filtering; replaced {len(removed) - filtered_count} with cluster representatives.")
    print(f"Wrote {len(representatives)} representatives and Mash-distance trees to {output_dir}")
