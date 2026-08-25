#!/usr/bin/env python

import argparse
import gzip
import os
import resource
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from time import perf_counter

import pandas as pd


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract per-cell pairs files based on barcode metadata."
    )
    parser.add_argument("--indir", dest="iput", required=True, help="Input directory")
    parser.add_argument(
        "--cluster",
        required=True,
        help="Four-column metadata: barcode, library, sample.pairs.gz, cluster",
    )
    parser.add_argument(
        "--max",
        dest="maxf",
        type=positive_int,
        default=1000,
        help="Maximum number of output files open at once",
    )
    parser.add_argument("--outdir", dest="oput", required=True, help="Output directory")
    parser.add_argument(
        "--threads",
        type=positive_int,
        default=4,
        help="Number of input files to process concurrently",
    )
    return parser.parse_args()


class OutputWriter:
    """Thread-safe, bounded LRU cache of append-only output handles."""

    def __init__(self, output_paths, max_open_files):
        self.output_paths = output_paths
        self.max_open_files = max_open_files
        self.handles = OrderedDict()
        self.lock = Lock()

    def write(self, cluster_id, line):
        with self.lock:
            handle = self.handles.pop(cluster_id, None)
            if handle is None:
                if len(self.handles) >= self.max_open_files:
                    _, oldest_handle = self.handles.popitem(last=False)
                    oldest_handle.close()
                handle = open(self.output_paths[cluster_id], "a", encoding="utf-8")
            self.handles[cluster_id] = handle
            handle.write(line)

    def close(self):
        with self.lock:
            while self.handles:
                _, handle = self.handles.popitem(last=False)
                handle.close()


def open_pairs(filename, mode="rt"):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode, encoding="utf-8")
    return open(filename, mode, encoding="utf-8")


def safe_open_file_limit(requested):
    soft_limit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    if soft_limit == resource.RLIM_INFINITY:
        return requested
    return max(1, min(requested, max(1, soft_limit - 32)))


def load_metadata(cluster_file):
    metadata = pd.read_csv(cluster_file, sep="\t", dtype=str, keep_default_na=False)
    required_columns = {"barcode", "library", "sample", "cluster"}
    missing = required_columns.difference(metadata.columns)
    if missing:
        raise ValueError(
            "Cluster metadata is missing columns: " + ", ".join(sorted(missing))
        )
    if metadata.empty:
        raise ValueError("Cluster metadata contains no cells")
    if (metadata[sorted(required_columns)] == "").to_numpy().any():
        raise ValueError("Cluster metadata contains empty required values")

    metadata["uniq_barcode"] = metadata["library"] + ":" + metadata["barcode"]
    conflicts = metadata.groupby("uniq_barcode")["cluster"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        raise ValueError(
            "A library/barcode pair maps to multiple clusters: "
            + ", ".join(conflicts.index[:5])
        )
    return metadata


def read_header(sample_path):
    header_lines = []
    with open_pairs(sample_path) as infile:
        for line in infile:
            if not line.startswith("#"):
                break
            header_lines.append(line)
    return "".join(header_lines)


def split_contacts_worker(barcode_to_cluster, sample, library, input_dir, writer):
    sample_path = os.path.join(input_dir, sample)
    matched = 0
    malformed = 0
    cross_barcode = 0

    with open_pairs(sample_path) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                malformed += 1
                continue
            if fields[-1] != fields[-2]:
                cross_barcode += 1
                continue

            cell_id = f"{library}:{fields[-1]}"
            cluster_id = barcode_to_cluster.get(cell_id)
            if cluster_id is None:
                continue

            fields[-2:] = [cell_id, cell_id]
            writer.write(cluster_id, "\t".join(fields) + "\n")
            matched += 1

    return sample, matched, malformed, cross_barcode


def split_contacts(input_dir, output_dir, cluster_file, max_open_files, num_threads):
    metadata = load_metadata(cluster_file)
    barcode_to_cluster = (
        metadata.drop_duplicates("uniq_barcode")
        .set_index("uniq_barcode")["cluster"]
        .to_dict()
    )

    os.makedirs(output_dir, exist_ok=True)
    unique_clusters = metadata["cluster"].drop_duplicates().tolist()
    output_paths = {
        cluster_id: os.path.join(output_dir, f"{cluster_id}.pairs")
        for cluster_id in unique_clusters
    }

    samples = metadata["sample"].drop_duplicates().tolist()
    first_sample_path = os.path.join(input_dir, samples[0])
    header = read_header(first_sample_path)
    for output_path in output_paths.values():
        with open(output_path, "w", encoding="utf-8") as outfile:
            outfile.write(header)

    sample_libraries = {}
    for sample in samples:
        libraries = metadata.loc[metadata["sample"] == sample, "library"].unique()
        if len(libraries) != 1:
            raise ValueError(
                f"Input pairs file {sample!r} is assigned to multiple libraries: "
                + ", ".join(libraries)
            )
        sample_libraries[sample] = libraries[0]

    handle_limit = safe_open_file_limit(max_open_files)
    writer = OutputWriter(output_paths, handle_limit)
    worker_count = min(num_threads, len(samples))
    try:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [
                executor.submit(
                    split_contacts_worker,
                    barcode_to_cluster,
                    sample,
                    sample_libraries[sample],
                    input_dir,
                    writer,
                )
                for sample in samples
            ]
            for future in as_completed(futures):
                sample, matched, malformed, cross_barcode = future.result()
                print(
                    f"{sample}: wrote {matched} pairs; skipped "
                    f"{cross_barcode} cross-barcode and {malformed} malformed records"
                )
    finally:
        writer.close()


def main():
    args = parse_args()
    start_time = perf_counter()
    print("Splitting contacts...")
    split_contacts(args.iput, args.oput, args.cluster, args.maxf, args.threads)
    print(f"Used (secs): {perf_counter() - start_time:.2f}")


if __name__ == "__main__":
    main()
