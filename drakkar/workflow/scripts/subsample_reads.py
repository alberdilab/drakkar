"""Stream deterministic paired subsampling through seqkit."""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from drakkar.subsample import subsample_fraction


def seqkit_stats(path):
    result = subprocess.run(
        ["seqkit", "stats", "-T", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one seqkit stats row for {path}")
    row = rows[0]
    return int(row["num_seqs"].replace(",", "")), int(row["sum_len"].replace(",", ""))


def stage_copy(source, destination):
    """Link unchanged input where possible, otherwise copy it to a staged path."""
    try:
        os.link(source, destination)
    except OSError:
        shutil.copyfile(source, destination)


def run(args):
    r1 = Path(args.r1)
    r2 = Path(args.r2)
    out_r1 = Path(args.out_r1)
    out_r2 = Path(args.out_r2)
    stats = Path(args.stats)
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    stats.parent.mkdir(parents=True, exist_ok=True)

    before_r1, bases_r1 = seqkit_stats(r1)
    before_r2, bases_r2 = seqkit_stats(r2)
    if before_r1 != before_r2:
        raise ValueError(f"{args.sample}: R1 and R2 have different record counts")
    fraction = subsample_fraction(args.target, before_r1, bases_r1 + bases_r2)

    staged_r1 = Path(str(out_r1) + ".tmp.fq.gz")
    staged_r2 = Path(str(out_r2) + ".tmp.fq.gz")
    staged_stats = Path(str(stats) + ".tmp")
    for source, staged in ((r1, staged_r1), (r2, staged_r2)):
        if staged.exists():
            staged.unlink()
        if fraction == 1:
            stage_copy(source, staged)
        else:
            subprocess.run(
                ["seqkit", "sample", "-p", format(fraction, ".17g"),
                 "-s", str(args.seed), "-2", "-o", str(staged), str(source)],
                check=True,
            )

    after_r1, _ = seqkit_stats(staged_r1)
    after_r2, _ = seqkit_stats(staged_r2)
    if after_r1 != after_r2:
        raise ValueError(f"{args.sample}: subsampling produced unequal mate counts")

    with staged_stats.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "reads_before", "reads_subsampled", "pairs_before", "pairs_after"])
        writer.writerow([args.sample, 2 * before_r1, 2 * after_r1, before_r1, after_r1])

    os.replace(staged_r1, out_r1)
    os.replace(staged_r2, out_r2)
    os.replace(staged_stats, stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--out-r1", required=True)
    parser.add_argument("--out-r2", required=True)
    parser.add_argument("--stats", required=True)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
