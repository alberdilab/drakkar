import json
import os
from pathlib import Path

import pandas as pd

from drakkar.cli_context import ERROR, RESET
from drakkar.output import print

def normalize_genome_name(name):
    if not name:
        return ""
    base = os.path.basename(str(name).strip())
    if base.endswith(".gz"):
        base = base[:-3]
    for ext in (".fa", ".fna", ".fasta"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base

def load_bins_map(output_dir):
    bins_path = Path(output_dir) / "data" / "bins_to_files.json"
    if not bins_path.exists():
        return {}
    with open(bins_path, "r") as f:
        return json.load(f)

def validate_and_write_quality_file(quality_path, output_dir):
    if not quality_path:
        return False
    if not os.path.isfile(quality_path):
        print(f"{ERROR}ERROR:{RESET} Quality file not found: {quality_path}")
        return False

    df = pd.read_csv(quality_path, sep=None, engine="python", encoding="utf-8-sig")
    col_map = {c: str(c).strip().lstrip("\ufeff").lower() for c in df.columns}
    df.rename(columns=col_map, inplace=True)
    required = {"genome", "completeness", "contamination"}
    if not required.issubset(set(df.columns)):
        print(f"{ERROR}ERROR:{RESET} Quality file must contain columns: genome, completeness, contamination")
        return False

    bins_map = load_bins_map(output_dir)
    if not bins_map:
        print(f"{ERROR}ERROR:{RESET} bins_to_files.json not found; cannot validate quality file.")
        return False

    alias_to_base = {}
    for bin_id, path in bins_map.items():
        base = os.path.basename(path)
        base_without_gz = base[:-3] if base.endswith(".gz") else base
        canonical = base_without_gz
        aliases = {
            str(bin_id).strip(),
            normalize_genome_name(bin_id),
            base,
            base_without_gz,
            normalize_genome_name(base),
            normalize_genome_name(base_without_gz),
        }
        for alias in aliases:
            if alias:
                alias_to_base[str(alias)] = canonical

    mapped = []
    missing = []
    for name in df["genome"].astype(str):
        key = os.path.basename(name.strip())
        candidates = [
            key,
            key[:-3] if key.endswith(".gz") else None,
            normalize_genome_name(key),
        ]
        mapped_base = None
        for candidate in candidates:
            if candidate and candidate in alias_to_base:
                mapped_base = alias_to_base[candidate]
                break
        if mapped_base:
            mapped.append(mapped_base)
        else:
            mapped.append(key)
            missing.append(key)

    if missing:
        print(f"{ERROR}ERROR:{RESET} Quality file missing bins: {missing}")
        return False

    out_dir = Path(output_dir) / "cataloging" / "final"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "all_bin_metadata.csv"
    out_df = df.copy()
    out_df["genome"] = mapped
    out_df = out_df[["genome", "completeness", "contamination"]]
    out_df.to_csv(out_path, index=False)
    return True
