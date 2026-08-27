#!/usr/bin/env python3
"""Read and optionally validate the release version embedded in card.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def card_version(path: str | Path) -> str:
    with Path(path).open(encoding="utf-8") as handle:
        payload = json.load(handle)
    version = str(payload.get("_version") or "").strip()
    if not version:
        raise ValueError(f"CARD database has no _version value: {path}")
    return version


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("card_json")
    parser.add_argument("--expect")
    args = parser.parse_args()

    version = card_version(args.card_json)
    if args.expect and version != args.expect:
        parser.error(
            f"CARD archive contains version {version}, but version {args.expect} was requested"
        )
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
