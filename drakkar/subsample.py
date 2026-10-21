"""Parse read subsampling targets and calculate a non-upsampling fraction."""

import argparse
import re
from decimal import Decimal


_PAIRS = re.compile(r"([0-9]+)([KMG]?)")
_BASES = re.compile(r"([0-9]+(?:\.[0-9]+)?)([MG]b)")
# Snakemake --config turns 0.00001 into the float 1e-05, so exponents are accepted.
_FRACTION = re.compile(r"[0-9]*\.[0-9]+(?:e-[0-9]+)?|[0-9]+e-[0-9]+")


def parse_subsample_value(value):
    """Return (fraction|pairs|bases, numeric value) for a CLI target."""
    value = str(value)
    if _FRACTION.fullmatch(value):
        fraction = float(value)
        if 0 < fraction < 1:
            return "fraction", fraction
        raise argparse.ArgumentTypeError(
            "--subsample fraction must be greater than 0 and less than 1; "
            "omit --subsample to keep every read"
        )
    match = _BASES.fullmatch(value)
    if match:
        amount = Decimal(match.group(1)) * {"Mb": 10**6, "Gb": 10**9}[match.group(2)]
        if amount > 0 and amount == amount.to_integral_value():
            return "bases", int(amount)
    match = _PAIRS.fullmatch(value)
    if match:
        amount = int(match.group(1)) * {"": 1, "K": 10**3, "M": 10**6, "G": 10**9}[match.group(2)]
        if amount > 0:
            return "pairs", amount
    raise argparse.ArgumentTypeError(
        "--subsample must be a fraction between 0 and 1 (e.g. 0.25), "
        "a positive pair count (e.g. 50000000 or 50M), or a positive "
        "total-base target (e.g. 10Gb)"
    )


def subsample_fraction(value, pairs, bases):
    """Fraction of pairs to retain, capped at one for count/base targets."""
    kind, target = parse_subsample_value(value)
    if pairs < 0 or bases < 0:
        raise ValueError("read and base counts must be nonnegative")
    if kind == "fraction":
        return target
    total = pairs if kind == "pairs" else bases
    return min(1.0, target / total) if total else 1.0
