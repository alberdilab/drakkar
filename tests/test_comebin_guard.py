"""Guards that keep COMEBin from killing a run on assemblies it cannot bin.

COMEBin seeds its clustering with single-copy marker genes. Assemblies too
small or too fragmented to carry them make that stage fail, which used to abort
the whole workflow. The rule now skips such assemblies up front and, when the
marker stage still fails, exports an empty contig-to-bin table instead of dying.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULES = ROOT / "drakkar" / "workflow" / "rules" / "cataloging.smk"
CONFIG = ROOT / "drakkar" / "workflow" / "config.yaml"

# Stand-ins for the run_comebin.sh outcomes the rule has to tell apart.
STUB_SUCCESS = """#!/usr/bin/env bash
while getopts a:p:t:o: OPT; do case $OPT in o) OUTDIR=$OPTARG;; a) ASM=$OPTARG;; esac; done
printf "contig\\n" > "$ASM.bacar_marker.2quarter_lencutoff_1001.seed"
mkdir -p "$OUTDIR/comebin_res"
printf "c1\\tbin1\\n" > "$OUTDIR/comebin_res/comebin_res.tsv"
"""

# COMEBin 1.1.0 aborts the stage and propagates a non-zero status.
STUB_MARKER_ABORT = """#!/usr/bin/env bash
echo "RuntimeError: HMMsearch completed but found zero marker hits: markers.hmmout" >&2
exit 1
"""

# COMEBin 1.0.4 logs the failed marker step and exits 0 without a result file.
STUB_MARKER_SILENT = """#!/usr/bin/env bash
while getopts a:p:t:o: OPT; do case $OPT in a) ASM=$OPTARG;; esac; done
echo "Hmmsearch failed! Not exist: $ASM.bacar_marker.hmmout"
exit 0
"""

# Markers ran but produced no seeds, so clustering crashes on something else.
STUB_EMPTY_SEED = """#!/usr/bin/env bash
while getopts a:p:t:o: OPT; do case $OPT in a) ASM=$OPTARG;; esac; done
: > "$ASM.bacar_marker.2quarter_lencutoff_1001.seed"
echo "IndexError in leiden clustering" >&2
exit 1
"""

STUB_OOM = """#!/usr/bin/env bash
echo "RuntimeError: CUDA out of memory" >&2
exit 137
"""

STUB_NO_OUTPUT = """#!/usr/bin/env bash
exit 0
"""


class _Namespace:
    """Attribute holder standing in for snakemake's input/output/params."""

    def __init__(self, **fields):
        self.__dict__.update(fields)


def comebin_shell() -> str:
    """Return the shell body of the `comebin` rule."""
    text = RULES.read_text(encoding="utf-8")
    rule = text.split("rule comebin:", 1)[1].split("rule comebin_table:", 1)[0]
    return textwrap.dedent(rule.split("shell:", 1)[1].split('"""', 2)[1])


class ComebinGuardTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.assembly = self.dir / "megahit" / "asm.fna"
        self.assembly.parent.mkdir(parents=True)
        self.outdir = self.dir / "comebin" / "asm"
        self.output = self.outdir / "comebin_res" / "comebin_res.tsv"

    def run_rule(self, stub: str, assembly_mb: int, skip_reason: str = ""):
        """Run the rule's shell body with a stubbed run_comebin.sh on PATH."""
        self.assembly.write_bytes(b">c1\n" + b"A" * (assembly_mb * 1024 * 1024))
        bam = self.dir / "sample.bam"
        bam.write_bytes(b"bam")

        bindir = self.dir / "bin"
        bindir.mkdir(exist_ok=True)
        launcher = bindir / "run_comebin.sh"
        launcher.write_text(stub)
        launcher.chmod(0o755)

        body = comebin_shell().format(
            input=_Namespace(assembly=str(self.assembly), bam=str(bam)),
            output=str(self.output),
            threads=8,
            wildcards=_Namespace(assembly="asm"),
            params=_Namespace(
                bamdir=str(self.dir / "comebin" / "asm_bams"),
                outdir=str(self.outdir),
                skip_reason=skip_reason,
            ),
        )
        # Snakemake runs shell bodies under `set -euo pipefail`.
        return subprocess.run(
            ["bash", "-c", "set -euo pipefail; " + body],
            capture_output=True,
            text=True,
            env=dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ['PATH']}"),
        )

    def test_a_skip_reason_shortcuts_the_run(self) -> None:
        result = self.run_rule(
            STUB_SUCCESS, assembly_mb=3, skip_reason="Assembly is smaller than 10 MB"
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.output.exists())
        self.assertEqual(self.output.stat().st_size, 0)
        self.assertIn("Assembly is smaller than 10 MB, skipping comebin", result.stdout)

    def test_large_enough_assembly_is_binned(self) -> None:
        result = self.run_rule(STUB_SUCCESS, assembly_mb=12)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), "c1\tbin1\n")

    def test_marker_abort_exports_an_empty_table(self) -> None:
        result = self.run_rule(STUB_MARKER_ABORT, assembly_mb=12)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.output.exists())
        self.assertEqual(self.output.stat().st_size, 0)

    def test_silent_marker_failure_exports_an_empty_table(self) -> None:
        result = self.run_rule(STUB_MARKER_SILENT, assembly_mb=12)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.output.exists())

    def test_empty_seed_file_exports_an_empty_table(self) -> None:
        result = self.run_rule(STUB_EMPTY_SEED, assembly_mb=12)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.output.exists())

    def test_resource_failures_still_propagate(self) -> None:
        """An OOM must keep failing so snakemake retries it with more memory."""
        result = self.run_rule(STUB_OOM, assembly_mb=12)

        self.assertEqual(result.returncode, 137)
        self.assertFalse(self.output.exists())

    def test_unexplained_missing_output_fails_loudly(self) -> None:
        result = self.run_rule(STUB_NO_OUTPUT, assembly_mb=12)

        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())
        self.assertIn("wrote no", result.stderr)

    def test_empty_marker_byproducts_are_cleared_before_running(self) -> None:
        """A stale empty byproduct would otherwise be reused on every retry."""
        self.assembly.write_bytes(b">c1\n" + b"A" * (12 * 1024 * 1024))
        stale = Path(str(self.assembly) + ".frag.faa")
        stale.touch()

        self.run_rule(STUB_SUCCESS, assembly_mb=12)

        self.assertFalse(stale.exists())


def load_skip_policy(min_binning_assembly_mb: int = 10) -> callable:
    """Exec `binning_skip_reason` out of the Snakemake rule file."""
    text = RULES.read_text(encoding="utf-8")
    block = text.split("def _file_size_mb(path):", 1)[1].split("\nrule ", 1)[0]
    namespace: dict = {"Path": Path, "MIN_BINNING_ASSEMBLY_MB": min_binning_assembly_mb}
    exec("def _file_size_mb(path):" + block, namespace)
    return namespace["binning_skip_reason"]


class BinningSkipPolicyTests(unittest.TestCase):
    """One policy decides which assemblies the binners skip."""

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.skip_reason = load_skip_policy()

    def write(self, size_mb: float) -> Path:
        path = self.dir / "assembly.fna"
        path.write_bytes(b"A" * int(size_mb * 1024 * 1024))
        return path

    def test_large_assembly_is_binned(self) -> None:
        self.assertEqual(self.skip_reason(self.write(25)), "")

    def test_small_assembly_names_the_threshold(self) -> None:
        self.assertEqual(
            self.skip_reason(self.write(4)), "Assembly is smaller than 10 MB"
        )

    def test_empty_assembly_is_reported_as_empty(self) -> None:
        self.assertEqual(self.skip_reason(self.write(0)), "Assembly is empty")

    def test_missing_assembly_does_not_raise(self) -> None:
        self.assertEqual(self.skip_reason(self.dir / "absent.fna"), "Assembly is empty")

    def test_threshold_is_configurable(self) -> None:
        """A 20 MB assembly is skipped once the threshold is raised above it."""
        skip_reason = load_skip_policy(min_binning_assembly_mb=50)

        self.assertEqual(skip_reason(self.write(20)), "Assembly is smaller than 50 MB")


class BinningThresholdConfigTests(unittest.TestCase):
    def test_every_binner_shares_one_policy(self) -> None:
        """All four binners gate on the same decision, none of them on its own copy."""
        text = RULES.read_text(encoding="utf-8")

        self.assertIn('MIN_BINNING_ASSEMBLY_MB = int(config.get("MIN_BINNING_ASSEMBLY_MB", 10))', text)
        self.assertEqual(
            text.count("skip_reason=lambda wildcards, input: binning_skip_reason(input.assembly)"), 4
        )
        self.assertEqual(text.count('if [ -n "{params.skip_reason}" ]; then'), 4)
        for binner in ("metabat2", "maxbin2", "semibin2", "comebin"):
            self.assertIn(f'echo "{{params.skip_reason}}, skipping {binner}..."', text)
        # No rule may hardcode the threshold it is meant to read from the policy.
        self.assertNotIn("< 10 )); then", text)

    def test_threshold_is_exposed_in_the_shipped_config(self) -> None:
        self.assertIn("MIN_BINNING_ASSEMBLY_MB: 10", CONFIG.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
