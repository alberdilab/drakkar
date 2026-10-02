from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "validate_tigrfam_database.py"


def load_module():
    spec = importlib.util.spec_from_file_location("validate_tigrfam_database", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class ValidateTIGRFAMDatabaseTests(unittest.TestCase):
    def test_exact_unversioned_accessions_and_complete_cutoffs_are_required(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            library = Path(tmpdir) / "tigrfams"
            library.write_text(
                "HMMER3/f\nNAME  merR\nACC   TIGR02053\nTC    100.0 80.0;\n//\n",
                encoding="utf-8",
            )
            self.assertEqual(module.read_library_accessions(library), {"TIGR02053"})

            library.write_text(
                "HMMER3/f\nNAME  versioned\nACC   TIGR02053.1\nTC    100.0 80.0;\n//\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "non-legacy accession"):
                module.read_library_accessions(library)

    def test_profile_without_cutoffs_is_rejected(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            library = Path(tmpdir) / "tigrfams"
            library.write_text(
                "HMMER3/f\nNAME  merR\nACC   TIGR02053\n//\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "has no complete trusted cutoff"):
                module.read_library_accessions(library)


if __name__ == "__main__":
    unittest.main()
