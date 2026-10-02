from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "validate_ncbifam_database.py"


def load_module():
    spec = importlib.util.spec_from_file_location("validate_ncbifam_database", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class ValidateNCBIfamDatabaseTests(unittest.TestCase):
    def test_profile_without_tc_is_rejected_during_install(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            library = Path(tmpdir) / "hmm_PGAP.LIB"
            library.write_text(
                "HMMER3/f\nNAME  missing_tc\nACC   NF040708.3\nLENG  144\n//\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "has no trusted cutoff.*does not use an E-value fallback"
            ):
                module.read_library_accessions(library)

    def test_exact_versioned_tigr_accession_is_valid_metadata(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata = Path(tmpdir) / "hmm_PGAP.tsv"
            metadata.write_text(
                "#ncbi_accession\tlabel\tsequence_cutoff\tdomain_cutoff\t"
                "hmm_length\tfamily_type\tsource\n"
                "TIGR04545.1\trSAM_ahbD_hemeb\t500\t500\t339\tequivalog\tJCVI\n",
                encoding="utf-8",
            )
            records = module.read_profile_cutoffs(metadata)

        self.assertEqual(list(records), ["TIGR04545.1"])
        self.assertNotIn("TIGR04545", records)


if __name__ == "__main__":
    unittest.main()
