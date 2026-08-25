from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "get_fvid.py"


def load_vfdb_mapping_module():
    spec = importlib.util.spec_from_file_location("get_fvid", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class VfdbMappingTests(unittest.TestCase):
    def test_current_set_b_header_uses_vfc_category_not_organism(self) -> None:
        module = load_vfdb_mapping_module()
        description = (
            "(plc1) phospholipase C "
            "[Phospholipase C (VF0470) - Exotoxin (VFC0235)] "
            "[Acinetobacter baumannii 1656-2]"
        )

        parsed = module.parse_header("VFG037170(gb|WP_001081754) " + description)

        self.assertEqual(parsed["entry"], "VFG037170(gb|WP_001081754)")
        self.assertEqual(parsed["vf"], description)
        self.assertEqual(parsed["vfc"], "VFC0235")
        self.assertEqual(parsed["vf_type"], "Exotoxin")

    def test_simple_vfc_section_is_supported(self) -> None:
        module = load_vfdb_mapping_module()

        parsed = module.parse_header(
            "VFG000001(gb|ABC123) adhesin [Adherence (VFC0001)] [Escherichia coli]"
        )

        self.assertEqual(parsed["vfc"], "VFC0001")
        self.assertEqual(parsed["vf_type"], "Adherence")

    def test_organism_only_header_does_not_create_a_vf_type(self) -> None:
        module = load_vfdb_mapping_module()

        parsed = module.parse_header(
            "VFG000002(gb|ABC124) hypothetical protein [Escherichia coli]"
        )

        self.assertEqual(parsed["vfc"], "")
        self.assertEqual(parsed["vf_type"], "")

    def test_nonstandard_vfc_text_preserves_id_without_guessing_category(self) -> None:
        module = load_vfdb_mapping_module()

        parsed = module.parse_header(
            "VFG000003(gb|ABC125) unclassified virulence factor VFC0999 [Organism strain]"
        )

        self.assertEqual(parsed["vfc"], "VFC0999")
        self.assertEqual(parsed["vf_type"], "")

    def test_mapping_output_declares_v2_schema(self) -> None:
        module = load_vfdb_mapping_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "vfdb.faa"
            output = Path(tmpdir) / "vfdb.tsv"
            fasta.write_text(
                ">VFG000001 adhesin [Adherence (VFC0001)] [Escherichia coli]\nMAAA\n",
                encoding="utf-8",
            )

            module.write_mapping(fasta, output)

            fields = output.read_text(encoding="utf-8").splitlines()[1].split("\t")

        self.assertEqual(fields[-1], "drakkar-vfdb-v2")


if __name__ == "__main__":
    unittest.main()
