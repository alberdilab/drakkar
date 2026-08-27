import gzip
import json
import tempfile
import unittest
from pathlib import Path

from drakkar.amr_inputs import (
    assemblies_from_directory,
    assemblies_from_table,
    inspect_assembly,
    write_amr_assemblies_manifest,
)
from drakkar.input_errors import InputFileError


class AmrInputTests(unittest.TestCase):
    @staticmethod
    def write_fasta(path, contig="contig_1", sequence="ACGTACGT"):
        path.parent.mkdir(parents=True, exist_ok=True)
        text = f">{contig} description\n{sequence}\n"
        if path.suffix == ".gz":
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(text)
        else:
            path.write_text(text, encoding="utf-8")

    def test_flat_and_nested_directories_are_resolved_without_recursive_guessing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            flat = Path(tmpdir) / "flat"
            self.write_fasta(flat / "sample_a.fna")
            self.write_fasta(flat / "sample_b.fa.gz", contig="b")
            resolved = assemblies_from_directory(flat)
            self.assertEqual(set(resolved), {"sample_a", "sample_b"})
            self.assertEqual(resolved["sample_a"]["contig_count"], 1)
            self.assertEqual(resolved["sample_a"]["total_length"], 8)

            nested = Path(tmpdir) / "nested"
            self.write_fasta(nested / "assembly_1" / "assembly_1.fna")
            self.write_fasta(nested / "assembly_2" / "final.contigs.fa")
            resolved = assemblies_from_directory(nested, default_type="isolate")
            self.assertEqual(set(resolved), {"assembly_1", "assembly_2"})
            self.assertEqual(resolved["assembly_1"]["assembly_type"], "isolate")

    def test_drakkar_cataloging_megahit_layout_is_detected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self.write_fasta(
                root / "cataloging" / "megahit" / "sample" / "sample.fna"
            )
            resolved = assemblies_from_directory(root)
            self.assertEqual(list(resolved), ["sample"])
            self.assertEqual(
                resolved["sample"]["source_layout"], "drakkar_cataloging_megahit"
            )

    def test_manifest_paths_are_relative_to_manifest_and_nan_cells_use_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self.write_fasta(root / "assemblies" / "a.fna")
            table = root / "assemblies.tsv"
            table.write_text(
                "assembly_id\tassembly_path\tassembly_type\torganism\n"
                "001\tassemblies/a.fna\t\t\n",
                encoding="utf-8",
            )
            resolved = assemblies_from_table(table, default_type="metagenome")
            self.assertEqual(resolved["001"]["assembly_type"], "metagenome")
            self.assertEqual(resolved["001"]["organism"], "")
            self.assertEqual(resolved["001"]["path"], str((root / "assemblies" / "a.fna").resolve()))

            output = root / "output"
            write_amr_assemblies_manifest(output=output, table=table)
            written = json.loads(
                (output / "data" / "amr_assemblies.json").read_text(encoding="utf-8")
            )
            self.assertEqual(written["001"]["sha256"], resolved["001"]["sha256"])

    def test_duplicate_contig_ids_and_ambiguous_nested_fastas_fail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            duplicated = root / "duplicate.fna"
            duplicated.write_text(">same one\nAAAA\n>same two\nTTTT\n", encoding="utf-8")
            with self.assertRaisesRegex(InputFileError, "duplicate first-token"):
                inspect_assembly(duplicated)

            collection = root / "collection"
            self.write_fasta(collection / "ambiguous" / "one.fa")
            self.write_fasta(collection / "ambiguous" / "two.fasta", contig="two")
            with self.assertRaisesRegex(InputFileError, "multiple FASTA"):
                assemblies_from_directory(collection)

            invalid = root / "protein.fna"
            invalid.write_text(">protein\nMPEPTIDE\n", encoding="utf-8")
            with self.assertRaisesRegex(InputFileError, "non-IUPAC nucleotide"):
                inspect_assembly(invalid)


if __name__ == "__main__":
    unittest.main()
