import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts" / "prepare_amrfinder_gff.py"
SPEC = importlib.util.spec_from_file_location("prepare_amrfinder_gff", SCRIPT)
gff = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gff)


class PrepareAmrfinderGffTests(unittest.TestCase):
    def test_name_is_added_from_id_without_changing_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "prodigal.gff"
            output = Path(tmpdir) / "amrfinder.gff"
            source.write_text(
                "##gff-version  3\n"
                "contig_1\tProdigal_v2.6.3\tCDS\t4\t99\t.\t+\t0\tID=contig_1_1;partial=00\n",
                encoding="utf-8",
            )
            self.assertEqual(gff.convert(source, output), 1)
            converted = output.read_text(encoding="utf-8")
            self.assertIn("ID=contig_1_1;partial=00;Name=contig_1_1", converted)
            self.assertIn("\t4\t99\t", converted)


if __name__ == "__main__":
    unittest.main()
