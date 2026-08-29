import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts" / "prepare_amrfinder_gff.py"
SPEC = importlib.util.spec_from_file_location("prepare_amrfinder_gff", SCRIPT)
gff = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gff)


PRODIGAL_GFF = (
    "##gff-version  3\n"
    "# Sequence Data: seqnum=1;seqlen=4993\n"
    "k141_100002-flag=0-multi=8.9637-len=4993\tProdigal_v2.6.3\tCDS\t4\t99\t.\t+\t0\t"
    "ID=1_1;partial=00;conf=99.99\n"
    "k141_100002-flag=0-multi=8.9637-len=4993\tProdigal_v2.6.3\tCDS\t120\t400\t.\t-\t0\t"
    "ID=1_2;partial=00;conf=99.99\n"
    "k141_2\tProdigal_v2.6.3\tCDS\t1\t80\t.\t+\t0\tID=2_1;partial=10;conf=99.99\n"
)

PRODIGAL_FAA = (
    ">k141_100002-flag=0-multi=8.9637-len=4993_1 # 4 # 99 # 1 # ID=1_1;partial=00\nMK*\n"
    ">k141_100002-flag=0-multi=8.9637-len=4993_2 # 120 # 400 # -1 # ID=1_2;partial=00\nMR*\n"
    ">k141_2_1 # 1 # 80 # 1 # ID=2_1;partial=10\nMA*\n"
)


class PrepareAmrfinderGffTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmpdir = Path(self._tmpdir.name)
        self.source = self.tmpdir / "prodigal.gff"
        self.output = self.tmpdir / "amrfinder.gff"
        self.proteins = self.tmpdir / "prodigal.faa"
        self.source.write_text(PRODIGAL_GFF, encoding="utf-8")
        self.proteins.write_text(PRODIGAL_FAA, encoding="utf-8")

    def test_identifiers_match_the_protein_fasta(self):
        self.assertEqual(gff.convert(self.source, self.output, self.proteins), 3)
        converted = self.output.read_text(encoding="utf-8")
        self.assertIn(
            "ID=k141_100002-flag=0-multi=8.9637-len=4993_1;partial=00;conf=99.99;"
            "Name=k141_100002-flag=0-multi=8.9637-len=4993_1",
            converted,
        )
        self.assertIn("ID=k141_2_1;partial=10;conf=99.99;Name=k141_2_1", converted)
        self.assertNotIn("ID=1_1;", converted)

    def test_coordinates_and_comments_are_preserved(self):
        gff.convert(self.source, self.output, self.proteins)
        converted = self.output.read_text(encoding="utf-8")
        self.assertIn("##gff-version  3\n", converted)
        self.assertIn("# Sequence Data: seqnum=1;seqlen=4993\n", converted)
        self.assertIn("\t120\t400\t.\t-\t", converted)

    def test_protein_absent_from_the_gff_is_reported(self):
        self.proteins.write_text(
            PRODIGAL_FAA + ">k141_2_2 # 200 # 280 # 1 # ID=2_2;partial=00\nMG*\n",
            encoding="utf-8",
        )
        with self.assertRaises(ValueError) as error:
            gff.convert(self.source, self.output, self.proteins)
        self.assertIn("k141_2_2", str(error.exception))


if __name__ == "__main__":
    unittest.main()
