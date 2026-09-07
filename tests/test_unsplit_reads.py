from __future__ import annotations

import email.message
import gzip
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.request import Request

from drakkar.fastq_split import read_index, record_stem, split_unsplit_fastq
from drakkar.input_errors import FastqSplitError
from drakkar.utils import file_samples_to_json


class FakeResponse(io.BytesIO):
    """A urlopen response that can advertise headers, like a real HTTPResponse."""

    def __init__(self, payload: bytes = b"", headers: dict | None = None) -> None:
        super().__init__(payload)
        self.headers = email.message.Message()
        for name, value in (headers or {}).items():
            self.headers[name] = str(value)

    def getheader(self, name, default=None):
        return self.headers.get(name, default)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def _record(name: str, index: int, length: int) -> str:
    """One FASTQ record in the shape the archive publishes for an unsplit run."""
    return f"@{name} {name.rsplit('.', 1)[-1]}/{index}\n{'A' * length}\n+\n{'J' * length}\n"


def _unsplit_payload(pairs: int, first_index: int = 1, accession: str = "SRR9851002") -> bytes:
    """A flat FASTQ of two halves, the ``first_index`` half written first.

    Reads are variable length and spot numbers run straight through both
    halves, exactly as the SRA loader writes them when it stores every read as
    its own single-read spot.
    """
    second_index = 2 if first_index == 1 else 1
    records = []
    spot = 1
    for index in (first_index, second_index):
        for offset in range(pairs):
            records.append(_record(f"{accession}.{spot}", index, 50 + offset))
            spot += 1
    return gzip.compress("".join(records).encode("utf-8"))


def _headers(payload: bytes) -> dict:
    return {"Content-Length": len(payload)}


def _read_names(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.strip() for i, line in enumerate(handle) if i % 4 == 0]


def _sheet(tmpdir: Path, rows: str, header: str) -> Path:
    infofile = tmpdir / "info.tsv"
    infofile.write_text(header + rows, encoding="utf-8")
    return infofile


UNSPLIT_URL = (
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR985/002/SRR9851002/SRR9851002.fastq.gz"
)
UNSPLIT_HTTPS_URL = UNSPLIT_URL.replace("ftp://", "https://", 1)


class ReadNameParsingTests(unittest.TestCase):
    def test_read_index_is_taken_from_the_ena_mate_field(self) -> None:
        self.assertEqual(read_index(b"@SRR9851002.1 1/1"), 1)
        self.assertEqual(read_index(b"@SRR9851002.7 7/2"), 2)

    def test_read_index_falls_back_to_a_suffix_on_the_read_name(self) -> None:
        self.assertEqual(read_index(b"@V300026712L2C001R0010000372/1"), 1)
        self.assertEqual(read_index(b"@V300026712L2C001R0010000372/2"), 2)

    def test_read_index_is_none_without_a_mate_suffix(self) -> None:
        self.assertIsNone(read_index(b"@SRR9851002.1 1"))
        self.assertIsNone(read_index(b"@SRR9851002.1 length=53"))
        self.assertIsNone(read_index(b"@"))

    def test_record_stem_drops_the_spot_number_and_mate_suffix(self) -> None:
        self.assertEqual(record_stem(b"@SRR9851002.1 1/1"), b"SRR9851002")
        self.assertEqual(record_stem(b"@V300026712L2C001R0010000372/1"), b"V300026712L2C001R0010000372")


class SplitUnsplitFastqTests(unittest.TestCase):
    def _split(self, tmpdir: Path, payload: bytes, **kwargs):
        source = tmpdir / "SRR9851002.fastq.gz"
        source.write_bytes(payload)
        read1 = tmpdir / "sample_1.fq.gz"
        read2 = tmpdir / "sample_2.fq.gz"
        count = split_unsplit_fastq(str(source), str(read1), str(read2), **kwargs)
        return count, read1, read2

    def test_records_are_routed_on_their_read_index_when_read_one_comes_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            count, read1, read2 = self._split(tmpdir, _unsplit_payload(3, first_index=1))

            self.assertEqual(count, 3)
            self.assertEqual(
                _read_names(read1),
                ["@SRR9851002.1 1/1", "@SRR9851002.2 2/1", "@SRR9851002.3 3/1"],
            )
            self.assertEqual(
                _read_names(read2),
                ["@SRR9851002.1 1/2", "@SRR9851002.2 2/2", "@SRR9851002.3 3/2"],
            )

    def test_records_are_routed_on_their_read_index_when_read_two_comes_first(self) -> None:
        # Which half leads varies from run to run, so position must not decide.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            count, read1, read2 = self._split(tmpdir, _unsplit_payload(3, first_index=2))

            self.assertEqual(count, 3)
            with gzip.open(read1, "rt") as handle:
                forward = handle.read()
            with gzip.open(read2, "rt") as handle:
                reverse = handle.read()

            self.assertEqual(forward.count("/1"), 3)
            self.assertNotIn("/2", forward)
            self.assertEqual(reverse.count("/2"), 3)
            self.assertNotIn("/1", reverse)

    def test_mates_are_renumbered_so_both_halves_share_a_read_id(self) -> None:
        # The archive numbers every read with its own spot, which leaves the two
        # halves with no name in common; seqkit pair would match nothing.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _, read1, read2 = self._split(tmpdir, _unsplit_payload(4, first_index=2))

            ids1 = [name.split()[0] for name in _read_names(read1)]
            ids2 = [name.split()[0] for name in _read_names(read2)]
            self.assertEqual(ids1, ids2)
            self.assertEqual(ids1, [f"@SRR9851002.{n}" for n in (1, 2, 3, 4)])

    def test_sequences_and_qualities_travel_with_their_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _, read1, read2 = self._split(tmpdir, _unsplit_payload(2, first_index=1))

            with gzip.open(read1, "rt") as handle:
                lines = handle.read().splitlines()
            self.assertEqual(lines[1], "A" * 50)
            self.assertEqual(lines[3], "J" * 50)
            self.assertEqual(lines[5], "A" * 51)
            with gzip.open(read2, "rt") as handle:
                reverse = handle.read().splitlines()
            self.assertEqual(reverse[1], "A" * 50)
            self.assertEqual(reverse[3], "J" * 50)
            self.assertEqual(reverse[5], "A" * 51)

    def test_a_plain_uncompressed_source_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            payload = gzip.decompress(_unsplit_payload(2))
            count, read1, read2 = self._split(tmpdir, payload)

            self.assertEqual(count, 2)
            self.assertEqual(len(_read_names(read1)), 2)
            self.assertEqual(len(_read_names(read2)), 2)

    def test_a_crlf_source_does_not_leak_carriage_returns_into_the_reads(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            text = gzip.decompress(_unsplit_payload(2)).decode("utf-8")
            payload = gzip.compress(text.replace("\n", "\r\n").encode("utf-8"))
            count, read1, read2 = self._split(tmpdir, payload)

            self.assertEqual(count, 2)
            with gzip.open(read1, "rb") as handle:
                forward = handle.read()
            self.assertNotIn(b"\r", forward)
            self.assertEqual(forward.splitlines()[1], b"A" * 50)

    def test_unequal_halves_are_rejected_and_leave_no_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            payload = gzip.compress(
                (
                    _record("SRR9851002.1", 1, 50)
                    + _record("SRR9851002.2", 1, 50)
                    + _record("SRR9851002.3", 2, 50)
                ).encode("utf-8")
            )

            with self.assertRaises(FastqSplitError) as raised:
                self._split(tmpdir, payload)

            self.assertIn("2 forward and 1 reverse", str(raised.exception))
            self.assertFalse((tmpdir / "sample_1.fq.gz").exists())
            self.assertFalse((tmpdir / "sample_2.fq.gz").exists())

    def test_a_record_without_a_mate_suffix_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            payload = gzip.compress(
                (
                    _record("SRR9851002.1", 1, 50)
                    + "@SRR9851002.2 2\nAAAA\n+\nJJJJ\n"
                    + _record("SRR9851002.3", 2, 50)
                ).encode("utf-8")
            )

            with self.assertRaises(FastqSplitError) as raised:
                self._split(tmpdir, payload)

            self.assertIn("carries no /1 or /2 read index", str(raised.exception))
            self.assertFalse((tmpdir / "sample_1.fq.gz").exists())

    def test_a_file_holding_only_one_mate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            payload = gzip.compress(
                (_record("SRR9851002.1", 1, 50) + _record("SRR9851002.2", 1, 50)).encode("utf-8")
            )

            with self.assertRaises(FastqSplitError) as raised:
                self._split(tmpdir, payload)

            self.assertIn("holds only one mate", str(raised.exception))
            self.assertFalse((tmpdir / "sample_1.fq.gz").exists())
            self.assertFalse((tmpdir / "sample_2.fq.gz").exists())

    def test_a_truncated_final_record_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            payload = gzip.compress(
                (_record("SRR9851002.1", 1, 50) + "@SRR9851002.2 2/2\nAAAA\n").encode("utf-8")
            )

            with self.assertRaises(FastqSplitError) as raised:
                self._split(tmpdir, payload)

            self.assertIn("truncated", str(raised.exception))

    def test_an_empty_source_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with self.assertRaises(FastqSplitError) as raised:
                self._split(tmpdir, gzip.compress(b""))

            self.assertIn("no FASTQ records", str(raised.exception))

    def test_partial_output_is_written_to_a_configured_temporary_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            scratch = tmpdir / "scratch"

            count, read1, read2 = self._split(
                tmpdir, _unsplit_payload(2), tmp_dir=str(scratch)
            )

            self.assertEqual(count, 2)
            self.assertTrue(read1.exists())
            self.assertTrue(read2.exists())
            # Promoted, so nothing is left behind in the scratch directory.
            self.assertEqual(sorted(os.listdir(scratch)), [])
            self.assertFalse((tmpdir / "sample_1.fq.gz.tmp").exists())

    def test_drakkar_tmpdir_selects_the_temporary_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            scratch = tmpdir / "elsewhere"
            source = tmpdir / "SRR9851002.fastq.gz"
            source.write_bytes(_unsplit_payload(2))
            infofile = _sheet(
                tmpdir,
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{source}\n",
                UnsplitSampleSheetTests.HEADER,
            )

            with patch.dict(os.environ, {"DRAKKAR_TMPDIR": str(scratch)}):
                file_samples_to_json(str(infofile), str(tmpdir))

            self.assertTrue(scratch.is_dir(), "the configured directory must be used")
            self.assertEqual(sorted(os.listdir(scratch)), [])
            cache = tmpdir / "data" / "reads_cache"
            self.assertTrue((cache / "SRR9851002_SRR9851002_1.fq.gz").exists())
            self.assertEqual([p.name for p in cache.glob("*.tmp")], [])

    def test_a_failed_split_leaves_the_configured_temporary_directory_clean(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            scratch = tmpdir / "scratch"
            payload = gzip.compress(
                (_record("SRR9851002.1", 1, 50) + _record("SRR9851002.2", 1, 50)).encode("utf-8")
            )

            with self.assertRaises(FastqSplitError):
                self._split(tmpdir, payload, tmp_dir=str(scratch))

            self.assertEqual(sorted(os.listdir(scratch)), [])


class UnsplitSampleSheetTests(unittest.TestCase):
    HEADER = "sample\trawreads1\trawreads2\treference_name\treference_path\trawreads_unsplit\n"

    def test_an_unsplit_row_is_downloaded_split_and_fed_in_as_a_read_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = _unsplit_payload(3, first_index=2)
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(headers=_headers(payload)),
                    FakeResponse(payload, headers=_headers(payload)),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            cache = Path(tmpdir) / "data" / "reads_cache"
            expected1 = cache / "SRR9851002_SRR9851002_1.fq.gz"
            expected2 = cache / "SRR9851002_SRR9851002_2.fq.gz"

            reads1 = json.loads((Path(tmpdir) / "data" / "sample_to_reads1.json").read_text())
            reads2 = json.loads((Path(tmpdir) / "data" / "sample_to_reads2.json").read_text())
            self.assertEqual(reads1["SRR9851002"], [str(expected1)])
            self.assertEqual(reads2["SRR9851002"], [str(expected2)])
            self.assertEqual(len(_read_names(expected1)), 3)
            self.assertEqual(len(_read_names(expected2)), 3)

    def test_the_downloaded_flat_file_is_removed_once_the_halves_are_written(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = _unsplit_payload(2)
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(headers=_headers(payload)),
                    FakeResponse(payload, headers=_headers(payload)),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            cache = Path(tmpdir) / "data" / "reads_cache"
            self.assertFalse((cache / "SRR9851002_SRR9851002.fastq.gz").exists())
            self.assertEqual(
                sorted(p.name for p in cache.iterdir()),
                ["SRR9851002_SRR9851002_1.fq.gz", "SRR9851002_SRR9851002_2.fq.gz"],
            )

    def test_existing_halves_are_reused_without_downloading_or_splitting_again(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = _unsplit_payload(2)
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(headers=_headers(payload)),
                    FakeResponse(payload, headers=_headers(payload)),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            # A resumed run must not touch the network at all.
            with patch("drakkar.utils.urlopen", side_effect=AssertionError("downloaded again")):
                file_samples_to_json(str(infofile), tmpdir)

            reads1 = json.loads((Path(tmpdir) / "data" / "sample_to_reads1.json").read_text())
            expected1 = Path(tmpdir) / "data" / "reads_cache" / "SRR9851002_SRR9851002_1.fq.gz"
            self.assertEqual(reads1["SRR9851002"], [str(expected1)])

    def test_an_ena_ftp_link_is_fetched_over_https(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = _unsplit_payload(2)
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(headers=_headers(payload)),
                    FakeResponse(payload, headers=_headers(payload)),
                ],
            ) as fake_urlopen:
                file_samples_to_json(str(infofile), tmpdir)

            requested = [
                call.args[0].full_url
                if isinstance(call.args[0], Request)
                else call.args[0]
                for call in fake_urlopen.call_args_list
            ]
            self.assertEqual(requested, [UNSPLIT_HTTPS_URL, UNSPLIT_HTTPS_URL])

    def test_a_local_unsplit_file_is_split_and_kept(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "SRR9851002.fastq.gz"
            source.write_bytes(_unsplit_payload(2))
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{source}\n",
                self.HEADER,
            )

            file_samples_to_json(str(infofile), tmpdir)

            self.assertTrue(source.exists(), "a user's own file must never be removed")
            cache = Path(tmpdir) / "data" / "reads_cache"
            self.assertTrue((cache / "SRR9851002_SRR9851002_1.fq.gz").exists())

    def test_a_sheet_may_mix_unsplit_rows_with_ordinary_read_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            read1 = Path(tmpdir) / "S001_1.fq.gz"
            read2 = Path(tmpdir) / "S001_2.fq.gz"
            read1.write_bytes(gzip.compress(_record("S001.1", 1, 50).encode("utf-8")))
            read2.write_bytes(gzip.compress(_record("S001.1", 2, 50).encode("utf-8")))
            payload = _unsplit_payload(2)
            infofile = _sheet(
                Path(tmpdir),
                f"S001\t{read1}\t{read2}\tref\t/refs/ref.fa\t\n"
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(headers=_headers(payload)),
                    FakeResponse(payload, headers=_headers(payload)),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            reads1 = json.loads((Path(tmpdir) / "data" / "sample_to_reads1.json").read_text())
            # The ordinary row keeps resolving to the user's own file, untouched.
            self.assertEqual(reads1["S001"], [str(read1.resolve())])
            self.assertIn("SRR9851002", reads1)

    def test_a_sheet_without_the_column_is_unaffected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            read1 = Path(tmpdir) / "S001_1.fq.gz"
            read2 = Path(tmpdir) / "S001_2.fq.gz"
            read1.write_bytes(gzip.compress(b"@S001.1 1/1\nAAAA\n+\nJJJJ\n"))
            read2.write_bytes(gzip.compress(b"@S001.1 1/2\nAAAA\n+\nJJJJ\n"))
            infofile = _sheet(
                Path(tmpdir),
                f"S001\t{read1}\t{read2}\tref\t/refs/ref.fa\n",
                "sample\trawreads1\trawreads2\treference_name\treference_path\n",
            )

            file_samples_to_json(str(infofile), tmpdir)

            reads1 = json.loads((Path(tmpdir) / "data" / "sample_to_reads1.json").read_text())
            reads2 = json.loads((Path(tmpdir) / "data" / "sample_to_reads2.json").read_text())
            self.assertEqual(reads1["S001"], [str(read1.resolve())])
            self.assertEqual(reads2["S001"], [str(read2.resolve())])

    def test_an_unsplit_value_alongside_rawreads1_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            read1 = Path(tmpdir) / "S001_1.fq.gz"
            read2 = Path(tmpdir) / "S001_2.fq.gz"
            read1.write_bytes(gzip.compress(b"@S001.1 1/1\nAAAA\n+\nJJJJ\n"))
            read2.write_bytes(gzip.compress(b"@S001.1 1/2\nAAAA\n+\nJJJJ\n"))
            infofile = _sheet(
                Path(tmpdir),
                f"S001\t{read1}\t{read2}\tref\t/refs/ref.fa\t{UNSPLIT_URL}\n",
                self.HEADER,
            )

            with self.assertRaises(SystemExit):
                file_samples_to_json(str(infofile), tmpdir)

    def test_a_failed_split_stops_the_run_with_the_sample_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "SRR9851002.fastq.gz"
            source.write_bytes(
                gzip.compress(
                    (
                        _record("SRR9851002.1", 1, 50) + _record("SRR9851002.2", 1, 50)
                    ).encode("utf-8")
                )
            )
            infofile = _sheet(
                Path(tmpdir),
                f"SRR9851002\t\t\tref\t/refs/ref.fa\t{source}\n",
                self.HEADER,
            )

            with self.assertRaises(SystemExit):
                file_samples_to_json(str(infofile), tmpdir)


if __name__ == "__main__":
    unittest.main()
