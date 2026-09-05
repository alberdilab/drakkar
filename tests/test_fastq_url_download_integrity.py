from __future__ import annotations

import email.message
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.request import Request

from drakkar import downloads as downloads_module
from drakkar.downloads import (
    _download_urlopen,
    _normalize_ena_fastq_url,
    _upgrade_ena_ftp_url,
    download_to_cache,
)
from drakkar.input_errors import DownloadError
from drakkar.utils import DEFAULT_DOWNLOAD_RETRIES, file_samples_to_json


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


class DisconnectingResponse(FakeResponse):
    """Delivers one chunk, then fails the way a dropped connection does."""

    def __init__(self, payload: bytes, headers: dict | None = None) -> None:
        super().__init__(payload, headers)
        self._delivered = False

    def read(self, size=-1):
        if self._delivered:
            raise ConnectionResetError("connection reset by peer")
        self._delivered = True
        return super().read(size)


def _called_url(call) -> str:
    target = call.args[0]
    return target.full_url if isinstance(target, Request) else target


def _called_method(call) -> str:
    target = call.args[0]
    return target.get_method() if isinstance(target, Request) else "GET"


class ContentLengthVerificationTests(unittest.TestCase):
    """A body shorter than Content-Length must never be accepted."""

    def test_download_urlopen_rejects_body_shorter_than_content_length(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = str(Path(tmpdir) / "reads.fq.gz.tmp")

            with patch(
                "drakkar.downloads.urlopen",
                return_value=FakeResponse(b"1" * 40, {"Content-Length": 100}),
            ):
                with self.assertRaises(DownloadError) as ctx:
                    _download_urlopen("https://example.org/reads.fq.gz", tmp_path)

            self.assertIn("truncated", str(ctx.exception))
            self.assertIn("100", str(ctx.exception))
            self.assertIn("40", str(ctx.exception))

    def test_download_urlopen_accepts_body_matching_content_length(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = str(Path(tmpdir) / "reads.fq.gz.tmp")

            with patch(
                "drakkar.downloads.urlopen",
                return_value=FakeResponse(b"1" * 100, {"Content-Length": 100}),
            ):
                _download_urlopen("https://example.org/reads.fq.gz", tmp_path)

            self.assertEqual(Path(tmp_path).stat().st_size, 100)

    def test_download_urlopen_ignores_content_length_of_encoded_body(self) -> None:
        # Content-Length then describes the encoded stream, not the bytes on disk.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = str(Path(tmpdir) / "reads.fq.gz.tmp")

            with patch(
                "drakkar.downloads.urlopen",
                return_value=FakeResponse(
                    b"1" * 40, {"Content-Length": 100, "Content-Encoding": "gzip"}
                ),
            ):
                _download_urlopen("https://example.org/reads.fq.gz", tmp_path)

            self.assertEqual(Path(tmp_path).stat().st_size, 40)

    def test_download_urlopen_verifies_without_any_caller_supplied_size(self) -> None:
        # The backstop is universal: no expected_size is passed anywhere here.
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "drakkar.downloads.urlopen",
                side_effect=[
                    FakeResponse(b"1" * 40, {"Content-Length": 100}),
                    FakeResponse(b"1" * 100, {"Content-Length": 100}),
                ],
            ) as urlopen_mock, patch("drakkar.downloads.time.sleep"):
                dest = download_to_cache(
                    "https://example.org/reads.fq.gz", "sample1", "rawreads1", tmpdir
                )

            self.assertEqual(urlopen_mock.call_count, 2)
            self.assertEqual(Path(dest).stat().st_size, 100)


class TemporaryFilePromotionTests(unittest.TestCase):
    """A .tmp file from a failed attempt must never reach the destination."""

    def test_truncated_download_never_promotes_tmp_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "drakkar.downloads.urlopen",
                side_effect=[
                    FakeResponse(b"1" * 40, {"Content-Length": 100})
                    for _ in range(DEFAULT_DOWNLOAD_RETRIES)
                ],
            ), patch("drakkar.downloads.time.sleep"):
                with self.assertRaises(DownloadError):
                    download_to_cache(
                        "https://example.org/reads.fq.gz", "sample1", "rawreads1", tmpdir
                    )

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertFalse((cache_dir / "sample1_reads.fq.gz").exists())
            self.assertEqual(list(cache_dir.glob("*.tmp")), [])

    def test_midstream_disconnect_never_promotes_tmp_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "drakkar.downloads.urlopen",
                side_effect=[
                    DisconnectingResponse(b"1" * 100, {"Content-Length": 100})
                    for _ in range(DEFAULT_DOWNLOAD_RETRIES)
                ],
            ), patch("drakkar.downloads.time.sleep"):
                with self.assertRaises(DownloadError):
                    download_to_cache(
                        "https://example.org/reads.fq.gz", "sample1", "rawreads1", tmpdir
                    )

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertFalse((cache_dir / "sample1_reads.fq.gz").exists())
            self.assertEqual(list(cache_dir.glob("*.tmp")), [])

    def test_stale_tmp_file_is_discarded_rather_than_resumed(self) -> None:
        real_download_once = downloads_module._download_once
        tmp_existed_at_download = []

        def recording_download_once(url, tmp_path):
            tmp_existed_at_download.append(Path(tmp_path).exists())
            return real_download_once(url, tmp_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            cache_dir.mkdir(parents=True)
            # A partial file left behind by a process that was killed mid-transfer.
            (cache_dir / "sample1_reads.fq.gz.tmp").write_bytes(b"9" * 60)

            with patch(
                "drakkar.downloads.urlopen",
                return_value=FakeResponse(b"1" * 100, {"Content-Length": 100}),
            ), patch("drakkar.downloads._download_once", recording_download_once):
                dest = download_to_cache(
                    "https://example.org/reads.fq.gz", "sample1", "rawreads1", tmpdir
                )

            # The transfer starts from nothing, so no stale byte can be promoted.
            self.assertEqual(tmp_existed_at_download, [False])
            self.assertEqual(Path(dest).read_bytes(), b"1" * 100)
            self.assertEqual(list(cache_dir.glob("*.tmp")), [])


class EnaFtpUrlUpgradeTests(unittest.TestCase):
    def test_ena_ftp_urls_are_upgraded_to_https(self) -> None:
        self.assertEqual(
            _upgrade_ena_ftp_url(
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR818/DRR818157/DRR818157_1.fastq.gz"
            ),
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR818/DRR818157/DRR818157_1.fastq.gz",
        )
        self.assertEqual(
            _upgrade_ena_ftp_url("ftp://ftp.ebi.ac.uk/pub/databases/ena/reads.fastq.gz"),
            "https://ftp.ebi.ac.uk/pub/databases/ena/reads.fastq.gz",
        )

    def test_non_ena_and_already_https_urls_are_left_alone(self) -> None:
        self.assertEqual(
            _upgrade_ena_ftp_url("ftp://ftp.example.org/reads.fastq.gz"),
            "ftp://ftp.example.org/reads.fastq.gz",
        )
        self.assertEqual(
            _upgrade_ena_ftp_url("https://ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz"),
            "https://ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz",
        )

    def test_normalize_ena_fastq_url_upgrades_schemed_and_schemeless_paths(self) -> None:
        self.assertEqual(
            _normalize_ena_fastq_url("ftp://ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz"),
            "https://ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz",
        )
        self.assertEqual(
            _normalize_ena_fastq_url("ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz"),
            "https://ftp.sra.ebi.ac.uk/vol1/reads.fastq.gz",
        )

    def test_url_read_column_is_downloaded_over_https(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\t"
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR818/DRR818157/DRR818157_1.fastq.gz\t"
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR818/DRR818157/DRR818157_2.fastq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b"", {"Content-Length": 100}),
                    FakeResponse(b"1" * 100, {"Content-Length": 100}),
                    FakeResponse(b"", {"Content-Length": 100}),
                    FakeResponse(b"2" * 100, {"Content-Length": 100}),
                ],
            ) as urlopen_mock:
                file_samples_to_json(str(infofile), tmpdir)

            requested = [_called_url(call) for call in urlopen_mock.call_args_list]
            self.assertTrue(
                all(url.startswith("https://ftp.sra.ebi.ac.uk/") for url in requested),
                requested,
            )
            self.assertEqual(
                [_called_method(call) for call in urlopen_mock.call_args_list],
                ["HEAD", "GET", "HEAD", "GET"],
            )


class UrlReadExpectedSizeTests(unittest.TestCase):
    """URL columns must get an expected size, exactly as the accession path does."""

    def test_head_probe_supplies_expected_size_and_rejects_short_download(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b"", {"Content-Length": 8}),
                    # A truncated body with no Content-Length of its own: only the
                    # HEAD-derived expected size can catch this one.
                    FakeResponse(b"bad"),
                    FakeResponse(b"12345678"),
                    FakeResponse(b"", {"Content-Length": 8}),
                    FakeResponse(b"87654321"),
                ],
            ) as urlopen_mock, patch("drakkar.utils.time.sleep"):
                file_samples_to_json(str(infofile), tmpdir)

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertEqual(urlopen_mock.call_count, 5)
            self.assertEqual(
                (cache_dir / "sample1_sample1_1.fq.gz").read_bytes(), b"12345678"
            )
            self.assertEqual(
                (cache_dir / "sample1_sample1_2.fq.gz").read_bytes(), b"87654321"
            )

    def test_head_probe_degrades_gracefully_when_server_reports_no_length(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    OSError("HEAD not allowed"),
                    FakeResponse(b"1" * 100),
                    FakeResponse(b"", {}),
                    FakeResponse(b"2" * 100),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            reads1 = json.loads(
                (Path(tmpdir) / "data" / "sample_to_reads1.json").read_text(encoding="utf-8")
            )
            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertEqual(reads1["sample1"], [str(cache_dir / "sample1_sample1_1.fq.gz")])
            self.assertEqual((cache_dir / "sample1_sample1_1.fq.gz").stat().st_size, 100)

    def test_cached_url_read_is_redownloaded_when_size_disagrees_with_server(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )
            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            cache_dir.mkdir(parents=True)
            # A truncated file left behind by the buggy download path.
            (cache_dir / "sample1_sample1_1.fq.gz").write_bytes(b"trunc")
            (cache_dir / "sample1_sample1_2.fq.gz").write_bytes(b"87654321")

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b"", {"Content-Length": 8}),
                    FakeResponse(b"12345678", {"Content-Length": 8}),
                    FakeResponse(b"", {"Content-Length": 8}),
                ],
            ) as urlopen_mock:
                file_samples_to_json(str(infofile), tmpdir)

            self.assertEqual(urlopen_mock.call_count, 3)
            self.assertEqual(
                (cache_dir / "sample1_sample1_1.fq.gz").read_bytes(), b"12345678"
            )
            self.assertEqual(
                (cache_dir / "sample1_sample1_2.fq.gz").read_bytes(), b"87654321"
            )


class UrlPairedSizeBalanceTests(unittest.TestCase):
    """The accession path's fallback guardrail must also cover URL columns."""

    def test_imbalanced_url_pair_is_rejected_and_partial_files_removed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b""),
                    FakeResponse(b"1" * 100),
                    FakeResponse(b""),
                    FakeResponse(b"2" * 80),
                ],
            ):
                with self.assertRaises(SystemExit):
                    file_samples_to_json(str(infofile), tmpdir)

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertFalse((cache_dir / "sample1_sample1_1.fq.gz").exists())
            self.assertFalse((cache_dir / "sample1_sample1_2.fq.gz").exists())
            self.assertFalse((Path(tmpdir) / "data" / "sample_to_reads1.json").exists())

    def test_balanced_url_pair_without_expected_sizes_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b""),
                    FakeResponse(b"1" * 100),
                    FakeResponse(b""),
                    FakeResponse(b"2" * 96),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertTrue((cache_dir / "sample1_sample1_1.fq.gz").exists())
            self.assertTrue((cache_dir / "sample1_sample1_2.fq.gz").exists())

    def test_exact_expected_sizes_skip_the_balance_guardrail(self) -> None:
        # Legitimately uneven mates must not be rejected once both sizes are known.
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                "sample1\thttps://example.org/sample1_1.fq.gz\t"
                "https://example.org/sample1_2.fq.gz\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b"", {"Content-Length": 100}),
                    FakeResponse(b"1" * 100, {"Content-Length": 100}),
                    FakeResponse(b"", {"Content-Length": 50}),
                    FakeResponse(b"2" * 50, {"Content-Length": 50}),
                ],
            ):
                file_samples_to_json(str(infofile), tmpdir)

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertEqual((cache_dir / "sample1_sample1_1.fq.gz").stat().st_size, 100)
            self.assertEqual((cache_dir / "sample1_sample1_2.fq.gz").stat().st_size, 50)

    def test_local_mate_is_never_deleted_when_the_guardrail_fires(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            local_read2 = Path(tmpdir) / "sample1_2.fq.gz"
            local_read2.write_bytes(b"2" * 80)
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\trawreads1\trawreads2\n"
                f"sample1\thttps://example.org/sample1_1.fq.gz\t{local_read2}\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[FakeResponse(b""), FakeResponse(b"1" * 100)],
            ):
                with self.assertRaises(SystemExit):
                    file_samples_to_json(str(infofile), tmpdir)

            cache_dir = Path(tmpdir) / "data" / "reads_cache"
            self.assertFalse((cache_dir / "sample1_sample1_1.fq.gz").exists())
            self.assertTrue(local_read2.exists())
            self.assertEqual(local_read2.stat().st_size, 80)


if __name__ == "__main__":
    unittest.main()
