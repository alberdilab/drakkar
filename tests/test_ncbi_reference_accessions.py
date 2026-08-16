from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.input_errors import DownloadError
from drakkar.utils import (
    argument_references_to_json,
    file_references_to_json,
    is_ncbi_assembly_accession,
    resolve_ncbi_assembly_accession,
)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def listing_response(*assembly_dirs: str) -> FakeResponse:
    rows = "\n".join(
        f'<a href="{assembly_dir}/">{assembly_dir}/</a>' for assembly_dir in assembly_dirs
    )
    return FakeResponse(f"<html><body>{rows}</body></html>".encode("utf-8"))


class NcbiAssemblyAccessionTests(unittest.TestCase):
    def test_recognises_assembly_accessions(self) -> None:
        for accession in ("GCF_000001405.40", "GCA_000001405", "gcf_000001405.4"):
            self.assertTrue(is_ncbi_assembly_accession(accession), accession)

    def test_rejects_non_assembly_accessions(self) -> None:
        for value in (
            "",
            "host.fna",
            "/data/host.fna.gz",
            "https://example.org/ref.fna.gz",
            "ERR4303216",
            "GCF_00000140.40",
            "GCF_0000014050",
        ):
            self.assertFalse(is_ncbi_assembly_accession(value), value)

    def test_resolves_versioned_accession_to_genomic_fasta_url(self) -> None:
        with patch(
            "drakkar.utils.urlopen",
            return_value=listing_response(
                "GCF_000001405.39_GRCh38.p13", "GCF_000001405.40_GRCh38.p14"
            ),
        ):
            url = resolve_ncbi_assembly_accession("GCF_000001405.40")

        self.assertEqual(
            url,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        )

    def test_resolves_unversioned_accession_to_latest_version(self) -> None:
        with patch(
            "drakkar.utils.urlopen",
            return_value=listing_response(
                "GCA_000001405.9_GRCh38.p3", "GCA_000001405.29_GRCh38.p14"
            ),
        ):
            url = resolve_ncbi_assembly_accession("gca_000001405")

        self.assertEqual(
            url,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
            "GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz",
        )

    def test_reports_available_versions_when_requested_version_is_missing(self) -> None:
        with patch(
            "drakkar.utils.urlopen",
            return_value=listing_response("GCF_000001405.40_GRCh38.p14"),
        ):
            with self.assertRaises(DownloadError) as context:
                resolve_ncbi_assembly_accession("GCF_000001405.99")

        self.assertIn("GCF_000001405.40", str(context.exception))

    def test_reports_missing_assembly(self) -> None:
        with patch("drakkar.utils.urlopen", return_value=FakeResponse(b"<html></html>")):
            with self.assertRaises(DownloadError):
                resolve_ncbi_assembly_accession("GCF_123456789.1")

    def test_argument_references_to_json_downloads_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_reads = Path(tmpdir) / "sample_to_reads1.json"
            sample_reads.write_text('{"sample1": ["reads_1.fq.gz"]}', encoding="utf-8")

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    listing_response("GCF_000001405.40_GRCh38.p14"),
                    FakeResponse(b">ref\nACGT\n"),
                ],
            ):
                argument_references_to_json("GCF_000001405.40", str(sample_reads), tmpdir)

            output_json = Path(tmpdir) / "data" / "reference_to_file.json"
            reference_to_file = json.loads(output_json.read_text(encoding="utf-8"))
            expected_path = (
                Path(tmpdir)
                / "data"
                / "references_cache"
                / "reference_GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
            )

            self.assertEqual(reference_to_file["reference"], [str(expected_path)])
            self.assertTrue(expected_path.exists())

            sample_to_reference = json.loads(
                (Path(tmpdir) / "data" / "sample_to_reference.json").read_text(encoding="utf-8")
            )
            self.assertEqual(sample_to_reference, {"sample1": "reference"})

    def test_file_references_to_json_downloads_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\treference_name\treference_path\n"
                "sample1\thost\tGCF_000001405.40\n"
                "sample2\thost\tGCF_000001405.40\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    listing_response("GCF_000001405.40_GRCh38.p14"),
                    FakeResponse(b">ref\nACGT\n"),
                ],
            ) as fake_urlopen:
                file_references_to_json(str(infofile), tmpdir)

            # The listing is fetched once and the FASTA downloaded once, even
            # though both rows name the same accession.
            self.assertEqual(fake_urlopen.call_count, 2)

            reference_to_file = json.loads(
                (Path(tmpdir) / "data" / "reference_to_file.json").read_text(encoding="utf-8")
            )
            expected_path = (
                Path(tmpdir)
                / "data"
                / "references_cache"
                / "host_GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
            )

            self.assertEqual(reference_to_file["host"], str(expected_path))
            self.assertTrue(expected_path.exists())

            sample_to_reference = json.loads(
                (Path(tmpdir) / "data" / "sample_to_reference.json").read_text(encoding="utf-8")
            )
            self.assertEqual(sample_to_reference, {"sample1": "host", "sample2": "host"})

    def test_file_references_to_json_mixes_accessions_paths_and_urls(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            local_reference = Path(tmpdir) / "local.fna"
            local_reference.write_text(">local\nACGT\n", encoding="utf-8")

            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\treference_name\treference_path\n"
                f"sample1\tlocal\t{local_reference}\n"
                "sample2\tremote\thttps://example.org/reference.fna.gz\n"
                "sample3\tncbi\tGCA_000001635\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                side_effect=[
                    FakeResponse(b">ref\nACGT\n"),
                    listing_response("GCA_000001635.9_GRCm39"),
                    FakeResponse(b">mouse\nACGT\n"),
                ],
            ):
                file_references_to_json(str(infofile), tmpdir)

            reference_to_file = json.loads(
                (Path(tmpdir) / "data" / "reference_to_file.json").read_text(encoding="utf-8")
            )
            cache_dir = Path(tmpdir) / "data" / "references_cache"

            self.assertEqual(reference_to_file["local"], str(local_reference.resolve()))
            self.assertEqual(
                reference_to_file["remote"], str(cache_dir / "remote_reference.fna.gz")
            )
            self.assertEqual(
                reference_to_file["ncbi"],
                str(cache_dir / "ncbi_GCA_000001635.9_GRCm39_genomic.fna.gz"),
            )

    def test_file_references_to_json_reports_unresolvable_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\treference_name\treference_path\n"
                "sample1\thost\tGCF_000001405.99\n",
                encoding="utf-8",
            )

            with patch(
                "drakkar.utils.urlopen",
                return_value=listing_response("GCF_000001405.40_GRCh38.p14"),
            ):
                with self.assertRaises(SystemExit):
                    file_references_to_json(str(infofile), tmpdir)

    def test_reference_cli_validation_accepts_accessions_when_allowed(self) -> None:
        self.assertTrue(
            cli_module.validate_path("GCF_000001405.40", "Reference", allow_accession=True)
        )
        self.assertFalse(
            cli_module.validate_path(
                "GCF_000001405.40", "Reference index tarball", allow_url=True
            )
        )


if __name__ == "__main__":
    unittest.main()
