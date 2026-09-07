"""Split a flat FASTQ holding both mates of a paired run into R1 and R2.

Some runs archived as PAIRED are published as a single ``<run>.fastq.gz``
instead of a ``_1``/``_2`` pair. When a submitter uploads reads that were
already quality-trimmed the two files no longer line up read-for-read, so the
SRA loader gives up on pairing them and stores every read as its own
single-read spot: all of one file's reads, then all of the other's. ENA mirrors
that object and publishes no per-mate URL, and NCBI keeps the original
submitted files in a requester-pays bucket, so the only place the mates can be
recovered is here, after the download.

Every read still carries the read index it was loaded under as a ``/1`` or
``/2`` suffix in its name, and records are routed on that index rather than on
position: which half comes first varies from run to run.

The split is a single streaming pass. Nothing is buffered in memory beyond one
record and nothing is staged uncompressed on disk, because these files run to
2-4 GB compressed each and a study can carry a hundred of them.
"""

import gzip
import os
import re
import shutil

from drakkar.input_errors import FastqSplitError

# Level 9 roughly triples the time spent compressing for a few percent of size,
# which is a poor trade on files this large.
DEFAULT_COMPRESS_LEVEL = 6

GZIP_MAGIC = b"\x1f\x8b"

READ_BUFFER_BYTES = 4 * 1024 * 1024

# The mate suffix as ENA writes it, on the tail of a whitespace field.
READ_INDEX_SUFFIX_PATTERN = re.compile(rb"/([12])$")

# Trailing spot number in a read name, e.g. the ".1" of "SRR9851002.1".
SPOT_NUMBER_SUFFIX_PATTERN = re.compile(rb"\.\d+$")


def configured_tmp_dir():
    """Directory for in-progress split output, or None to write beside the outputs.

    ``DRAKKAR_TMPDIR`` redirects the partial files off the output filesystem,
    which matters when a study's worth of unsplit runs would not otherwise fit
    alongside the finished halves.
    """
    value = (os.environ.get("DRAKKAR_TMPDIR") or "").strip()
    return value or None


def read_index(header):
    """Return 1 or 2 for a FASTQ header line, or None when it carries no mate suffix.

    ENA writes the index as the tail of the last whitespace field
    (``@SRR9851002.1 1/1``); older Illumina and BGI names carry it on the read
    name itself (``@READ/1``). Both are read, in that order.
    """
    fields = header.split()
    if not fields:
        return None
    for field in (fields[-1], fields[0]):
        match = READ_INDEX_SUFFIX_PATTERN.search(field)
        if match:
            return int(match.group(1))
    return None


def record_stem(header):
    """Base name shared by every read of a run, e.g. b'SRR9851002'.

    Both the mate suffix and the per-read spot number are stripped, so the stem
    is what remains once everything that varies between reads is removed.
    """
    fields = header.split()
    if not fields:
        return None
    name = fields[0].lstrip(b"@")
    name = READ_INDEX_SUFFIX_PATTERN.sub(b"", name)
    name = SPOT_NUMBER_SUFFIX_PATTERN.sub(b"", name)
    return name or None


def _looks_gzipped(path):
    with open(path, "rb") as handle:
        return handle.read(2) == GZIP_MAGIC


def _open_source(path):
    if _looks_gzipped(path):
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=READ_BUFFER_BYTES)


def _tmp_path(dest_path, tmp_dir):
    if tmp_dir:
        os.makedirs(tmp_dir, exist_ok=True)
        return os.path.join(tmp_dir, f"{os.path.basename(dest_path)}.tmp")
    return f"{dest_path}.tmp"


def _promote(tmp_path, dest_path):
    try:
        os.replace(tmp_path, dest_path)
    except OSError:
        # A configured DRAKKAR_TMPDIR may sit on another filesystem, where
        # rename cannot work and the bytes have to be copied across.
        shutil.move(tmp_path, dest_path)


def _remove_quietly(paths):
    for path in paths:
        try:
            if path and os.path.isfile(path):
                os.remove(path)
        except OSError:
            pass


def _terminated(line):
    """The line with exactly one trailing newline and no stray carriage return.

    A CRLF-terminated source would otherwise carry the ``\r`` into the sequence
    and quality strings, where it is not a legal character.
    """
    return line.rstrip(b"\r\n") + b"\n"


def split_unsplit_fastq(
    source_path,
    read1_path,
    read2_path,
    tmp_dir=None,
    compresslevel=DEFAULT_COMPRESS_LEVEL,
    fallback_stem=None,
):
    """Stream ``source_path`` into gzipped ``read1_path`` and ``read2_path``.

    Returns the number of records written to each output. Raises
    ``FastqSplitError``, leaving no output behind, when the file does not split
    into two equal, non-empty halves -- a mis-split silently produces a garbage
    library, so it must never be handed downstream.

    Records are renumbered per half, so mate *k* of both outputs shares the name
    ``<stem>.<k>``. The archive numbers every read of an unsplit run with its
    own spot, which leaves the two halves with no name in common and is useless
    as a pairing key; renumbering reproduces the naming an ordinary ENA
    ``_1``/``_2`` pair would have carried, which is what the tools downstream
    pair on.
    """
    tmp1 = _tmp_path(read1_path, tmp_dir)
    tmp2 = _tmp_path(read2_path, tmp_dir)
    counts = {1: 0, 2: 0}
    stem = None

    try:
        _remove_quietly([tmp1, tmp2])
        with _open_source(source_path) as source, gzip.open(
            tmp1, "wb", compresslevel=compresslevel
        ) as out1, gzip.open(tmp2, "wb", compresslevel=compresslevel) as out2:
            handles = {1: out1, 2: out2}
            record_number = 0
            while True:
                header = source.readline()
                if not header:
                    break
                sequence = source.readline()
                separator = source.readline()
                quality = source.readline()
                record_number += 1

                if not quality:
                    raise FastqSplitError(
                        f"record {record_number} is truncated: a FASTQ record needs four lines."
                    )
                header = header.rstrip(b"\r\n")
                if not header.startswith(b"@") or not separator.startswith(b"+"):
                    raise FastqSplitError(
                        f"record {record_number} is not a well-formed FASTQ record."
                    )

                index = read_index(header)
                if index is None:
                    raise FastqSplitError(
                        f"record {record_number} carries no /1 or /2 read index in its name "
                        f"({header.decode('utf-8', 'replace')}), so its mate cannot be determined."
                    )

                if stem is None:
                    stem = record_stem(header) or _encoded_fallback_stem(fallback_stem)

                counts[index] += 1
                position = counts[index]
                handle = handles[index]
                handle.write(b"@%s.%d %d/%d\n" % (stem, position, position, index))
                handle.write(_terminated(sequence))
                handle.write(b"+\n")
                handle.write(_terminated(quality))

        _validate_counts(counts, source_path)
        _promote(tmp1, read1_path)
        _promote(tmp2, read2_path)
    except BaseException:
        _remove_quietly([tmp1, tmp2])
        raise

    return counts[1]


def _encoded_fallback_stem(fallback_stem):
    return (fallback_stem or "read").encode("utf-8", "replace")


def _validate_counts(counts, source_path):
    if not counts[1] and not counts[2]:
        raise FastqSplitError(f"{source_path} contains no FASTQ records.")
    for index, other in ((1, 2), (2, 1)):
        if not counts[index]:
            raise FastqSplitError(
                f"every one of the {counts[other]} records in {source_path} carries read index "
                f"{other}, so the file holds only one mate and cannot be split into a pair."
            )
    if counts[1] != counts[2]:
        raise FastqSplitError(
            f"{source_path} splits into {counts[1]} forward and {counts[2]} reverse records. "
            "The two mates must contain the same number of reads."
        )
