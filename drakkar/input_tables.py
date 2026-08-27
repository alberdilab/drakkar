import pandas as pd

from drakkar.input_errors import report_input_resolution_errors
from drakkar.output import print

# Candidate delimiters, in priority order. Tab wins ties so that historical
# tab-separated tables keep behaving exactly as before.
INPUT_TABLE_DELIMITERS = ("\t", ",", ";")

DELIMITER_LABELS = {"\t": "tab-separated", ",": "comma-separated", ";": "semicolon-separated"}

# Each table is read several times while the run is prepared; the delimiter is
# announced only once per file to keep the launch output readable.
_ANNOUNCED_DELIMITERS = set()

def detect_table_delimiter(table_path):
    """Guess the delimiter of a user-provided table from its header line.

    Returns a tab when the file is empty or the header contains none of the
    supported delimiters, so single-column tables keep working.
    """
    header = ""
    try:
        with open(table_path, "r", encoding="utf-8-sig") as handle:
            for line in handle:
                if line.strip():
                    header = line
                    break
    except OSError:
        return "\t"

    counts = {delimiter: header.count(delimiter) for delimiter in INPUT_TABLE_DELIMITERS}
    delimiter = max(INPUT_TABLE_DELIMITERS, key=lambda candidate: counts[candidate])
    return delimiter if counts[delimiter] else "\t"

def read_input_table(table_path, label="input table", dtype=None):
    """Read a user-provided table, accepting tab, comma, or semicolon delimiters.

    Column names are stripped of surrounding whitespace, and unreadable files
    are reported as input errors instead of raising a pandas traceback.
    """
    delimiter = detect_table_delimiter(table_path)
    if delimiter != "\t" and str(table_path) not in _ANNOUNCED_DELIMITERS:
        _ANNOUNCED_DELIMITERS.add(str(table_path))
        print(f"Reading {label} {table_path} as a {DELIMITER_LABELS[delimiter]} table.")

    try:
        df = pd.read_csv(
            table_path, sep=delimiter, encoding="utf-8-sig", dtype=dtype
        )
    except Exception as exc:
        report_input_resolution_errors([f"Could not read the {label} {table_path}: {exc}"])

    df.columns = [str(column).strip() for column in df.columns]
    return df
