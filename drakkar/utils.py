import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
import yaml
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode, urlparse
from urllib.request import urlopen

import pandas as pd

from drakkar import __version__
from drakkar.ascii import (
    DRAKKAR_LOGO_ART,
    DRAKKAR_SHIP_ART,
    END_ART,
    INTRO_ROWS,
    UNLOCK_ART,
    UPDATE_SUCCESS_ASCII_TEMPLATE,
)
from drakkar.output import Text as RichText, print, prompt

from drakkar import display as _display
from drakkar import downloads as _downloads
from drakkar import input_errors as _input_errors
from drakkar import input_manifests as _input_manifests
from drakkar import system_checks as _system_checks

from drakkar.display import (
    BANNER_DELAY_SECONDS,
    DRAKKAR_INTRO_STYLE,
    DRAKKAR_LOGO_STYLE,
    DRAKKAR_SHIP_STYLE,
    DRAKKAR_VERSION_BADGE_STYLE,
    _ascii_block_width,
    _ascii_block_with_bottom_right_badge,
    _center_block,
    _format_update_success_ascii,
    _intro_box,
    _styled_ascii_art,
    _version_badge_lines,
    get_drakkar_banner_blocks,
    get_drakkar_banner_renderables,
)
from drakkar.downloads import (
    DEFAULT_DOWNLOAD_RETRIES,
    READ1_BASENAME_PATTERN,
    READ2_BASENAME_PATTERN,
    REMOTE_URL_SCHEMES,
    _download_once,
    _download_sftp,
    _download_urlopen,
    _has_value,
    _normalize_ena_fastq_url,
    _normalized_value,
    _resolve_reads_path,
    _retry_delay,
    _split_paired_fastq_urls,
    _validate_paired_read_maps,
)
from drakkar.input_errors import DownloadError, InputFileError, require_non_empty_file
from drakkar.input_manifests import ASSEMBLY_COLUMN_CANDIDATES

_ORIGINAL_BANNER_ANIMATION_ENABLED = _display._banner_animation_enabled


def _sync_download_dependencies():
    _downloads.urlopen = urlopen
    _downloads.time = time
    _downloads.print = print
    _input_errors.print = print


def _sync_display_dependencies():
    _display.print = print
    _display.time = time
    banner_animation_enabled = globals()["_banner_animation_enabled"]
    if (
        getattr(banner_animation_enabled, "__module__", None) == __name__
        and getattr(banner_animation_enabled, "__name__", None) == "_banner_animation_enabled"
    ):
        banner_animation_enabled = _ORIGINAL_BANNER_ANIMATION_ENABLED
    _display._banner_animation_enabled = banner_animation_enabled


def _sync_system_dependencies():
    _system_checks.print = print
    _system_checks.prompt = prompt


def report_input_resolution_errors(errors):
    _input_errors.print = print
    return _input_errors.report_input_resolution_errors(errors)


def is_url(value):
    return _downloads.is_url(value)


def download_to_cache(*args, **kwargs):
    _sync_download_dependencies()
    return _downloads.download_to_cache(*args, **kwargs)


def resolve_input_manifest(*args, **kwargs):
    _sync_download_dependencies()
    return _downloads.resolve_input_manifest(*args, **kwargs)


def resolve_accession_to_reads(*args, **kwargs):
    _sync_download_dependencies()
    return _downloads.resolve_accession_to_reads(*args, **kwargs)


def resolve_sample_read_lists(*args, **kwargs):
    _sync_download_dependencies()
    return _downloads.resolve_sample_read_lists(*args, **kwargs)


def resolve_preprocessed_read_lists(*args, **kwargs):
    _sync_download_dependencies()
    return _downloads.resolve_preprocessed_read_lists(*args, **kwargs)


def _banner_animation_enabled(*args, **kwargs):
    return _display._banner_animation_enabled(*args, **kwargs)


def display_banner_sequence(*args, **kwargs):
    _sync_display_dependencies()
    return _display.display_banner_sequence(*args, **kwargs)


def display_drakkar(*args, **kwargs):
    _sync_display_dependencies()
    return _display.display_drakkar(*args, **kwargs)


def display_unlock(*args, **kwargs):
    _sync_display_dependencies()
    return _display.display_unlock(*args, **kwargs)


def display_end(*args, **kwargs):
    _sync_display_dependencies()
    return _display.display_end(*args, **kwargs)


def display_update_success(*args, **kwargs):
    _sync_display_dependencies()
    return _display.display_update_success(*args, **kwargs)


def is_snakemake_locked(*args, **kwargs):
    return _system_checks.is_snakemake_locked(*args, **kwargs)


def check_screen_session(*args, **kwargs):
    _sync_system_dependencies()
    return _system_checks.check_screen_session(*args, **kwargs)


def check_reference_columns(*args, **kwargs):
    return _input_manifests.check_reference_columns(*args, **kwargs)


def get_assembly_column_name(*args, **kwargs):
    return _input_manifests.get_assembly_column_name(*args, **kwargs)


def check_assembly_column(*args, **kwargs):
    return _input_manifests.check_assembly_column(*args, **kwargs)


def file_samples_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_samples_to_json(*args, **kwargs)


def file_references_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_references_to_json(*args, **kwargs)


def argument_samples_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.argument_samples_to_json(*args, **kwargs)


def argument_references_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.argument_references_to_json(*args, **kwargs)


def file_preprocessed_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_preprocessed_to_json(*args, **kwargs)


def argument_preprocessed_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.argument_preprocessed_to_json(*args, **kwargs)


def file_assemblies_to_json(*args, **kwargs):
    return _input_manifests.file_assemblies_to_json(*args, **kwargs)


def file_coverages_to_json(*args, **kwargs):
    return _input_manifests.file_coverages_to_json(*args, **kwargs)


def file_bins_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_bins_to_json(*args, **kwargs)


def path_bins_to_json(*args, **kwargs):
    return _input_manifests.path_bins_to_json(*args, **kwargs)


def file_mags_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_mags_to_json(*args, **kwargs)


def path_mags_to_json(*args, **kwargs):
    return _input_manifests.path_mags_to_json(*args, **kwargs)


def microdiversity_selection_to_json(*args, **kwargs):
    return _input_manifests.microdiversity_selection_to_json(*args, **kwargs)


def preprocessing_summary(*args, **kwargs):
    return _input_manifests.preprocessing_summary(*args, **kwargs)


def file_transcriptome_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.file_transcriptome_to_json(*args, **kwargs)


def argument_transcriptome_to_json(*args, **kwargs):
    _sync_download_dependencies()
    return _input_manifests.argument_transcriptome_to_json(*args, **kwargs)
