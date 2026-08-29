import subprocess
import sys

try:
    from importlib.metadata import PackageNotFoundError, version as get_distribution_version
except ImportError:  # pragma: no cover - Python < 3.8 fallback
    try:
        from importlib_metadata import PackageNotFoundError, version as get_distribution_version
    except ImportError:  # pragma: no cover - fallback if backport is absent
        PackageNotFoundError = Exception
        get_distribution_version = None

from drakkar import __version__
from drakkar.config_persistence import (
    print_reconcile_report,
    reconcile_config,
    snapshot_config,
)
from drakkar.display import display_update_success
from drakkar.output import print

def get_installed_drakkar_version():
    if get_distribution_version is None:
        return __version__
    try:
        return get_distribution_version("drakkar")
    except PackageNotFoundError:
        return __version__
    except Exception:
        return __version__

def run_update(skip_deps=False):
    # workflow/config.yaml ships inside the package, so the reinstall below
    # replaces it with the defaults of the new release. Save the site-specific
    # values first, then write them back once the new config is in place.
    try:
        saved = snapshot_config()
    except OSError as exc:
        print(f"Current config.yaml values could not be saved: {exc}", file=sys.stderr, flush=True)
        saved = None
    pip_cmd = [
        sys.executable, "-m", "pip", "install",
        "--upgrade", "--force-reinstall",
    ]
    if skip_deps:
        pip_cmd.append("--no-deps")
    pip_cmd.append("git+https://github.com/alberdilab/drakkar.git")
    try:
        update_result = subprocess.run(pip_cmd)
    except Exception as exc:
        print(f"Update failed: {exc}", file=sys.stderr, flush=True)
        return 1
    if update_result.returncode != 0:
        return update_result.returncode
    if saved is not None:
        restore_config_values(saved)
    display_update_success(get_installed_drakkar_version())
    return 0


def restore_config_values(saved):
    """Put the pre-update database values back into the freshly installed config."""
    try:
        resolutions = reconcile_config(saved=saved)
    except OSError as exc:
        print(f"Database paths could not be restored into config.yaml: {exc}", file=sys.stderr, flush=True)
        print("Run 'drakkar config --restore' once the file is writable.", file=sys.stderr, flush=True)
        return
    print_reconcile_report(resolutions)
