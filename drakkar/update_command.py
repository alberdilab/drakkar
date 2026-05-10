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
    display_update_success(get_installed_drakkar_version())
    return 0
