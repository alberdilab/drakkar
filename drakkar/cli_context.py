from pathlib import Path

import yaml

PACKAGE_DIR = Path(__file__).parent
CONFIG_PATH = PACKAGE_DIR / "workflow" / "config.yaml"


def load_config():
    """Load fixed variables from config.yaml."""
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH, "r") as f:
            return yaml.safe_load(f)
    return {}


config_vars = load_config()

HEADER1 = "\033[1;95m"
ERROR = "\033[1;31m"
INFO = "\033[1;34m"
RESET = "\033[0m"

WORKFLOW_RUN_COMMANDS = {
    "complete",
    "preprocessing",
    "cataloging",
    "profiling",
    "dereplicating",
    "annotating",
    "inspecting",
    "expressing",
    "database",
    "environments",
}

READ_ONLY_COMMANDS = {"config", "logging"}
CATALOGING_BINNER_ORDER = ("metabat", "maxbin", "semibin", "comebin")
DEFAULT_CATALOGING_BINNERS = ",".join(CATALOGING_BINNER_ORDER)
CATALOGING_BINNER_ALIASES = {
    "metabat": "metabat",
    "metabat2": "metabat",
    "maxbin": "maxbin",
    "maxbin2": "maxbin",
    "semibin": "semibin",
    "semibin2": "semibin",
    "comebin": "comebin",
}
