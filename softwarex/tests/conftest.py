import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


@pytest.fixture
def cert_path() -> Path:
    return ROOT / "examples" / "certs" / "cert_demo.json"
