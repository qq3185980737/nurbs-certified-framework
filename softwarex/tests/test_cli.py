import json
from pathlib import Path

from ssi_cert import cli


def test_accept(cert_path, capsys):
    rc = cli.main(["verify", str(cert_path), "--quiet"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "ACCEPT" in out


def test_reject(cert_path, tmp_path, capsys):
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    # A root box outside the parameter domain must not be silently accepted.
    cert["tree"]["box"][1] = "0x1.0000000000000p+1"
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(cert), encoding="utf-8")
    rc = cli.main(["verify", str(bad), "--quiet"])
    assert rc == 1
    assert "REJECT" in capsys.readouterr().out


def test_malformed_input(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    rc = cli.main(["verify", str(bad)])
    assert rc == 2
    assert "ERROR" in capsys.readouterr().err


def test_info(cert_path, capsys):
    rc = cli.main(["info", str(cert_path)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "sha256:" in out
    assert "keys:" in out


def test_verify_dir(cert_path, capsys):
    rc = cli.main(["verify-dir", str(cert_path.parent), "--quiet"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "SUMMARY" in out
    assert "reject=0" in out


def test_chunked():
    root = Path(__file__).resolve().parents[1]
    path = root / "examples" / "certs" / "cert_chunked.json"
    rc = cli.main(["verify", str(path), "--chunked", "--jobs", "1", "--quiet"])
    assert rc == 0
