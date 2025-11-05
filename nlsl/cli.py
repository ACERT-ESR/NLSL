#!/usr/bin/env python3
import argparse
import tempfile
import urllib.request
import zipfile
import os
from pathlib import Path

ZIP_URL = "https://github.com/acert-esr/nlsl/archive/refs/heads/master.zip"
DEST_DIR = Path("NLSL_examples")


def download_and_extract_examples():
    fd, tmp = tempfile.mkstemp(suffix=".zip")
    tmp_path = Path(tmp)
    try:
        with os.fdopen(fd, "wb") as f, urllib.request.urlopen(ZIP_URL) as r:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        with zipfile.ZipFile(tmp_path) as z:
            top = next(Path(n).parts[0] for n in z.namelist() if n)
            prefix = f"{top}/examples/"
            DEST_DIR.mkdir(exist_ok=True)
            for info in z.infolist():
                name = info.filename
                if not name.startswith(prefix) or name.endswith("/"):
                    continue
                rel = Path(*Path(name).parts[2:])  # strip "<top>/examples/"
                out_path = DEST_DIR / rel
                out_path.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as src, open(out_path, "wb") as dst:
                    dst.write(src.read())
    finally:
        tmp_path.unlink(missing_ok=True)
    for d in ("NLSL-master", "nlsl-master"):
        p = Path(d)
        if p.is_dir():
            try:
                p.rmdir()
            except OSError:
                pass


def exampledir(argv=None):
    parser = argparse.ArgumentParser(prog="nlsl")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser(
        "exampledir", help="unpack NLSL examples into ./NLSL_examples"
    )
    args = parser.parse_args(argv)
    if args.cmd == "exampledir":
        download_and_extract_examples()
