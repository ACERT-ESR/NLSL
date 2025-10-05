#!/usr/bin/env bash
# Build distribution artifacts and upload them to PyPI using local configuration.
set -euo pipefail

python -m pip install --upgrade pip
python -m pip install --upgrade build twine

rm -rf build dist
python -m build

twine upload dist/*
