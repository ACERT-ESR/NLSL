#!/usr/bin/env bash
# Build distribution artifacts and upload them to PyPI using local configuration.
set -euo pipefail

python -m pip install --upgrade pip
python -m pip install --upgrade build twine auditwheel patchelf

rm -rf build dist
python -m build

# Repair linux wheels to comply with PyPI's manylinux policy.
for whl in dist/*.whl; do
    if [[ "${whl}" == *linux_* ]]; then
        python -m auditwheel repair "${whl}" -w dist
        rm "${whl}"
    fi
done

twine upload dist/*
