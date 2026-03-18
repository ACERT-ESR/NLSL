# Agent Instructions

- When building the project locally, always run `pip install -e . --no-build-isolation` to ensure the package uses the expected toolchain.
- Before finishing work, run `pytest` to execute the full test suite.

- Do not preserve backwards compatibility when the user requests a refactor; prefer the cleanest API and remove compatibility shims.
