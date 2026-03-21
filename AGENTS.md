# Agent Instructions

- When running locally (e.g. within vscode), be sure to run pip and pytest within the `base` environment: `~/base/bin/activate`
- When building the project locally, always run `pip install -e . --no-build-isolation` to ensure the package uses the expected toolchain.
- Before finishing work, run `pytest` to execute the full test suite.
- Do not preserve backwards compatibility when the user requests a refactor; prefer the cleanest API and remove compatibility shims.
- Note that the tests can take several minutes (up to 10 min) to execute.  Do not give up on the tests prematurely.
