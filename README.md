# NLSL

The NLSL project provides a Python interface to the ACERT nonlinear least
squares (NLSL) spectroscopic fitting toolkit, bundling the original Fortran
core with a modern Python API for building, fitting, and analysing spin
resonance models.【F:pyproject.toml†L11-L24】【F:nlsl/__init__.py†L1-L104】

## Key features

- **Pythonic parameter access** – the :class:`nlsl.nlsl` class exposes the
  complete parameter set from the legacy NLSL engine through mapping-style
  accessors, making it easy to read or update magnetic, diffusion, and fitting
  options from Python code.【F:nlsl/__init__.py†L60-L170】
- **Programmable fitting control** – the :class:`nlsl.fit_params` helper keeps
  the low-level Levenberg–Marquardt settings in sync with the Fortran
  structures, so tolerances, iteration limits, and other solver options can be
  modified directly without crafting text run files.【F:nlsl/__init__.py†L15-L58】
- **Runfile compatibility** – existing NLSL ``.run`` scripts can be replayed via
  :meth:`nlsl.nlsl.procline`, preserving the familiar workflow while allowing
  incremental migration to the Python API.【F:nlsl/__init__.py†L82-L109】
- **Spectrum capture** – helper methods such as
  :meth:`nlsl.nlsl.current_spectrum` and :meth:`nlsl.nlsl.write_spc` let you
  evaluate the current model or persist generated spectra for downstream
  analysis and visualisation.【F:nlsl/__init__.py†L105-L158】

## Installation

Clone the repository and install the package in editable mode so the compiled
Fortran components stay in sync with your working tree:

```bash
pip install -e . --no-build-isolation
```

## Usage overview

Instantiate :class:`nlsl.nlsl` to work with parameters programmatically, or
stream traditional runfiles back into the engine.

```python
import nlsl

n = nlsl.nlsl()
n["nsite"] = 1
n.fit_params["maxitr"] = 40
n.procline("data sampl1 ascii nspline 200 bc 20 shift")
site_spectra = n.fit()
total_spectrum = n["weights"] @ site_spectra
```

The mapping interface mirrors the parameter names defined by the original
Fortran code (see ``nlshlp.txt`` for a full command reference), while
:meth:`nlsl.nlsl.fit` runs the nonlinear least-squares optimiser and returns the
latest site spectra, while the companion weight matrix stays accessible through
the mapping interface.【F:nlsl/__init__.py†L60-L204】

## Examples

A curated set of runnable examples is included under ``examples/`` to help you
get started with real datasets and typical NLSL workflows:

### Download the examples directory

You can download the full repository ZIP directly from GitHub and extract only the `examples/` folder:

[⬇️ Download ZIP (GitHub)](https://github.com/acert-esr/nlsl/archive/refs/heads/master.zip)

```bash
unzip NLSL-master.zip "NLSL-master/examples/*"
```

### Description of examples

- **Runfile suite (`sampl1`–`sampl5`)** – canonical `.run` scripts paired with
  `.dat` inputs that demonstrate anisotropic rotations, multi-spectrum fits,
  and other standard analyses. Run them with the classic command-line
  interface and compare the generated `.spc` and `.log` outputs to the provided
  reference files to verify your setup.【F:examples/README.rst†L1-L22】
- **`runexample.py`** – command-line helper that executes one of the numbered
  runfile scenarios, cascades any nested `call` directives, reports relative
  RMS errors, and plots the resulting spectra with experimental overlays for
  quick diagnostics.【F:examples/runexample.py†L1-L93】
- **`runexample_first.py`** – simplified variant that focuses on the first
  sample run, computing a global relative RMS metric while rendering the fitted
  spectrum.【F:examples/runexample_first.py†L1-L74】
- **`break_down_example.py`** – step-by-step walkthrough of the first sample
  showcasing programmatic parameter updates, basis set configuration, and fit
  execution directly from Python without relying on runfiles.【F:examples/break_down_example.py†L1-L83】
- **`BL05CHSA.py`** – reproduces the `BL05CHSA.run` workflow, highlighting how
  to plot individual spectral components alongside the overall fit and
  cumulative integrals.【F:examples/BL05CHSA.py†L1-L40】
- **`c16pc371e.py`** – demonstrates interacting with a custom runfile while
  inspecting the solver configuration before and after execution and plotting
  the normalised spectra.【F:examples/c16pc371e.py†L1-L46】

Each script expects to be run from within the `examples/` directory so that the
associated data and reference outputs are discovered automatically.【F:examples/runexample.py†L22-L41】

## Additional resources

- The ``examples/README.rst`` file documents the legacy testing workflow for
  the runfile suite and points to the reference outputs for verification.
- ``nlshlp.txt`` lists every command available to the original interpreter along
  with detailed descriptions of their arguments and expected effects.
