import numpy as np
import pyspecdata as psd
import re
from pathlib import Path
from numpy import r_
from pyspecdata.datadir import pyspec_config
import pygmo as pg
import os
from scipy.optimize import minimize

# Register the example directory with pyspecdata if it is not already present.
if not Path(psd.getDATADIR("nlsl_examples")).exists():
    pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(Path.cwd()))

# Heavy objects (nddata + compiled NLSL state) MUST NOT be created at import time,
# because multiprocessing (spawn) will import this module in each worker.
# Instead, we create exactly one NLSL instance per *process* lazily.


def _build_experimental_spectrum():
    """Load and preprocess the experimental spectrum (runs in each process once)."""
    d_local = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
    d_local.set_units("$B_0$", None)
    d_local = d_local.chunk_auto("harmonic")["harmonic", 0]["phase", 0]

    # We need a temporary NLSL instance only to read the expdat max_points limit.
    import nlsl as _nlsl
    n_tmp = _nlsl.nlsl()

    # ☐ TODO -- the max points should be an accessible property or attribute supplied by the
    # class -- we should not need to do the following in our script
    max_points = n_tmp._core.expdat.data.shape[0] // max(n_tmp._core.expdat.nft.shape[0], 1)

    # {{{ we use convolution to downsample the data
    if d_local.data.shape[0] > max_points:
        divisor = d_local.shape["$B_0$"] // max_points + 1
        dB = np.diff(d_local["$B_0$"][r_[0, 1]]).item()
        d_orig_max = d_local.data.max()
        d_local.convolve("$B_0$", dB / 6 * divisor)
        d_local = d_local["$B_0$", 0::divisor]
        d_local *= d_orig_max / d_local.data.max()
    # }}}

    # Normalize experimental data (keep your existing convention)
    d_local.data /= d_local.data.max() - d_local.data.min()

    return d_local


# Provide reasonable starting parameters
initial_params = {
    "in2": np.int32(2),
    "gxx": np.float64(2.0089),
    "gyy": np.float64(2.0058),
    "gzz": np.float64(2.0021),
    "axx": np.float64(4.9),
    "ayy": np.float64(4.9),
    "azz": np.float64(33.0),
    "lemx": np.int32(6),
    "lomx": np.int32(5),
    "kmx": np.int32(4),
    "mmx": np.int32(4),
    "ipnmx": np.int32(2),
    "rbar": np.float64(7.935044819885658),
    "n": np.float64(0.6932829362508512),
    "c20": np.float64(2.5061887785331742),
    "c22": np.float64(-1.2608085913642202),
    "betad": np.float64(29.942804557290664),
    "gib0": np.float64(2.0101124440329796),
    "gib2": np.float64(0.5050889997199121),
}

param_tokens = sorted(
    set(initial_params.keys())
    - set(["in2", "kmx", "mmx", "lemx", "lomx", "ipnmx", "b0"])
)
print(param_tokens)

# -------------------------
# Build bounds for all optimizable scalar parameters
# (computed from initial_params only; does not touch NLSL at import time)
# -------------------------

bounds = []
for k in param_tokens:
    v = float(initial_params[k])
    if k in ["gxx", "gyy", "gzz"]:
        # g-tensor is really a difference about 2
        bounds.append((((v - 2) * 0.6) + 2, ((v - 2) * 1.4) + 2))
    elif k == "betad":
        bounds.append((0.0, 180.0))
    else:
        # generic: allow roughly one order of magnitude around current value,
        # handling negative values safely
        low, high = v * 0.1, v * 2.0
        if low > high:
            low, high = high, low
        bounds.append((low, high))


# -------------------------
# Per-process NLSL + nddata cache (ONE instance per process)
# -------------------------

_WORKER = {"pid": None, "n": None, "d": None}


def _get_ctx():
    pid = os.getpid()
    if _WORKER["pid"] == pid and _WORKER["n"] is not None and _WORKER["d"] is not None:
        return _WORKER

    d_local = _build_experimental_spectrum()

    import nlsl as _nlsl
    n_local = _nlsl.nlsl()
    n_local.update(initial_params)
    n_local.load_nddata(d_local)

    _WORKER.update({"pid": pid, "n": n_local, "d": d_local})
    return _WORKER

def objective(x_vec):
    ctx = _get_ctx()
    n_local = ctx["n"]
    for value, token in zip(x_vec, param_tokens):
        n_local[token] = value
    site_spectra = n_local.current_spectrum
    simulated_total = np.squeeze(n_local.weights @ site_spectra)
    simulated_total /= simulated_total.max() - simulated_total.min()
    return simulated_total


def residual_norm(x_vec):
    ctx = _get_ctx()
    d_local = ctx["d"]
    sim = objective(x_vec)
    return np.linalg.norm(d_local.data - sim)


# Target residual (L2 norm) for guidance / restarts
target_residual = 5e-2


class NLSLFitness:
    """Pagmo problem wrapper around the NLSL spectrum fit."""

    def __init__(self, param_tokens, bounds):
        self.param_tokens = list(param_tokens)
        self.bounds = list(bounds)

    def fitness(self, x):
        return [residual_norm(x)]

    def get_bounds(self):
        lb = [b[0] for b in self.bounds]
        ub = [b[1] for b in self.bounds]
        return (lb, ub)


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    # Build pagmo problem and archipelago (islands + migration)
    prob = pg.problem(NLSLFitness(param_tokens, bounds))

    pop_size = 40
    num_islands = 8

    algo = pg.algorithm(pg.de1220(gen=20))
    algo.set_verbosity(1)

    archi = pg.archipelago(
        n=num_islands,
        algo=algo,
        prob=prob,
        pop_size=pop_size,
        udi=pg.mp_island(),
    )

    max_epochs = 50

    for epoch in range(max_epochs):
        archi.evolve()
        archi.wait_check()

        # Gather full population across all islands to measure spread
        all_pop = [isl.get_population().get_x() for isl in archi]
        all_x = np.vstack(all_pop)
        mean_vec = np.mean(all_x, axis=0)
        std_vec = np.std(all_x, axis=0)
        spread = float(np.mean(std_vec / (np.abs(mean_vec) + 1e-12)))

        champs_f = archi.get_champions_f()
        best_resid = float(min(f[0] for f in champs_f))

        print(f"epoch {epoch}: best_resid={best_resid:.5g}, spread={spread:.3g}")

        # Diversity reinjection if population has collapsed but fit is still poor
        if spread < 1e-3 and best_resid > target_residual:
            for isl in archi:
                if np.random.rand() < 0.3:
                    isl.set_population(pg.population(prob, pop_size))


    # Extract best solution from the archipelago
    champ_f = archi.get_champions_f()
    champ_x = archi.get_champions_x()
    idx_best = np.argmin([f[0] for f in champ_f])
    best_x = np.array(champ_x[idx_best], dtype=float)

    # -------------------------
    # Local refinement with Powell (derivative-free), inspired by
    # GA+Powell hybrids used successfully for EPR spectral fitting.
    # -------------------------

    lb = np.array([b[0] for b in bounds], dtype=float)
    ub = np.array([b[1] for b in bounds], dtype=float)


    def residual_for_scipy(x_vec):
        # Project onto bounds to be robust against Powell overshooting
        x_vec = np.clip(x_vec, lb, ub)
        sim = objective(x_vec)
        ctx = _get_ctx()
        d_local = ctx["d"]
        res = d_local.data - sim
        return float(np.linalg.norm(res))


    powell_result = minimize(
        residual_for_scipy,
        best_x,
        method="Powell",
        options={"maxiter": 2000, "xtol": 1e-4, "ftol": 1e-4},
    )

    best_x_refined = np.clip(powell_result.x, lb, ub)

    # Apply optimized parameters in the main process context
    ctx = _get_ctx()
    n_local = ctx["n"]
    d_local = ctx["d"]

    for value, token in zip(best_x_refined, param_tokens):
        n_local[token] = value

    simulated_total = objective(best_x_refined)
    residual = d_local.data - simulated_total
    with psd.figlist_var() as fl:
        fl.next("Global+local DE/Powell fit on pyspecdata trace")
        fl.plot(d_local, alpha=0.6, label="experimental")

        field_axis = np.asarray(d_local[d_local.dimlabels[0]], dtype=float)
        fl.plot(field_axis, simulated_total, alpha=0.9, label="global+local fit")
        fl.plot(field_axis, residual, alpha=0.8, label="residual")
    print("final params", n_local.items())
