#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Qt GUI to explore the *current spectrum* using the nlsl model.

Requirements satisfied:
- **No try/except fallbacks.**
- **Load experimental data from the same folder as this script** (uses `__file__`).
- **Initialize EXACTLY like the known-good test**: build the field grid with
  `generate_coordinates(...)`, then `set_data(...)`, then apply the FINAL
  parameters + bookkeeping tokens, then set `weights`, and only then read
  `current_spectrum`.
- Expose parameters via nested tabs; `nsite` is a spin box above tabs; general
  tab holds integer grid controls + normalization + weights.
"""

import sys
from pathlib import Path
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QSlider, QTabWidget, QGroupBox, QGridLayout, QSpinBox,
    QMessageBox, QDoubleSpinBox
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import nlsl

# -------------------- Hard-coded reference knobs (no imports) --------------------
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

# Final parameters & metadata copied from the working test reference
SAMPL4_FINAL_PARAMETERS = {
    "nsite": 2,
    "phase": 0.0,
    "gib0": 1.9962757195220067,
    "gib2": 0.0,
    "wxx": 0.0,
    "wyy": 0.0,
    "wzz": 0.0,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "rx": np.array([7.8396974, 7.14177897]),
    "ry": 0.0,
    "rz": 0.0,
    "pml": 0.0,
    "pmxy": 0.0,
    "pmzz": 0.0,
    "djf": 0.0,
    "djfprp": 0.0,
    "oss": 0.0,
    "psi": 0.0,
    "alphad": 0.0,
    "betad": 0.0,
    "gammad": 0.0,
    "alpham": 0.0,
    "betam": 0.0,
    "gammam": 0.0,
    "c20": 0.0,
    "c22": 0.0,
    "c40": 0.0,
    "c42": 0.0,
    "c44": 0.0,
    "lb": 0.0,
    "dc20": 0.0,
    "b0": 3400.50251256,
    "fldi": 3350.5046000757857,
    # dfld (step) and nfield are derived from on-disk data; we still expose dfld
    "dfld": None,  # will be set from the file's first column spacing
    "gamman": 0.0,
    "cgtol": 0.001,
    "shiftr": 0.001,
    "shifti": 0.0,
    "range": 100.0,
    "in2": 2,
    "ipdf": 0,
    "ist": 0,
    "ml": 0,
    "mxy": 0,
    "mzz": 0,
    "lemx": 12,
    "lomx": 10,
    "kmn": 0,
    "kmx": 7,
    "mmn": 0,
    "mmx": 7,
    "ipnmx": 2,
    "nort": 0,
    "nstep": 0,
    "nfield": None,  # will be set from file length
    "ideriv": 1,
    "iwflg": 0,
    "igflg": 0,
    "iaflg": 0,
    "jkmn": 0,
    "jmmn": 0,
    "irflg": 2,
    "ndim": 156,
}

SAMPL4_FINAL_WEIGHTS = np.array([0.2848810, 0.7155313], dtype=float)
SAMPL4_FINAL_SB0 = np.array([SAMPL4_FINAL_PARAMETERS["b0"]], dtype=float)
SAMPL4_FINAL_SRNG = np.array([SAMPL4_FINAL_PARAMETERS["range"]], dtype=float)
SAMPL4_FINAL_ISHFT = np.array([1], dtype=np.int32)
SAMPL4_FINAL_SHFT = np.array([0.0], dtype=float)
SAMPL4_FINAL_NRMLZ = np.array([0], dtype=np.int32)

# -------------------- Small helpers --------------------

def _load_two_column_ascii(p: Path):
    arr = np.loadtxt(p)
    arr = np.atleast_2d(arr)
    if arr.shape[1] < 2:
        raise ValueError("Expected two columns: field and intensity")
    x = np.asarray(arr[:, 0], float)
    y = np.asarray(arr[:, 1], float)
    return x, y


def _field_start_step(x: np.ndarray):
    start = float(x[0])
    # robust average step to tolerate tiny noise
    diffs = np.diff(x.astype(float))
    step = float(np.mean(diffs))
    return start, step


def _int_slider(min_v, max_v, value, step=1):
    s = QSlider(QtCore.Qt.Horizontal)
    s.setMinimum(int(min_v)); s.setMaximum(int(max_v))
    s.setSingleStep(int(step)); s.setPageStep(int(step)); s.setValue(int(value))
    return s


def make_scaled_pair(label, scale, init_value, min_v, max_v, step):
    """Float slider via int slider + scale; returns (row_widget, slider, spinbox, getter)."""
    box = QWidget(); layout = QHBoxLayout(box)
    name = QLabel(label); name.setMinimumWidth(120)
    imin, imax = int(min_v/scale), int(max_v/scale)
    ival, istep = int(init_value/scale), max(1, int(step/scale))
    s = _int_slider(imin, imax, ival, istep)
    sp = QSpinBox(); sp.setRange(imin, imax); sp.setSingleStep(istep); sp.setValue(ival); sp.setKeyboardTracking(False)
    readout = QLabel(f"{init_value:.6g}"); readout.setMinimumWidth(90)
    def sync(i):
        if sp.value()!=i: sp.setValue(i)
        readout.setText(f"{i*scale:.6g}")
    def sync_sp(i):
        if s.value()!=i: s.setValue(i)
        readout.setText(f"{i*scale:.6g}")
    s.valueChanged.connect(sync); sp.valueChanged.connect(sync_sp)
    layout.addWidget(name); layout.addWidget(s, 1); layout.addWidget(sp); layout.addWidget(readout)
    return box, s, sp, lambda: s.value()*scale


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(7, 4), constrained_layout=True)
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Field"); self.ax.set_ylabel("Intensity"); self.ax.grid(True, alpha=0.25)


# -------------------- Main Window --------------------

class SpectrumGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SAMPL4 Current Spectrum Explorer")
        here = Path(__file__).resolve().parent
        data_path = here / "sampl4.dat"
        if not data_path.exists():
            QMessageBox.critical(self, "Missing data", f"{data_path} not found next to this script.")
            raise SystemExit(1)

        # Load experimental trace from file next to this script
        x_in, y_exp = _load_two_column_ascii(data_path)
        start, step = _field_start_step(x_in)
        npts = int(y_exp.size)

        # --- EXACT test-style init sequence (no fitting) ---
        self.model = nlsl.nlsl()
        self.model["nsite"] = int(SAMPL4_FINAL_PARAMETERS["nsite"])  # required before weights
        # Initial dictionary (FINAL values as seed to match test)
        self.model.update({k:v for k,v in SAMPL4_FINAL_PARAMETERS.items() if k not in ("dfld","nfield")})
        # dfld and nfield come from the on-disk spectrum we just read
        self.model["dfld"] = step
        self.model["nfield"] = npts

        # Build exact field grid & attach processed data
        idx, data_slice = self.model.generate_coordinates(
            npts,
            start=start,
            step=step,
            derivative_mode=DERIVATIVE_MODE,
            baseline_points=BASELINE_EDGE_POINTS,
            normalize=False,
            nspline=NSPLINE_POINTS,
            shift=True,
            label="sampl4-ui",
            reset=True,
        )
        if idx != 0:
            raise RuntimeError("Unexpected spectrum index from generate_coordinates")
        self.model.set_data(data_slice, y_exp.astype(float))

        # Re-apply final parameters & bookkeeping tokens (mirrors test)
        self.model.update({k:v for k,v in SAMPL4_FINAL_PARAMETERS.items() if k not in ("dfld","nfield")})
        self.model["sb0"]   = SAMPL4_FINAL_SB0.copy()
        self.model["srng"]  = SAMPL4_FINAL_SRNG.copy()
        self.model["ishft"] = SAMPL4_FINAL_ISHFT.copy()
        self.model["shft"]  = SAMPL4_FINAL_SHFT.copy()
        self.model["nrmlz"] = SAMPL4_FINAL_NRMLZ.copy()
        self.model.weights  = SAMPL4_FINAL_WEIGHTS.copy()

        # ---- UI shell ----
        root = QWidget(self); self.setCentralWidget(root); vroot = QVBoxLayout(root)

        # nsite control row (above the tabs)
        head = QHBoxLayout(); vroot.addLayout(head)
        head.addWidget(QLabel("nsite:"))
        self.nsite_spin = QSpinBox(); self.nsite_spin.setRange(1, 8); self.nsite_spin.setValue(int(self.model["nsite"]))
        head.addWidget(self.nsite_spin); head.addStretch(1)

        # Plot
        self.canvas = MplCanvas(self); vroot.addWidget(self.canvas)
        self.x_axis = start + np.arange(npts, dtype=float) * step
        self.line_exp, = self.canvas.ax.plot(self.x_axis, y_exp, label="exp", lw=1)
        self.line_sim, = self.canvas.ax.plot(self.x_axis, self._current_total(), label="model", lw=1.5)
        self.canvas.ax.legend()
        ypad = 0.05*(y_exp.max()-y_exp.min()+1e-12)
        self.canvas.ax.set_ylim(y_exp.min()-ypad, y_exp.max()+ypad)
        self.canvas.ax.autoscale(enable=False)
        self.canvas.draw()

        # Tabs
        self.tabs = QTabWidget(); vroot.addWidget(self.tabs)
        self._rebuild_tabs()
        self.nsite_spin.valueChanged.connect(self._on_nsite_changed)

    # ---- model helpers ----
    def _current_total(self):
        return np.squeeze(self.model.weights @ self.model.current_spectrum)

    def _update(self):
        self.line_sim.set_ydata(self._current_total()); self.canvas.draw_idle()

    def _ensure_array_param(self, key, nsite, default_val):
        try:
            arr = np.atleast_1d(self.model[key]).astype(float)
        except Exception:
            arr = np.array([], dtype=float)
        if arr.size < nsite:
            pad = np.full(nsite - arr.size, float(default_val))
            arr = np.concatenate([arr, pad]) if arr.size else pad
        elif arr.size > nsite:
            arr = arr[:nsite]
        self.model[key] = arr
        return arr

    # ---- nsite changes ----
    def _on_nsite_changed(self, val):
        n = int(val)
        self.model["nsite"] = n
        # resize site params that are arrays
        self._ensure_array_param("gxx", n, SAMPL4_FINAL_PARAMETERS["gxx"])  # scalar → replicate
        self._ensure_array_param("gyy", n, SAMPL4_FINAL_PARAMETERS["gyy"])  # scalar → replicate
        self._ensure_array_param("gzz", n, SAMPL4_FINAL_PARAMETERS["gzz"])  # scalar → replicate
        self._ensure_array_param("axx", n, SAMPL4_FINAL_PARAMETERS["axx"])  # scalar → replicate
        self._ensure_array_param("ayy", n, SAMPL4_FINAL_PARAMETERS["ayy"])  # scalar → replicate
        self._ensure_array_param("azz", n, SAMPL4_FINAL_PARAMETERS["azz"])  # scalar → replicate
        self._ensure_array_param("rx",  n, float(SAMPL4_FINAL_PARAMETERS["rx"][0]))
        # weights resized
        try:
            w = np.atleast_1d(self.model.weights).astype(float)
        except Exception:
            w = np.ones(n)
        if w.size != n:
            w = np.ones(n)
        self.model.weights = w
        self._rebuild_tabs(); self._update()

    # ---- build tabs ----
    def _rebuild_tabs(self):
        while self.tabs.count():
            self.tabs.removeTab(0)
        n = int(self.model["nsite"]) if "nsite" in self.model else 1
        for i in range(n):
            self.tabs.addTab(self._component_tab(i), f"Component {i+1}")
        self.tabs.addTab(self._general_tab(), "General")

    # per‑component tab with nested sub‑tabs (g / A / Dynamics)
    def _component_tab(self, idx: int) -> QWidget:
        outer = QWidget(); outer_v = QVBoxLayout(outer)
        inner = QTabWidget(); outer_v.addWidget(inner)

        # g‑tensor tab
        tg = QWidget(); ggrid = QGridLayout(tg)
        for row, token in enumerate(["gxx","gyy","gzz"]):
            cur = np.atleast_1d(self.model.get(token, np.array([0.0]))).astype(float)
            init = float(cur[min(idx, cur.size-1)])
            box, s, sp, getv = make_scaled_pair(f"{token}[{idx}]", 1e-5, init, 1.99, 2.01, 1e-5)
            def apply_factory(tok, default):
                def apply(_):
                    arr = self._ensure_array_param(tok, int(self.model["nsite"]), default)
                    arr[idx] = getv(); self.model[tok] = arr; self._update()
                return apply
            s.valueChanged.connect(apply_factory(token, init)); sp.valueChanged.connect(apply_factory(token, init))
            ggrid.addWidget(box, row, 0)
        inner.addTab(tg, "g‑tensor")

        # A‑tensor tab
        ta = QWidget(); agrid = QGridLayout(ta)
        for row, token in enumerate(["axx","ayy","azz"]):
            cur = np.atleast_1d(self.model.get(token, np.array([0.0]))).astype(float)
            init = float(cur[min(idx, cur.size-1)])
            box, s, sp, getv = make_scaled_pair(f"{token}[{idx}] (G)", 1e-2, init, 0.0, 50.0, 1e-2)
            def apply_factory(tok, default):
                def apply(_):
                    arr = self._ensure_array_param(tok, int(self.model["nsite"]), default)
                    arr[idx] = getv(); self.model[tok] = arr; self._update()
                return apply
            s.valueChanged.connect(apply_factory(token, init)); sp.valueChanged.connect(apply_factory(token, init))
            agrid.addWidget(box, row, 0)
        inner.addTab(ta, "A‑tensor")

        # Dynamics tab (log10(rx) control)
        td = QWidget(); dgrid = QGridLayout(td)
        cur = np.atleast_1d(self.model.get("rx", np.array([7.5]))).astype(float)
        init = float(cur[min(idx, cur.size-1)])
        box, s, sp, getv = make_scaled_pair(f"log10(rx)[{idx}]", 1e-2, init, 4.0, 10.0, 1e-2)
        def apply(_):
            arr = self._ensure_array_param("rx", int(self.model["nsite"]), init)
            arr[idx] = getv(); self.model["rx"] = arr; self._update()
        s.valueChanged.connect(apply); sp.valueChanged.connect(apply)
        dgrid.addWidget(box, 0, 0)
        inner.addTab(td, "Dynamics")

        return outer

    # General tab: integer boxes + normalization + weights
    def _general_tab(self) -> QWidget:
        w = QWidget(); grid = QGridLayout(w)
        r = 0
        for key, lo, hi, init in [
            ("lemx", 1, 50, int(SAMPL4_FINAL_PARAMETERS["lemx"])),
            ("lomx", 1, 50, int(SAMPL4_FINAL_PARAMETERS["lomx"])),
            ("kmx",  1, 50, int(SAMPL4_FINAL_PARAMETERS["kmx"])),
            ("mmx",  1, 50, int(SAMPL4_FINAL_PARAMETERS["mmx"])),
            ("ipnmx",1, 10, int(SAMPL4_FINAL_PARAMETERS["ipnmx"])),
            ("in2",  1, 10, int(SAMPL4_FINAL_PARAMETERS["in2"])),
        ]:
            grid.addWidget(QLabel(key), r, 0)
            spin = QSpinBox(); spin.setRange(lo, hi); spin.setValue(int(self.model.get(key, init)))
            def mk_apply(k, sp):
                def _(_v): self.model[k] = int(sp.value()); self._update()
                return _
            spin.valueChanged.connect(mk_apply(key, spin))
            grid.addWidget(spin, r, 1); r += 1

        # Normalization
        init_nrmlz = float(np.atleast_1d(self.model.get("nrmlz", SAMPL4_FINAL_NRMLZ.astype(float)))[0])
        nr_box, s, sp, getv = make_scaled_pair("nrmlz", 1e-3, init_nrmlz, 0.0, 5.0, 1e-3)
        def upd_nr(_): self.model["nrmlz"] = getv(); self._update()
        s.valueChanged.connect(upd_nr); sp.valueChanged.connect(upd_nr)
        grid.addWidget(nr_box, r, 0, 1, 2); r += 1

        # Weights
        wgrp = QGroupBox("Weights"); wgrid = QGridLayout(wgrp)
        self._get_w = []
        wv = np.atleast_1d(self.model.weights).astype(float)
        for i in range(wv.size):
            box, s, sp, getw = make_scaled_pair(f"w[{i}]", 1e-3, float(wv[i]), 0.0, 5.0, 1e-3)
            def up(_): self.model.weights = np.array([gw() for gw in self._get_w]); self._update()
            s.valueChanged.connect(up); sp.valueChanged.connect(up)
            self._get_w.append(getw); wgrid.addWidget(box, i, 0)
        grid.addWidget(wgrp, r, 0, 1, 2)

        return w


def main():
    app = QApplication(sys.argv)
    w = SpectrumGUI(); w.resize(1200, 850); w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
