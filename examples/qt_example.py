#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Qt GUI to explore a spectrum using the nlsl model.

• Loads two‑column ASCII experimental data from the same folder (e.g., sampl4.dat).
• Uses model.load_data(...) so the field axis is owned by the model.
• Plots experimental vs model; only the model Line2D ydata is updated on changes.
• nsite control sits above the tabs and determines how many top‑level component tabs appear.
• All parameters from INITIAL_PARAMETERS are exposed, except nsite.
  - Per‑component tabs (for each site): g‑tensor, A‑tensor, Dynamics (rx, rbar, gib0).
  - General tab: integer boxes for lemx/lomx/kmx/mmx/ipnmx and in2, plus
    normalization and per‑site weight sliders.
"""

import sys
from pathlib import Path
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QSlider, QTabWidget, QGroupBox, QGridLayout, QSpinBox,
    QMessageBox
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import nlsl

# ---------------- Data loader ----------------

def _load_experimental(dirpath: Path):
    """Return x,y from whitespace-delimited ASCII; prefer sampl4.dat if present."""
    last_err = None
    for fname in ["sampl4.dat", "experimental.dat", "experimental.txt", "experimental.csv"]:
        p = dirpath / fname
        if not p.exists():
            continue
        try:
            arr = np.loadtxt(p)
            arr = np.atleast_2d(arr)
            if arr.shape[1] == 1:
                y = arr[:, 0]
                x = np.arange(y.size, dtype=float)
            else:
                x, y = arr[:, 0], arr[:, 1]
            return np.asarray(x, float), np.asarray(y, float)
        except Exception as e:
            last_err = e
    raise FileNotFoundError(f"No valid experimental data found next to script: {last_err}")

# ---------------- Tiny UI helpers ----------------

def make_int_slider(min_v, max_v, value, step=1):
    s = QSlider(QtCore.Qt.Horizontal)
    s.setMinimum(int(min_v))
    s.setMaximum(int(max_v))
    s.setSingleStep(int(step))
    s.setPageStep(int(step))
    s.setValue(int(value))
    return s


def make_scaled_pair(label, scale, init_value, min_v, max_v, step):
    """Float slider via integer slider + scale. Returns (row_widget, slider, spinbox, getter)."""
    box = QWidget(); layout = QHBoxLayout(box)
    name = QLabel(label); name.setMinimumWidth(120)

    imin, imax = int(min_v/scale), int(max_v/scale)
    ival, istep = int(init_value/scale), max(1, int(step/scale))
    s = make_int_slider(imin, imax, ival, istep)
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
        self.fig = Figure(figsize=(7,4), constrained_layout=True)
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Field (arb)"); self.ax.set_ylabel("Intensity (arb)"); self.ax.grid(True, alpha=0.25)


# ---------------- Main window ----------------

class SpectrumGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Current Spectrum Explorer")
        here = Path(__file__).resolve().parent
        try:
            x_in, y_exp = _load_experimental(here)
        except Exception as e:
            QMessageBox.critical(self, "Data load error", str(e)); raise

        # ---- Model (reasonable initial parameters) ----
        self.model = nlsl.nlsl()
        self.INIT = {
            "nsite": 2,
            "in2": 2,
            # g and A tensors (TEMPO-like)
            "gxx": 2.0089, "gyy": 2.0063, "gzz": 2.0021,
            "axx": 5.0,   "ayy": 5.0,   "azz": 33.0,
            # grid / integration controls
            "lemx": 12, "lomx": 10, "kmx": 7, "mmx": 7, "ipnmx": 2,
            # dynamics-ish site params
            "gib0": 0.5,
            "rx": np.array([np.log10(3.0e8), np.log10(1.0e7)]),
            "rbar": np.array([2.0, 2.0]),
        }
        self.model.update(self.INIT)

        # prefer on-disk file; otherwise temp from x_in,y_exp
        data_path = here / "sampl4.dat" if (here/"sampl4.dat").exists() else None
        if data_path is None:
            tmp = here / "_ui_temp.dat"
            np.savetxt(tmp, np.column_stack([x_in, y_exp])); data_path = tmp

        self.model.load_data(
            data_path,
            nspline=200,
            bc_points=20,
            shift=True,
            normalize=False,
            derivative_mode=1,
        )
        try:
            self.model.fit(); self.model.fit()
        except Exception:
            pass
        # equal weights (size nsite)
        self.model.weights = np.ones(int(self.model["nsite"]))

        fields = np.asarray(self.model.field_axes[0], float)
        experimental = np.squeeze(self.model.experimental_data)

        # ---- UI shell ----
        root = QWidget(self); self.setCentralWidget(root); vroot = QVBoxLayout(root)

        # nsite control row
        head = QHBoxLayout(); vroot.addLayout(head)
        head.addWidget(QLabel("nsite:"))
        self.nsite_spin = QSpinBox(); self.nsite_spin.setRange(1, 8); self.nsite_spin.setValue(int(self.model["nsite"]))
        head.addWidget(self.nsite_spin); head.addStretch(1)

        # plot
        self.canvas = MplCanvas(self); vroot.addWidget(self.canvas)
        self.x_axis = fields
        self.line_exp, = self.canvas.ax.plot(self.x_axis, experimental, label="exp", lw=1)
        self.line_sim, = self.canvas.ax.plot(self.x_axis, self._current_total(), label="model", lw=1.5)
        self.canvas.ax.legend()
        ypad = 0.05*(experimental.max()-experimental.min()+1e-12)
        self.canvas.ax.set_ylim(experimental.min()-ypad, experimental.max()+ypad)
        self.canvas.ax.autoscale(enable=False)
        self.canvas.draw()

        # tabs
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
        # resize site params
        self._ensure_array_param("gxx", n, 2.0089)
        self._ensure_array_param("gyy", n, 2.0063)
        self._ensure_array_param("gzz", n, 2.0021)
        self._ensure_array_param("axx", n, 5.0)
        self._ensure_array_param("ayy", n, 5.0)
        self._ensure_array_param("azz", n, 33.0)
        self._ensure_array_param("rx",  n, np.log10(1e8))
        self._ensure_array_param("rbar",n, 2.0)
        self._ensure_array_param("gib0",n, 0.5)
        # weights
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

    # per‑component tab with nested sub‑tabs
    def _component_tab(self, idx: int) -> QWidget:
        outer = QWidget(); outer_v = QVBoxLayout(outer)
        inner = QTabWidget(); outer_v.addWidget(inner)

        # g‑tensor tab
        tg = QWidget(); ggrid = QGridLayout(tg)
        def add_g(token, row):
            cur = np.atleast_1d(self.model.get(token, np.array([0.0]))).astype(float)
            init = float(cur[min(idx, cur.size-1)])
            box, s, sp, getv = make_scaled_pair(f"{token}[{idx}]", 1e-5, init, 1.99, 2.01, 1e-5)
            def apply(_):
                arr = self._ensure_array_param(token, int(self.model["nsite"]), init)
                arr[idx] = getv(); self.model[token] = arr; self._update()
            s.valueChanged.connect(apply); sp.valueChanged.connect(apply)
            ggrid.addWidget(box, row, 0)
        add_g("gxx", 0); add_g("gyy", 1); add_g("gzz", 2)
        inner.addTab(tg, "g‑tensor")

        # A‑tensor tab
        ta = QWidget(); agrid = QGridLayout(ta)
        def add_a(token, row):
            cur = np.atleast_1d(self.model.get(token, np.array([0.0]))).astype(float)
            init = float(cur[min(idx, cur.size-1)])
            box, s, sp, getv = make_scaled_pair(f"{token}[{idx}] (G)", 1e-2, init, 0.0, 50.0, 1e-2)
            def apply(_):
                arr = self._ensure_array_param(token, int(self.model["nsite"]), init)
                arr[idx] = getv(); self.model[token] = arr; self._update()
            s.valueChanged.connect(apply); sp.valueChanged.connect(apply)
            agrid.addWidget(box, row, 0)
        add_a("axx", 0); add_a("ayy", 1); add_a("azz", 2)
        inner.addTab(ta, "A‑tensor")

        # Dynamics tab
        td = QWidget(); dgrid = QGridLayout(td)
        def add_dyn(token, label, scale, lo, hi, step, row, default_val):
            cur = np.atleast_1d(self.model.get(token, np.array([default_val]))).astype(float)
            init = float(cur[min(idx, cur.size-1)])
            box, s, sp, getv = make_scaled_pair(f"{label}[{idx}]", scale, init, lo, hi, step)
            def apply(_):
                arr = self._ensure_array_param(token, int(self.model["nsite"]), default_val)
                arr[idx] = getv(); self.model[token] = arr; self._update()
            s.valueChanged.connect(apply); sp.valueChanged.connect(apply)
            dgrid.addWidget(box, row, 0)
        add_dyn("rx",   "log10(rx)", 1e-2, 4.0, 10.0, 1e-2, 0, np.log10(1e8))
        add_dyn("rbar", "rbar (nm)",  1e-3, 0.1, 5.0,  1e-3, 1, 2.0)
        add_dyn("gib0", "gib0",       1e-3, 0.0, 5.0,  1e-3, 2, 0.5)
        inner.addTab(td, "Dynamics")

        return outer

    # General tab: integer boxes + normalization + weights
    def _general_tab(self) -> QWidget:
        w = QWidget(); grid = QGridLayout(w)
        r = 0
        # Integer boxes (leave these as int editors): lemx/lomx/kmx/mmx/ipnmx and in2
        for key, lo, hi, init in [
            ("lemx", 1, 50, int(self.INIT.get("lemx", 12))),
            ("lomx", 1, 50, int(self.INIT.get("lomx", 10))),
            ("kmx",  1, 50, int(self.INIT.get("kmx", 7))),
            ("mmx",  1, 50, int(self.INIT.get("mmx", 7))),
            ("ipnmx",1, 10, int(self.INIT.get("ipnmx", 2))),
            ("in2",  1, 10, int(self.INIT.get("in2", 2))),
        ]:
            grid.addWidget(QLabel(key), r, 0)
            spin = QSpinBox(); spin.setRange(lo, hi); spin.setValue(int(self.model.get(key, init)))
            def mk_apply(k, sp):
                def _(_v): self.model[k] = int(sp.value()); self._update()
                return _
            spin.valueChanged.connect(mk_apply(key, spin))
            grid.addWidget(spin, r, 1); r += 1

        # Normalization slider
        init_nrmlz = float(self.model.get("nrmlz", 1.0))
        nr_box, s, sp, getv = make_scaled_pair("nrmlz", 1e-3, init_nrmlz, 0.0, 5.0, 1e-3)
        def upd_nr(_): self.model["nrmlz"] = getv(); self._update()
        s.valueChanged.connect(upd_nr); sp.valueChanged.connect(upd_nr)
        grid.addWidget(nr_box, r, 0, 1, 2); r += 1

        # Weights group
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


if __name__=='__main__':
    main()
