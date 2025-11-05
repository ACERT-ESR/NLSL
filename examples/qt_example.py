from pathlib import Path
import re
import numpy as np
from PyQt5 import QtCore, QtWidgets
import nlsl
from nlsl.data import process_spectrum

# Use Qt backend BEFORE importing pyplot
import matplotlib

matplotlib.use("Qt5Agg")  # PyQt5
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
)


# --- Hard-coded SAMPL4 setup (no imports from tests or references) ---
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

# Final parameters (mirror of classic runfile-4 solution)
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
    "dfld": None,  # set from processed data below
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
    "nfield": NSPLINE_POINTS,
    "ideriv": DERIVATIVE_MODE,
    "iwflg": 0,
    "igflg": 0,
    "iaflg": 0,
    "jkmn": 0,
    "jmmn": 0,
    "irflg": 2,
    "ndim": 156,
}

# Site populations and spectral metadata for reproducing the converged
# simulation
SAMPL4_FINAL_WEIGHTS = np.array([0.2848810, 0.7155313])
SAMPL4_FINAL_ISHFT = np.array([1], dtype=np.int32)
SAMPL4_FINAL_SHFT = np.array([0.0])
SAMPL4_FINAL_NRMLZ = np.array([0], dtype=np.int32)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SAMPL4 components (Qt + Matplotlib)")

        # ---- Prepare model/data once ----
        (
            self.x,
            self.y_exp,
            self.model,
            self.site_spectra,
        ) = self.prepare_model()
        self.nsite = (
            int(self.model["nsite"])
            if "nsite" in self.model
            else len(self.site_spectra)
        )

        # ---- Build UI ----
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vbox = QtWidgets.QVBoxLayout(central)

        # Use a splitter so the plot keeps most of the space
        splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        vbox.addWidget(splitter)

        # Matplotlib figure/canvas (transparent figure, white axes)
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.fig.patch.set_alpha(0.0)  # figure transparent
        self.ax.set_facecolor("white")  # axes transparent
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setAutoFillBackground(False)
        self.canvas.setStyleSheet("background: transparent")
        splitter.addWidget(self.canvas)

        # Tab widget below plot
        self.tabs = QtWidgets.QTabWidget()
        # keep tabs compact
        self.tabs.setSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum
        )
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)

        # --- General tab: weight sliders (0..1, complementary) ---
        general = QtWidgets.QWidget()
        self.tabs.addTab(general, "general")
        form = QtWidgets.QFormLayout(general)

        self.weight_steps = 1000
        self.s1 = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.s1.setRange(0, self.weight_steps)
        self.s1.setValue(int(SAMPL4_FINAL_WEIGHTS[0] * self.weight_steps))
        self.l1 = QtWidgets.QLabel()

        self.s2 = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.s2.setRange(0, self.weight_steps)
        self.s2.setValue(int(SAMPL4_FINAL_WEIGHTS[1] * self.weight_steps))
        self.l2 = QtWidgets.QLabel()

        row1 = QtWidgets.QHBoxLayout()
        row1.addWidget(self.s1)
        row1.addWidget(self.l1)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.s2)
        row2.addWidget(self.l2)

        form.addRow("weight 1", row1)
        form.addRow("weight 2", row2)

        # --- Integer params (now inside GENERAL tab, 4 columns + scroll) ---
        int_group = QtWidgets.QGroupBox("integer params")
        int_group_layout = QtWidgets.QVBoxLayout(int_group)
        int_scroll = QtWidgets.QScrollArea()
        int_scroll.setWidgetResizable(True)
        int_group_layout.addWidget(int_scroll)
        int_inner = QtWidgets.QWidget()
        int_scroll.setWidget(int_inner)
        int_grid = QtWidgets.QGridLayout(int_inner)
        self.int_boxes = {}
        int_keys = [
            k
            for k, v in SAMPL4_FINAL_PARAMETERS.items()
            if isinstance(v, (int, np.integer))
        ]
        int_keys.sort()
        for idx, k in enumerate(int_keys):
            v = SAMPL4_FINAL_PARAMETERS[k]
            r = idx // 4
            c = (idx % 4) * 2
            lbl = QtWidgets.QLabel(k)
            spin = QtWidgets.QSpinBox()
            spin.setRange(-999999, 999999)
            spin.setValue(int(v))
            spin.valueChanged.connect(self._make_int_param_handler(k))
            int_grid.addWidget(lbl, r, c)
            int_grid.addWidget(spin, r, c + 1)
            self.int_boxes[k] = spin
        form.addRow(int_group)

        # ---- Site tabs with nested subtabs for tensor groups ----

        # Detect tensor groups by keys ending with xx/yy/zz
        self.tensor_groups = self._detect_tensor_groups(
            SAMPL4_FINAL_PARAMETERS
        )
        self.tensor_ranges = {
            "g": (1.95, 2.05, 1000, 4),  # (min,max,steps,decimals)
            "a": (0.0, 50.0, 1000, 3),
            "w": (0.0, 20.0, 1000, 3),
            "r": (0.0, 8.0, 1000, 3),
            "p": (0.0, 8.0, 1000, 3),
            "d": (0.0, 8.0, 1000, 3),
        }

        # list of dicts: per site -> {fam: (subwidget, index)}
        self.site_family_tabs = []
        for site in range(self.nsite):
            site_tab = QtWidgets.QWidget()
            self.tabs.addTab(site_tab, f"site {site + 1}")
            site_layout = QtWidgets.QVBoxLayout(site_tab)
            sub = QtWidgets.QTabWidget()
            site_layout.addWidget(sub)

            fam_map = {}
            for base, comps in self.tensor_groups.items():
                label = {
                    "g": "g tensor",
                    "a": "A tensor",
                    "w": "linewidth tensor",
                    "r": "R tensor (rotational diffusion)",
                    "p": "P tensor (non-Brownian)",
                    "d": "Dj tensor (anisotropic viscosity)",
                }.get(base, f"{base} tensor")
                page = QtWidgets.QWidget()
                grid = QtWidgets.QGridLayout(page)
                idx = sub.addTab(page, label)
                fam_map[base] = (sub, idx)

                vmin, vmax, steps, dec = self.tensor_ranges.get(
                    base, (0.0, 1.0, 1000, 3)
                )
                # choose which component suffix set exists for this family
                comp_set = [
                    s
                    for s in ("xx", "yy", "zz", "x", "y", "z")
                    if (base + s) in comps
                ]
                for r, comp in enumerate(comp_set):
                    key = base + comp
                    if key not in comps:
                        continue
                    slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
                    slider.setRange(0, steps)
                    current = self._site_value_from_model(key, site)
                    slider.setValue(
                        self._float_to_slider(current, vmin, vmax, steps)
                    )
                    label_val = QtWidgets.QLabel(f"{current:.{dec}f}")
                    slider.valueChanged.connect(
                        self._make_tensor_handler(
                            key, site, vmin, vmax, steps, label_val
                        )
                    )
                    grid.addWidget(QtWidgets.QLabel(key), r, 0)
                    grid.addWidget(slider, r, 1)
                    grid.addWidget(label_val, r, 2)
            self.site_family_tabs.append(fam_map)

        # After building site tabs, enforce initial visibility based on ipdf
        self._refresh_family_visibility()

        # ---- Initial plot (keep Line2D handles to update ydata only) ----
        colours = [
            "#1f77b4",
            "#2ca02c",
            "#9467bd",
            "#8c564b",
        ]  # stylistic example
        (self.exp_line,) = self.ax.plot(
            self.x,
            self.y_exp,
            color="black",
            linewidth=1.0,
            label="experimental",
        )

        # weighted components and total using current weights
        self.weights = SAMPL4_FINAL_WEIGHTS.astype(float).copy()
        comp = self.weights[:, None] * self.site_spectra
        total = self.weights @ self.site_spectra

        self.comp_lines = []
        for i in range(comp.shape[0]):
            (line,) = self.ax.plot(
                self.x,
                comp[i],
                color=colours[i % len(colours)],
                linewidth=1.2,
                alpha=0.7,
                label=f"component {i + 1}",
            )
            self.comp_lines.append(line)

        (self.total_line,) = self.ax.plot(
            self.x,
            total,
            color="#d62728",
            linewidth=2.0,
            alpha=0.8,
            label="sum",
        )

        self.ax.set_xlabel("Magnetic field (G)")
        self.ax.set_ylabel("Intensity (arb. units)")
        self.ax.set_title(
            "sampl4 components from model.current_spectrum (no fit)"
        )
        self.ax.legend(loc="upper right")
        self.ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        self.fig.tight_layout()
        self.canvas.draw()

        # connect signals AFTER first draw
        self.s1.valueChanged.connect(self.on_weight1_changed)
        self.s2.valueChanged.connect(self.on_weight2_changed)

        # set initial labels and keep sliders consistent to sum=1
        self.update_weight_labels()
        self._sync_partner_sliders(init=True)

    # ---- Data/model preparation (no fitting) ----
    def prepare_model(self):
        data_path = Path(__file__).resolve().parent / "sampl4.dat"
        proc = process_spectrum(
            data_path,
            NSPLINE_POINTS,
            BASELINE_EDGE_POINTS,
            derivative_mode=DERIVATIVE_MODE,
            normalize=False,
        )
        field_start = float(proc.start)
        field_step = float(proc.step)
        point_count = proc.y.size
        y_exp = proc.y.copy()

        params = dict(SAMPL4_FINAL_PARAMETERS)
        params["dfld"] = field_step
        sb0 = np.array([params["b0"]])
        srng = np.array([params["range"]])

        model = nlsl.nlsl()
        model["nsite"] = 2
        index, sl = model.generate_coordinates(
            point_count,
            start=field_start,
            step=field_step,
            derivative_mode=DERIVATIVE_MODE,
            baseline_points=BASELINE_EDGE_POINTS,
            normalize=False,
            nspline=NSPLINE_POINTS,
            shift=True,
            label=f"sampl4-eval-{point_count}",
            reset=True,
        )
        model.set_data(sl, y_exp[:point_count])
        model.update(params)
        model["sb0"] = sb0
        model["srng"] = srng
        model["ishft"] = SAMPL4_FINAL_ISHFT
        model["shft"] = SAMPL4_FINAL_SHFT
        model["nrmlz"] = SAMPL4_FINAL_NRMLZ
        model.weights = SAMPL4_FINAL_WEIGHTS

        site_spectra = model.current_spectrum  # (nsite, npts)
        x = field_start + field_step * np.arange(point_count)
        return x, y_exp, model, site_spectra

    # ---- Weight change handlers ----
    @QtCore.pyqtSlot()
    def on_weight1_changed(self, value=None):
        # s1 drove change; enforce w1 + w2 = 1 and update partner slider
        w1 = self.s1.value() / self.weight_steps
        w2 = max(0.0, 1.0 - w1)
        # block signals to avoid recursion
        self.s2.blockSignals(True)
        self.s2.setValue(int(round(w2 * self.weight_steps)))
        self.s2.blockSignals(False)
        self.weights = np.array([w1, w2], dtype=float)
        self.update_weight_labels()
        self._recompute_and_redraw()

    def on_weight2_changed(self, value=None):
        w2 = self.s2.value() / self.weight_steps
        w1 = max(0.0, 1.0 - w2)
        self.s1.blockSignals(True)
        self.s1.setValue(int(round(w1 * self.weight_steps)))
        self.s1.blockSignals(False)
        self.weights = np.array([w1, w2], dtype=float)
        self.update_weight_labels()
        self._recompute_and_redraw()

    # ---- Tensor slider handlers ----
    def _make_tensor_handler(self, key, site, vmin, vmax, steps, label_widget):
        def handler(value_int):
            val = self._slider_to_float(value_int, vmin, vmax, steps)
            # numeric label formatting by family
            if key.startswith("g"):
                label_widget.setText(f"{val:.4f}")
            else:
                label_widget.setText(f"{val:.3f}")
            self._set_site_param(key, site, val)
            # current_spectrum is a PROPERTY; re-read after param change
            self.site_spectra = self.model.current_spectrum
            self._recompute_and_redraw(update_components_only=True)

        return handler

    # ---- Helpers for model param array management ----
    def _site_value_from_model(self, key, site):
        n = self.nsite
        v = self.model[key]
        arr = (
            np.full(n, float(v))
            if np.ndim(v) == 0
            else np.array(v, dtype=float)
        )
        if arr.size != n:
            arr = np.resize(arr, n)
        return float(arr[site])

    def _set_site_param(self, key, site, value):
        n = self.nsite
        v = self.model[key]
        arr = (
            np.full(n, float(v))
            if np.ndim(v) == 0
            else np.array(v, dtype=float)
        )
        if arr.size != n:
            arr = np.resize(arr, n)
        arr[site] = float(value)
        self.model[key] = arr

    # ---- Integer param handler factory ----
    def _make_int_param_handler(self, key):
        def handler(val):
            self.model[key] = int(val)
            # Special logic: ipdf controls which diffusion families are active
            if key == "ipdf":
                self._refresh_family_visibility()
            # Re-evaluate spectrum after model integer change
            self.site_spectra = self.model.current_spectrum
            self._recompute_and_redraw(update_components_only=True)

        return handler

    # ---- Slider <-> float mapping ----
    @staticmethod
    def _slider_to_float(value_int, vmin, vmax, steps):
        return vmin + (vmax - vmin) * (value_int / steps)

    @staticmethod
    def _float_to_slider(value_float, vmin, vmax, steps):
        # clamp then map
        vf = max(min(value_float, vmax), vmin)
        return int(round((vf - vmin) / (vmax - vmin) * steps))

    # ---- Tensor-group detection ----
    @staticmethod
    def _detect_tensor_groups(param_dict):
        """
        Detect true tensor families only. Excludes integer/basis keys that
        begin with i/k/l/m (e.g., ipnmx, kmx, lomx, mxy, mzz, etc.). Only
        families starting with g, a, w, r, p, d are considered, and they may
        use xx/yy/zz or x/y/z suffixes.
        """
        groups = {}
        pat = re.compile(r"^(?P<fam>[gawrpd])[a-z]*(?:xx|yy|zz|x|y|z)$")
        for k in param_dict.keys():
            m = pat.match(k)
            if not m:
                continue
            fam = m.group("fam")  # one of g,a,w,r,p,d
            groups.setdefault(fam, set()).add(k)
        return groups

    # ---- Family visibility rules driven by ipdf ----
    def _refresh_family_visibility(self):
        ipdf = int(self.model.get("ipdf", 0))
        # Mapping per docs:
        # ipdf=0 → Brownian (R tensor active)
        # ipdf=1 → Non‑Brownian (P tensor active)
        # ipdf=2 → Anisotropic viscosity (Dj tensor active)
        enable_by_fam = {
            "g": True,
            "a": True,
            "w": True,
            "r": ipdf == 0,
            "p": ipdf == 1,
            "d": ipdf == 2,
        }
        for site, fam_map in enumerate(self.site_family_tabs):
            for fam, (sub, idx) in fam_map.items():
                if hasattr(sub, "setTabVisible"):
                    sub.setTabVisible(idx, bool(enable_by_fam.get(fam, True)))
                else:
                    sub.setTabEnabled(idx, bool(enable_by_fam.get(fam, True)))

    # ---- Recompute components & update line ydata only ----
    def _recompute_and_redraw(self, update_components_only=False):
        if not update_components_only:
            # If weights changed but site spectra didn't, reuse cached spectra
            self.site_spectra = self.site_spectra
        comp = self.weights[:, None] * self.site_spectra
        total = self.weights @ self.site_spectra
        # update component lines
        for i, line in enumerate(self.comp_lines):
            if i < comp.shape[0]:
                line.set_ydata(comp[i])
        # update total line
        self.total_line.set_ydata(total)
        self.canvas.draw_idle()

    def update_weight_labels(self):
        self.l1.setText(f"{self.weights[0]:.3f}")
        self.l2.setText(f"{self.weights[1]:.3f}")

    def _sync_partner_sliders(self, init=False):
        # Ensure sliders reflect exact complementarity
        w1 = self.s1.value() / self.weight_steps
        w2 = max(0.0, 1.0 - w1)
        self.s2.blockSignals(True)
        self.s2.setValue(int(round(w2 * self.weight_steps)))
        self.s2.blockSignals(False)
        if init:
            self.weights = np.array([w1, w2], dtype=float)
            self.update_weight_labels()
        self.l1.setText(f"{self.weights[0]:.3f}")
        self.l2.setText(f"{self.weights[1]:.3f}")


def main():
    app = QtWidgets.QApplication([])
    # Let the Qt window show its own background; make central widget default
    # palette
    app.setStyleSheet(
        "QMainWindow{background:#d6d6d6}"
    )  # subtle Qt grey (optional)
    w = MainWindow()
    w.resize(1100, 820)
    w.show()
    app.exec_()


if __name__ == "__main__":
    main()
