from pathlib import Path
import numpy as np

# Use Qt backend BEFORE importing pyplot
import matplotlib
matplotlib.use("Qt5Agg")  # PyQt5
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from PyQt5 import QtCore, QtWidgets

import nlsl
from nlsl.data import process_spectrum

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

# Site populations and spectral metadata for reproducing the converged simulation
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

        # ---- Build UI ----
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vbox = QtWidgets.QVBoxLayout(central)

        # Matplotlib figure/canvas (transparent background)
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.fig.patch.set_alpha(0.0)  # figure transparent
        self.ax.set_facecolor("none")  # axes transparent
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setAutoFillBackground(False)
        self.canvas.setStyleSheet("background: transparent")
        vbox.addWidget(self.canvas, stretch=1)

        # Tab widget below plot
        tabs = QtWidgets.QTabWidget()
        vbox.addWidget(tabs, stretch=0)

        # --- General tab: weight sliders ---
        general = QtWidgets.QWidget()
        tabs.addTab(general, "general")
        form = QtWidgets.QFormLayout(general)

        self.s1 = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.s1.setRange(0, 100)
        self.s1.setValue(int(SAMPL4_FINAL_WEIGHTS[0] * 100))
        self.l1 = QtWidgets.QLabel()

        self.s2 = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.s2.setRange(0, 100)
        self.s2.setValue(int(SAMPL4_FINAL_WEIGHTS[1] * 100))
        self.l2 = QtWidgets.QLabel()

        # rows
        row1 = QtWidgets.QHBoxLayout()
        row1.addWidget(self.s1)
        row1.addWidget(self.l1)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.s2)
        row2.addWidget(self.l2)

        form.addRow("weight 1", row1)
        form.addRow("weight 2", row2)

        # ---- Initial plot (keep Line2D handles to update ydata only) ----
        colours = ["#1f77b4", "#2ca02c"]  # stylistic example
        self.exp_line, = self.ax.plot(
            self.x, self.y_exp, color="black", linewidth=1.0, label="experimental"
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
                label=f"component {i+1}",
            )
            self.comp_lines.append(line)

        (self.total_line,) = self.ax.plot(
            self.x, total, color="#d62728", linewidth=2.0, alpha=0.8, label="sum"
        )

        self.ax.set_xlabel("Magnetic field (G)")
        self.ax.set_ylabel("Intensity (arb. units)")
        self.ax.set_title("sampl4 components from model.current_spectrum (no fit)")
        self.ax.legend(loc="upper right")
        self.ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        self.fig.tight_layout()
        self.canvas.draw()

        # connect signals AFTER first draw
        self.s1.valueChanged.connect(self.on_weights_changed)
        self.s2.valueChanged.connect(self.on_weights_changed)

        # set initial labels
        self.update_weight_labels()

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
            label="sampl4-single-eval",
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

    # ---- Weight change handler: update only line ydata ----
    @QtCore.pyqtSlot()
    def on_weights_changed(self):
        # raw sliders 0..100 -> normalize to sum=1
        w1 = float(self.s1.value())
        w2 = float(self.s2.value())
        s = max(w1 + w2, 1e-12)
        self.weights = np.array([w1 / s, w2 / s], dtype=float)
        self.update_weight_labels()

        # recompute components/total from cached site spectra
        comp = self.weights[:, None] * self.site_spectra
        total = self.weights @ self.site_spectra

        # update line ydata only (no re-plotting)
        for i, line in enumerate(self.comp_lines):
            line.set_ydata(comp[i])
        self.total_line.set_ydata(total)
        self.canvas.draw_idle()

    def update_weight_labels(self):
        self.l1.setText(f"{self.weights[0]:.3f}")
        self.l2.setText(f"{self.weights[1]:.3f}")


def main():
    app = QtWidgets.QApplication([])
    # Let the Qt window show its own background; make central widget default palette
    app.setStyleSheet("QMainWindow{background:#d6d6d6}")  # subtle Qt grey (optional)
    w = MainWindow()
    w.resize(1000, 720)
    w.show()
    app.exec_()


if __name__ == "__main__":
    main()
