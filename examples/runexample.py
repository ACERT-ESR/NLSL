#!/usr/bin/python
from pylab import *
from subprocess import Popen, PIPE, STDOUT
import os
import nlsl
import sys
import numpy as np

rc("font", size=18)


def read_column_data(filename):
    fp = open(filename, "r")
    data = []
    for j in fp.readlines():
        data.append(j.split())
    data = array(data, dtype=double)
    fp.close()
    return data


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise OSError("This function takes one argument -- the number of the example")
    print("about to run nlsl")
    filename_base = "sampl" + sys.argv[1]
    examples_dir = os.path.dirname(__file__)
    os.chdir(examples_dir)
    data_files_out = []
    n = nlsl.nlsl()

    def run_file(thisfp):
        for thisline in thisfp.readlines():
            if thisline[:5] == "call ":
                fp_called = open(thisline[5:].strip())
                run_file(fp_called)
                fp_called.close()
            elif thisline[:5] == "data ":
                n.procline(thisline)
                data_files_out.append(thisline[5:].strip().split(" ")[0])
            else:
                n.procline(thisline)
        thisfp.close()

    run_file(open(filename_base + ".run"))
    print("final parameters:")
    print("→" * 20)
    print("\n".join(f"{k:20s}{v}" for k, v in n.items()))
    print("final fit parameters:")
    print("→" * 20)
    print("\n".join(f"{k:20s}{v}" for k, v in n.fit_params.items()))

    # Compute relative RMS error for each output spectrum
    rms_values = []
    rms_sq_total = 0.0
    exp_sq_total = 0.0
    for thisdatafile in data_files_out:
        data_calc = read_column_data(thisdatafile + ".spc")
        exp_sq = np.sum(data_calc[:, 1] ** 2)
        rms_sq = np.sum((data_calc[:, 2] - data_calc[:, 1]) ** 2)
        exp_sq_total += exp_sq
        rms_sq_total += rms_sq
        if exp_sq > 0:
            rms = np.sqrt(rms_sq) / np.sqrt(exp_sq)
            rms_values.append((thisdatafile, rms))
            print(f"{thisdatafile}: relative rms = {rms:.5g}")
    if exp_sq_total > 0 and len(rms_values) > 1:
        rel_rms = np.sqrt(rms_sq_total) / np.sqrt(exp_sq_total)
        print(f"overall relative rms = {rel_rms:.5g}")

    # print "result:",nlsl.parameters.asdict
    fig = figure(figsize=(10, 5 * len(data_files_out)))
    fig.subplots_adjust(left=0.15)
    fig.subplots_adjust(right=0.6)
    fig.subplots_adjust(hspace=0.3)
    if len(data_files_out) == 1:
        fig.subplots_adjust(bottom=0.15)
    fig.subplots_adjust(top=0.95)
    for j, thisdatafile in enumerate(data_files_out):
        subplot(len(data_files_out), 1, j + 1)

        # fig.add_axes([0.1,0.1,0.6,0.8]) # l b w h
        def show_spc(filename, linestyles=["k", "r"], show_components=True):
            data = read_column_data(filename)
            fields = data[:, 0]
            experimental = data[:, 1]
            fit = data[:, 2]
            integral_of_spectrum = cumsum(experimental)
            normalization = abs(sum(integral_of_spectrum))
            if data.shape[1] > 3:
                components = data[:, 3:]
            else:
                components = None
            plot(
                fields,
                experimental / normalization,
                linestyles[0],
                linewidth=1,
                label="experimental",
            )
            plot(
                fields,
                fit / normalization,
                linestyles[1],
                alpha=0.5,
                linewidth=2,
                label="fit (" + filename + ")",
            )
            max_of_fit = max(fit) / normalization
            if components is not None and show_components:
                plot(
                    fields,
                    components / normalization,
                    alpha=0.3,
                    linewidth=1,
                    label="component",
                )
            ax = gca()
            ylims = ax.get_ylim()
            scale_integral = max_of_fit / max(integral_of_spectrum)
            plot(
                fields,
                integral_of_spectrum * scale_integral,
                "k:",
                alpha=0.5,
                linewidth=1,
                label="%0.2g $*\int \int$" % scale_integral,
            )
            legend(
                bbox_to_anchor=(1.05, 0, 0.5, 1),  # bounding box l b w h
                loc=2,  # upper left (of the bounding box)
                borderaxespad=0.0,
            )
            ax.set_ylim(ylims)
            setp(ax.get_xticklabels(), rotation=45)
            if show_components:  # if showing the components, also show rms
                rms = mean((fit / normalization - experimental / normalization) ** 2)
                ax.text(
                    0.75,
                    0.75,
                    "rms = %0.2g" % rms,
                    horizontalalignment="center",
                    verticalalignment="top",
                    transform=ax.transAxes,
                )

        show_spc(thisdatafile + ".spc")
        # show_spc(thisdatafile+'.spc_ref',linestyles = ['g--','b--'],show_components = False)
    show()
