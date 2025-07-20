#!/usr/bin/python
import os
import numpy as np
import nlsl


def read_column_data(filename):
    with open(filename, "r") as fp:
        data = [line.split() for line in fp]
    return np.array(data, dtype=np.double)


def run_example(example, allowed_rel_rms=None):
    """Run the numbered NLSL example and return the relative RMS error."""
    print(f"about to run nlsl example {example}")
    examples_dir = os.path.join(os.path.dirname(__file__), os.pardir, "examples")
    os.chdir(examples_dir)

    filename_base = f"sampl{example}"
    data_files_out = []

    def run_file(thisfp):
        for thisline in thisfp.readlines():
            if thisline[:5] == "call ":
                fp_called = open(thisline[5:].strip())
                run_file(fp_called)
                fp_called.close()
            elif thisline[:5] == "data ":
                nlsl.procline(thisline)
                data_files_out.append(thisline[5:].strip().split(' ')[0])
            else:
                nlsl.procline(thisline)
        thisfp.close()

    nlsl.nlsinit()
    run_file(open(filename_base + '.run'))

    rms_sq_total = 0.0
    exp_sq_total = 0.0
    for thisdatafile in data_files_out:
        data_calc = read_column_data(thisdatafile + '.spc')
        exp_sq_total += np.sum(data_calc[:, 1] ** 2)
        rms_sq_total += np.sum((data_calc[:, 2] - data_calc[:, 1]) ** 2)

    if exp_sq_total > 0:
        relative_rms = np.sqrt(rms_sq_total) / np.sqrt(exp_sq_total)
        if allowed_rel_rms is not None:
            assert relative_rms < allowed_rel_rms * 1.01, (
                'rms error / norm(experimental) = %0.3g' % relative_rms
            )
        return relative_rms


if __name__ == "__main__":
    run_example(1)
