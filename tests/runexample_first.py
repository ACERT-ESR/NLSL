#!/usr/bin/python
import os
import numpy as np
import nlsl


def read_column_data(filename):
    with open(filename, 'r') as fp:
        data = [line.split() for line in fp]
    return np.array(data, dtype=np.double)


if __name__ == "__main__":
    print("about to run nlsl example 1")
    # Ensure files are read relative to this directory
    os.chdir(os.path.dirname(__file__))

    filename_base = 'sampl1'
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
        print('rms error / norm(experimental) = %0.3g' % relative_rms)
