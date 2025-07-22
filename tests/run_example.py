#!/usr/bin/python
import os
import numpy as np
import nlsl


def read_column_data(filename):
    with open(filename, "r") as fp:
        data = [line.split() for line in fp]
    return np.array(data, dtype=np.double)


def run_example(example, allowed_rel_rms=None):
    """Run the numbered NLSL example and return list of relative RMS errors."""

    print(f"about to run nlsl example {example}")
    examples_dir = os.path.join(os.path.dirname(__file__), os.pardir, "examples")
    os.chdir(examples_dir)

    filename_base = f"sampl{example}"
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
                data_files_out.append(thisline[5:].strip().split(' ')[0])
            else:
                n.procline(thisline)
        thisfp.close()

    run_file(open(filename_base + '.run'))

    rel_rms_list = []
    for thisdatafile in data_files_out:
        data_calc = read_column_data(thisdatafile + '.spc')
        exp_sq = np.sum(data_calc[:, 1] ** 2)
        rms_sq = np.sum((data_calc[:, 2] - data_calc[:, 1]) ** 2)
        if exp_sq > 0:
            rel_rms = np.sqrt(rms_sq) / np.sqrt(exp_sq)
            rel_rms_list.append(rel_rms)

    if allowed_rel_rms is not None and rel_rms_list:
        assert len(rel_rms_list) == len(allowed_rel_rms)
        for rms, allowed in zip(rel_rms_list, allowed_rel_rms):
            assert rms < allowed * 1.01, (
                'rms error / norm(experimental) = %0.3g' % rms
            )
    return rel_rms_list


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run an NLSL example")
    parser.add_argument("example", type=int, nargs="?", default=1,
                        help="example number to run")
    parser.add_argument("--allowed-rel-rms", type=float, nargs="*", dest="allowed",
                        default=None,
                        help="fail if relative RMS exceeds these values")
    args = parser.parse_args()

    rms_list = run_example(args.example, allowed_rel_rms=args.allowed)
    if rms_list:
        for i, rms in enumerate(rms_list, 1):
            print(f"spectrum {i}: relative rms = {rms:.5g}")
