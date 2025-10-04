import pytest
import numpy as np
import os
import nlsl

EXAMPLES = [
    (1, [0.0404]),
    (2, [0.0331, 0.0513]),
    (3, [0.06113]),
    (4, [0.04001]),
    (5, [0.075, 0.1592]),
]

def read_column_data(filename):
    with open(filename, "r") as fp:
        data = [line.split() for line in fp]
    return np.array(data, dtype=np.double)


def run_example(example, allowed_rel_rms=None):
    """Run the numbered NLSL example and return list of relative RMS errors."""

    runfile_location  = os.path.dirname(__file__)
    print(f"about to run nlsl example {example} in location {runfile_location}")
    os.chdir(runfile_location)

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
    n.write_spc()

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
                f'rms error / norm(experimental) = {rms}, but only {allowed*1.01} allowed'
            )
    return rel_rms_list

@pytest.mark.parametrize("example,allowed", EXAMPLES)
def test_runexample(example, allowed):
    rel_rms = run_example(example, allowed_rel_rms=allowed)
    # the following error should never trigger (b/c it's in run example), but just in case
    assert rel_rms and all(
        r < a * 1.01 for r, a in zip(rel_rms, allowed)
    ), f"I was expecting errors of {allowed}, but got errors of {rel_rms}"
