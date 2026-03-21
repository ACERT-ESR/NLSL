import pytest
import numpy as np
import os
from pathlib import Path
import nlsl

from .runfile_helpers import run_runfile

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

    runfile_path = Path(__file__).with_name(f"sampl{example}.run")
    data_files_out = []
    n = nlsl.nlsl()

    def record_commands(command):
        if command.lower().startswith("data "):
            data_files_out.append(command[5:].strip().split(' ')[0])

    run_runfile(n, runfile_path, command_callback=record_commands)
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
