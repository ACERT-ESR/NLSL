# Pytest Full Suite Results

Running `pytest` after rebuilding the editable `nlsl` module with the system `ninja` executable currently fails. The failure stems from missing Fortran symbol `ipfind` in the `nlsl.fortrancore` module, which triggers 20 alias-related tests and one sample test to error out. See the latest `pytest` run for details.
