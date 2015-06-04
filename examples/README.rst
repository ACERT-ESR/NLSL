        NLSL test files
====================================

The `.run` files `sampl1` through `sampl5` give examples of use of the
program nlsl in typical experimental situations.  They read in the
`.dat` files, and the results of the fits are stored in the `.log`
files and the `.spc` files which may be compared with the original
`.dat` files, or with the `.*_ref` files.

The original usage was:

- Run all the runfiles (assuming you are in the "examples" directory):
  - one usage:

    ../nlsl
    NLSL> read sampl1.run
    "lots of output here"
    NLSL> quit
  - another usage:

    ../nlsl < sampl1.run

- Check the .spc and .log files created against the ones contained in
  this directory.

