project('NLSL', 'c', 'fortran', version: '0.1.0')

py = import('python').find_installation(pure: false)
numpy_inc = run_command(py, '-c', 'import numpy,sys;sys.stdout.write(numpy.get_include())').stdout().strip()
fortranobject_src = run_command(py, '-c', 'import numpy,os,sys;import os;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src","fortranobject.c"))').stdout().strip()
numpy_f2py_src = run_command(py, '-c', 'import numpy,os,sys;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src"))').stdout().strip()
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
python_dep = py.dependency()

# Fortran sources grouped by subsystem
cli_sources = files(
    'src/fortran/cli/addprm.f90',
    'src/fortran/cli/assgnc.f90',
    'src/fortran/cli/confc.f90',
    'src/fortran/cli/convtc.f90',
    'src/fortran/cli/datac.f90',
    'src/fortran/cli/fitc.f90',
    'src/fortran/cli/fixc.f90',
    'src/fortran/cli/helpc.f90',
    'src/fortran/cli/letc.f90',
    'src/fortran/cli/nlsl.f90',
    'src/fortran/cli/ordrpr.f90',
    'src/fortran/cli/parc.f90',
    'src/fortran/cli/rmvprm.f90',
    'src/fortran/cli/scalec.f90',
    'src/fortran/cli/series.f90',
    'src/fortran/cli/setprm.f90',
    'src/fortran/cli/shiftc.f90',
    'src/fortran/cli/sitec.f90',
    'src/fortran/cli/srchc.f90',
    'src/fortran/cli/statc.f90',
    'src/fortran/cli/stdio.f90',
    'src/fortran/cli/strutl1.f90',
    'src/fortran/cli/strutl2.f90',
    'src/fortran/cli/varyc.f90',
)

least_squares_sources = files(
    'src/fortran/least_squares/brent.f90',
    'src/fortran/least_squares/caldlm.f90',
    'src/fortran/least_squares/covar.f90',
    'src/fortran/least_squares/covrpt.f90',
    'src/fortran/least_squares/fitl.f90',
    'src/fortran/least_squares/lmnls.f90',
    'src/fortran/least_squares/lmpar.f90',
    'src/fortran/least_squares/maxl.f90',
    'src/fortran/least_squares/mnbrak.f90',
    'src/fortran/least_squares/qrsolv.f90',
    'src/fortran/least_squares/sscale.f90',
    'src/fortran/least_squares/sshift.f90',
    'src/fortran/least_squares/stats.f90',
)

matrix_model_sources = files(
    'src/fortran/matrix_model/anxlk.f90',
    'src/fortran/matrix_model/basis.f90',
    'src/fortran/matrix_model/basisc.f90',
    'src/fortran/matrix_model/bessel.f90',
    'src/fortran/matrix_model/bincom.f90',
    'src/fortran/matrix_model/cd2km.f90',
    'src/fortran/matrix_model/dfunc.f90',
    'src/fortran/matrix_model/fmomnt.f90',
    'src/fortran/matrix_model/l1pfun.f90',
    'src/fortran/matrix_model/lbasix.f90',
    'src/fortran/matrix_model/lcheck.f90',
    'src/fortran/matrix_model/lfun.f90',
    'src/fortran/matrix_model/lgrint.f90',
    'src/fortran/matrix_model/matrll.f90',
    'src/fortran/matrix_model/momdls.f90',
    'src/fortran/matrix_model/pmatrl.f90',
    'src/fortran/matrix_model/setmts.f90',
    'src/fortran/matrix_model/symdef.f90',
    'src/fortran/matrix_model/tdchek.f90',
    'src/fortran/matrix_model/tdsqz.f90',
    'src/fortran/matrix_model/tensym.f90',
)

krylov_solver_sources = files(
    'src/fortran/krylov_solvers/ccrints.f90',
    'src/fortran/krylov_solvers/cfs.f90',
    'src/fortran/krylov_solvers/cgltri.f90',
    'src/fortran/krylov_solvers/cscg.f90',
    'src/fortran/krylov_solvers/eprls.f90',
    'src/fortran/krylov_solvers/eprmat.f90',
    'src/fortran/krylov_solvers/eprprm.f90',
    'src/fortran/krylov_solvers/fz.f90',
    'src/fortran/krylov_solvers/pstvec.f90',
    'src/fortran/krylov_solvers/scmvm.f90',
    'src/fortran/krylov_solvers/stvect.f90',
    'src/fortran/krylov_solvers/tridag.f90',
)

io_preprocessing_sources = files(
    'src/fortran/io_preprocessing/correl.f90',
    'src/fortran/io_preprocessing/gconvl.f90',
    'src/fortran/io_preprocessing/getdat.f90',
    'src/fortran/io_preprocessing/setnm.f90',
    'src/fortran/io_preprocessing/setspc.f90',
    'src/fortran/io_preprocessing/writec.f90',
    'src/fortran/io_preprocessing/writr.f90',
)

core_sources = files(
    'src/fortran/core/errmsg.f90',
    'src/fortran/core/expdat.f90',
    'src/fortran/core/ipar.f90',
    'src/fortran/core/ipsfind.f90',
    'src/fortran/core/iterat.f90',
    'src/fortran/core/lmcom.f90',
    'src/fortran/core/lpnam.f90',
    'src/fortran/core/mspctr.f90',
    'src/fortran/core/mtsdef.f90',
    'src/fortran/core/nlsdim.f90',
    'src/fortran/core/nlsnam.f90',
    'src/fortran/core/parcom.f90',
    'src/fortran/core/parsav.f90',
    'src/fortran/core/physcn.f90',
    'src/fortran/core/pidef.f90',
    'src/fortran/core/rnddbl.f90',
    'src/fortran/core/timer.f90',
)

numerics_sources = files(
    'src/fortran/numerics/daxpy.f90',
    'src/fortran/numerics/dchex.f90',
    'src/fortran/numerics/dcopy.f90',
    'src/fortran/numerics/ddot.f90',
    'src/fortran/numerics/drotg.f90',
    'src/fortran/numerics/enorm.f90',
    'src/fortran/numerics/ftfuns.f90',
    'src/fortran/numerics/ftwork.f90',
    'src/fortran/numerics/plgndr.f90',
    'src/fortran/numerics/qrfac.f90',
    'src/fortran/numerics/qrutil.f90',
    'src/fortran/numerics/w3j.f90',
    'src/fortran/numerics/zaxpy.f90',
    'src/fortran/numerics/zaypx.f90',
    'src/fortran/numerics/zdotu.f90',
)

fortran_sources = cli_sources + least_squares_sources + matrix_model_sources + \
    krylov_solver_sources + io_preprocessing_sources + core_sources + numerics_sources

# generate f2py wrappers
c_wrappers = custom_target(
    'f2py_wrappers',
    # interface blocks require the modules containing fit variables
    input: [
        'pynlsl.pyf',
        'src/fortran/krylov_solvers/eprprm.f90',
        'src/fortran/core/lmcom.f90',
        'src/fortran/core/parcom.f90',
        'src/fortran/core/lpnam.f90',
        'src/fortran/core/expdat.f90',
        'src/fortran/core/iterat.f90',
        'src/fortran/core/ipsfind.f90',
        'src/fortran/least_squares/fitl.f90',
        'src/fortran/matrix_model/lfun.f90',
        'src/fortran/cli/stdio.f90',
        'src/fortran/cli/addprm.f90',
        'src/fortran/core/mspctr.f90',
        'src/fortran/cli/nlsl.f90',
        'src/fortran/cli/setprm.f90',
        'src/fortran/cli/rmvprm.f90',
        'src/fortran/io_preprocessing/setspc.f90',
        'src/fortran/io_preprocessing/writr.f90',
    ],
    output: ['fortrancoremodule.c', 'fortrancore-f2pywrappers.f', 'fortrancore-f2pywrappers2.f90'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fortrancore', '--lower'],
    build_by_default: true,
)

fixed_lib = static_library(
    'nlsl_fixed',
    fortran_sources,
    fortran_args: ['-O2','-g','-fopenmp','-std=gnu','-ffixed-form','-ffixed-line-length-none'],
    pic: true,
    install: false,
)

mod_dir = join_paths(meson.current_build_dir(), 'libnlsl_fixed.a.p')

wrapper_lib = static_library(
    'nlsl_wrapper',
    [c_wrappers[1], c_wrappers[2]],
    fortran_args: ['-O2','-g','-fopenmp','-std=gnu', '-I' + mod_dir],
    dependencies: [python_dep],
    pic: true,
    install: false,
)

c_sources = files(
    'catch.c',
    'genio.c',
    'lprmpt.c',
    'pltx_dummy.c',
    fortranobject_src,
)

c_lib = static_library(
    'nlsl_c',
    c_sources + [c_wrappers[0]],
    c_args: ['-DADD_UNDERSCORE', '-I' + numpy_inc, '-I' + numpy_f2py_src],
    dependencies: [python_dep],
    pic: true,
    install: false,
)

py.extension_module(
    'fortrancore',
    [],
    objects: [
        fixed_lib.extract_all_objects(),
        wrapper_lib.extract_all_objects(),
        c_lib.extract_all_objects(),
    ],
    c_args: ['-I'+numpy_inc, '-I'+numpy_f2py_src],
    link_args: ['-fopenmp'],
    subdir: 'nlsl',
    install: true,
)

py.install_sources('nlsl/__init__.py', 'nlsl/data.py', 'nlsl/cli.py', subdir: 'nlsl')

# Install the Python test suite alongside the package so editable installs
# can locate the extension module without modifying sys.path in tests.
install_subdir('tests',
    install_dir: join_paths(py.get_install_dir(), 'nlsl', 'tests'))
