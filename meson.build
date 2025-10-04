project('NLSL', 'c', 'fortran', version: '0.1.0')

py = import('python').find_installation(pure: false)
numpy_inc = run_command(py, '-c', 'import numpy,sys;sys.stdout.write(numpy.get_include())').stdout().strip()
fortranobject_src = run_command(py, '-c', 'import numpy,os,sys;import os;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src","fortranobject.c"))').stdout().strip()
numpy_f2py_src = run_command(py, '-c', 'import numpy,os,sys;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src"))').stdout().strip()
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
python_dep = py.dependency()

# Fortran sources
fortran_sources = files(
    'src/fortran/addprm.f90',
    'src/fortran/anxlk.f90',
    'src/fortran/assgnc.f90',
    'src/fortran/basis.f90',
    'src/fortran/basisc.f90',
    'src/fortran/bessel.f90',
    'src/fortran/bincom.f90',
    'src/fortran/brent.f90',
    'src/fortran/caldlm.f90',
    'src/fortran/ccrints.f90',
    'src/fortran/cd2km.f90',
    'src/fortran/cfs.f90',
    'src/fortran/cgltri.f90',
    'src/fortran/confc.f90',
    'src/fortran/convtc.f90',
    'src/fortran/correl.f90',
    'src/fortran/covar.f90',
    'src/fortran/covrpt.f90',
    'src/fortran/cscg.f90',
    'src/fortran/datac.f90',
    'src/fortran/daxpy.f90',
    'src/fortran/dchex.f90',
    'src/fortran/dcopy.f90',
    'src/fortran/ddot.f90',
    'src/fortran/dfunc.f90',
    'src/fortran/drotg.f90',
    'src/fortran/enorm.f90',
    'src/fortran/eprls.f90',
    'src/fortran/eprmat.f90',
    'src/fortran/eprprm.f90',
    'src/fortran/errmsg.f90',
    'src/fortran/expdat.f90',
    'src/fortran/fitc.f90',
    'src/fortran/fitl.f90',
    'src/fortran/fixc.f90',
    'src/fortran/fmomnt.f90',
    'src/fortran/ftfuns.f90',
    'src/fortran/ftwork.f90',
    'src/fortran/fz.f90',
    'src/fortran/gconvl.f90',
    'src/fortran/getdat.f90',
    'src/fortran/helpc.f90',
    'src/fortran/ipar.f90',
    'src/fortran/ipsfind.f90',
    'src/fortran/iterat.f90',
    'src/fortran/l1pfun.f90',
    'src/fortran/lbasix.f90',
    'src/fortran/lcheck.f90',
    'src/fortran/letc.f90',
    'src/fortran/lfun.f90',
    'src/fortran/lgrint.f90',
    'src/fortran/lmcom.f90',
    'src/fortran/lmnls.f90',
    'src/fortran/lmpar.f90',
    'src/fortran/lpnam.f90',
    'src/fortran/matrll.f90',
    'src/fortran/maxl.f90',
    'src/fortran/mnbrak.f90',
    'src/fortran/momdls.f90',
    'src/fortran/mspctr.f90',
    'src/fortran/mtsdef.f90',
    'src/fortran/nlsdim.f90',
    'src/fortran/nlsl.f90',
    'src/fortran/nlsnam.f90',
    'src/fortran/ordrpr.f90',
    'src/fortran/parc.f90',
    'src/fortran/parcom.f90',
    'src/fortran/parsav.f90',
    'src/fortran/physcn.f90',
    'src/fortran/pidef.f90',
    'src/fortran/plgndr.f90',
    'src/fortran/pmatrl.f90',
    'src/fortran/pstvec.f90',
    'src/fortran/qrfac.f90',
    'src/fortran/qrsolv.f90',
    'src/fortran/qrutil.f90',
    'src/fortran/rmvprm.f90',
    'src/fortran/rnddbl.f90',
    'src/fortran/scalec.f90',
    'src/fortran/scmvm.f90',
    'src/fortran/series.f90',
    'src/fortran/setmts.f90',
    'src/fortran/setnm.f90',
    'src/fortran/setprm.f90',
    'src/fortran/setspc.f90',
    'src/fortran/shiftc.f90',
    'src/fortran/sitec.f90',
    'src/fortran/srchc.f90',
    'src/fortran/sscale.f90',
    'src/fortran/sshift.f90',
    'src/fortran/statc.f90',
    'src/fortran/stats.f90',
    'src/fortran/stdio.f90',
    'src/fortran/strutl1.f90',
    'src/fortran/strutl2.f90',
    'src/fortran/stvect.f90',
    'src/fortran/symdef.f90',
    'src/fortran/tdchek.f90',
    'src/fortran/tdsqz.f90',
    'src/fortran/tensym.f90',
    'src/fortran/timer.f90',
    'src/fortran/tridag.f90',
    'src/fortran/varyc.f90',
    'src/fortran/w3j.f90',
    'src/fortran/writec.f90',
    'src/fortran/writr.f90',
    'src/fortran/zaxpy.f90',
    'src/fortran/zaypx.f90',
    'src/fortran/zdotu.f90'
)

# generate f2py wrappers
c_wrappers = custom_target(
    'f2py_wrappers',
    # interface blocks require the modules containing fit variables
    input: [
        'pynlsl.pyf',
        'src/fortran/eprprm.f90',
        'src/fortran/lmcom.f90',
        'src/fortran/parcom.f90',
        'src/fortran/lpnam.f90',
        'src/fortran/expdat.f90',
        'src/fortran/iterat.f90',
        'src/fortran/ipsfind.f90',
        'src/fortran/fitl.f90',
        'src/fortran/lfun.f90',
        'src/fortran/stdio.f90',
        'src/fortran/mspctr.f90',
        'src/fortran/nlsl.f90',
        'src/fortran/setprm.f90',
        'src/fortran/setspc.f90',
        'src/fortran/writr.f90',
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

py.install_sources('nlsl/__init__.py', subdir: 'nlsl')

# Install the Python test suite alongside the package so editable installs
# can locate the extension module without modifying sys.path in tests.
install_subdir('tests',
    install_dir: join_paths(py.get_install_dir(), 'nlsl', 'tests'))
