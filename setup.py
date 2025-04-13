import logging
from datetime import datetime
import os
import shutil
import sys
import subprocess
import tempfile
from setuptools import setup, Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler


# Formatter for the logger
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)

# Install packages from pip ==============================================================
def install_with_pip(pack, vers=None, log=None, namepkg=None):

    # Update pip
    p = subprocess.Popen([sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()

    # sys.executable gives the path of the python interpreter
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if vers is None:
        m = "{}: ** {}: Installing {}".format(now, namepkg, pack)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}".format(pack)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
    else:
        m = "{}: ** {}: Installing {}=={}".format(now, namepkg, pack, vers)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers), " &>install.log"])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()


def is_msd_fftw3_installed(log=None, namepkg=None):

    """
    Installing the Mean Square Displacement program from github.
    This function compiles and install the msd program under ~/.local/bin

    This code has been forked from:
        Raul P. Pelaez 2019. Fast Mean Square Displacement
        This code computes the MSD of a list of trajectories in spunto format using
        the FFT in O(N) time. It can be compiled/run with or without GPU support.
    """

    import git

    giturl = "https://github.com/jrdcasa/MeanSquareDisplacement.git"
    start_dir = os.getcwd()
    install_dir = 'polyanagro/ext_libc/MeanSquareDisplacement'
    executable_file_compiled = 'polyanagro/ext_libc/MeanSquareDisplacement/build/bin/msd'
    executable_file_installed = os.path.expanduser('~/.local/bin/msd')

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    msg = "\n\t\t INSTALLING MSD_FFTW3 PROGRAM INTO THE PYTHON ENVIRONMENT FROM GITHUB\n\n"

    if os.path.isfile(executable_file_installed):
        msg += "{}: ** {}: MSD_FFTW3 program is already installed in your system. {}\n".format(now, namepkg, giturl)
        msg += "{}".format(os.path.join(os.getcwd(), executable_file_installed))
        print(msg) if log is None else log.info(msg)
    else:
        msg += "{}: ** {}: MSD_FFTW3 program is not installed in your system\n".format(now, namepkg)
        msg += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        msg += "{}".format(os.path.join(os.getcwd(), install_dir))
        print(msg) if log is None else log.info(msg)

        fullpath_install = os.path.join(os.path.abspath(install_dir))

        # Look at thirdparty directory
        if os.path.isdir(fullpath_install):
            shutil.rmtree(fullpath_install)
        os.makedirs(fullpath_install)

        try:
            git.Repo.clone_from(giturl, fullpath_install)
        except git.GitCommandError as e:
            msg = "================= ERROR INSTALL ================\n"
            msg += "** {}: The github repository for MSD_FFTW3 program is not valid or not exists.!!!\n".format(namepkg)
            msg += "** {}: giturl     : {}\n".format(namepkg, giturl)
            msg += "** {}: install_dir: {}\n".format(namepkg, fullpath_install)
            msg += "** {}: MSD_FFTW3 program cannot be installed\n".format(namepkg)
            msg += "** {}: The installation is aborted\n".format(namepkg)
            msg += "** Error: {}\n".format(e.stderr)
            msg += "\n================= ERROR INSTALL ================"
            print(msg) if log is None else log.info(msg)
            exit()

        os.chdir(fullpath_install)
        os.makedirs(os.path.join(fullpath_install, "build"))
        os.chdir("./build")

        command = ["cmake", "..", "-DDONT_USE_CUDA=ON", "-DUSE_BOOST=OFF",  "-DUSE_CPU=ON" ,"-DUSE_MKL=OFF"]
        subprocess.run(command, check=True)
        command = ["make"]
        subprocess.run(command, check=True)
        command = ["cp", "bin/msd", os.path.expanduser("~/.local/bin/msd")]
        subprocess.run(command, stderr=subprocess.PIPE, text=True)

        if os.path.isfile(os.path.join(os.getcwd(), executable_file_installed)):
            msg += "{}: ** {}: MSD_FFTW3 program has been installed in your system. {}\n".format(now, namepkg, giturl)
            msg += "{}".format(os.path.join(os.getcwd(), executable_file_installed))
            print(msg) if log is None else log.info(msg)
        else:
            msg = "================= ERROR INSTALL ================\n"
            msg += "** {}: giturl     : {}\n".format(namepkg, giturl)
            msg += "** {}: install_dir: {}\n".format(namepkg, fullpath_install)
            msg += "** {}: MSD_FFTW3 program has not been installed\n".format(namepkg)
            msg += "** {}: The installation is aborted\n".format(namepkg)
            msg += "** File does not exists: {}\n".format(executable_file_installed)
            msg += "\n================= ERROR INSTALL ================"
            print(msg) if log is None else log.info(msg)
            exit()

    os.chdir(start_dir)

# ================================================================================================
def install_topology_library(log=None, namepkg=None):
    """
    Installing the python topology library if is not present in the python environment.
    """

    import git

    giturl = "https://github.com/jrdcasa/topology.git"
    install_dir = 'topology'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    msg = "\n\t\t INSTALLING TOPOLOGY LIBRARY INTO THE PYTHON ENVIRONMENT FROM GITHUB\n\n"

    if os.path.isdir(os.path.join(os.getcwd(), install_dir)):
        msg += "{}: ** {}: topology library is already installed in your system. {}".format(now, namepkg, giturl)
        print(msg) if log is None else log.info(msg)
    else:
        msg += "{}: ** {}: topology library is not installed in your system\n".format(now, namepkg)
        msg += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        print(msg) if log is None else log.info(msg)

        fullpath_install = os.path.abspath(install_dir)

        # Look at thirdparty directory
        if os.path.isdir(fullpath_install):
            shutil.rmtree(fullpath_install)
        os.makedirs(fullpath_install)

        try:
            git.Repo.clone_from(giturl, fullpath_install)
        except git.GitCommandError as e:
            msg = "================= ERROR INSTALL ================\n"
            msg += "** {}: The github repository for topology is not valid or not exists.!!!\n".format(namepkg)
            msg += "** {}: giturl     : {}\n".format(namepkg, giturl)
            msg += "** {}: install_dir: {}\n".format(namepkg, fullpath_install)
            msg += "** {}: Topology library cannot be installed\n".format(namepkg)
            msg += "** {}: The installation is aborted\n".format(namepkg)
            msg += "** Error: {}\n".format(e.stderr)
            msg += "\n================= ERROR INSTALL ================"
            print(msg) if log is None else log.info(msg)
            exit()

        os.chdir(fullpath_install)
        #print(os.getcwd())
        subprocess.call(["python", "setup.py", "install"])
        os.chdir("..")

# Disabling-output-when-compiling-with-distutil =================================================
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as fout:
                if include is not None:
                    fout.write('#include {0!s}\n'.format(include))
                fout.write('int main(void) {\n')
                fout.write('    {0!s};\n'.format(funcname))
                fout.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)

# Does this compiler support OpenMP parallelization?""" ==============================================================
def detect_openmp():
    print("POLYANAGRO: Attempting to autodetect OpenMP support... ")
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("POLYANAGRO: Compiler supports OpenMP")
    else:
        print("POLYANAGRO: Did not detect OpenMP support.")

    return hasopenmp

# Setup external extensions ==============================================================
def setup_external_extensions(debug_cflags=False, use_openmp=True):
    has_openmp = detect_openmp()

    # parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    mathlib = ['m']
    define_macros = []
    extra_compile_args = ['-std=c99', '-ffast-math', '-O3', '-funroll-loops', '-Wno-cpp']
    if debug_cflags:
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    extensions_install = [
        Extension("ext_libc.c_rg_openmp", ["polyanagro/ext_libc/c_rg_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),
        Extension("ext_libc.c_unit_bond_vectors", ["polyanagro/ext_libc/c_unit_bond_vectors.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),
        Extension("ext_libc.c_rdf_openmp", ["polyanagro/ext_libc/c_rdf_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=['-lgomp']),
        Extension("ext_libc.c_distC", ["polyanagro/ext_libc/c_distC.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=['-lgomp']),
        Extension("ext_libc.c_acf_openmp", ["polyanagro/ext_libc/c_acf_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=['-lgomp']),
        Extension("ext_libc.c_internal_distances", ["polyanagro/ext_libc/c_internal_distances.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  ),
        Extension("ext_libc.c_msd_openmp", ["polyanagro/ext_libc/c_msd_openmp.pyx"],
                  libraries=mathlib + parallel_libraries + ["fftw3"],
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=["-lfftw3", "-lm"]
                  ),
    ]

    return extensions_install


# Requeriments to be manually installed  ===========================================================================
def check_requirements_outside():
    # Check for swig (http://www.swig.org)
    if os.system("which swig"):
        msg = "ERROR. Please install SWIG in your system (http://www.swig.org)\n"
        msg += "ERROR. Ubuntu: apt get install swig"
        print(msg) if logger is None else logger.info(msg)
        exit()

    # Check for cmake
    if os.system("which cmake"):
        msg = "ERROR. Please install CMAKE in your system\n"
        msg += "ERROR. Ubuntu: apt get install cmake"
        print(msg) if logger is None else logger.info(msg)
        exit()

    # Check for pygraphviz
    try:
        import pygraphviz as pag
    except ImportError:
        msg = "ERROR. Please install PYGRAPHVIZ in your system"

        print(msg) if logger is None else logger.info(msg)
        exit()


# Main setup
if __name__ == '__main__':

    # Creating the logger to install.log file ===================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    # SWIG, cmake and pygraphviz are needed for topology library
    check_requirements_outside()

    # Print sys path ===================================
    m1 = "\t\t SYS PATH\n"
    for item in sys.path:
        m1 += item + "\n"
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 += "\n\t\t INSTALLING PIP PACKAGES ({})\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    # Install requirements ===================================
    # with open('requirements.txt') as f:
    #     required = f.read().splitlines()
    # for ipack in required:
    #     try:
    #         pkg, version = ipack.split(">=")[0:2]
    #         if pkg[0] == "#":
    #             continue
    #         install_with_pip(pkg, vers=version, log=logger)
    #     except ValueError:
    #         pkg = ipack
    #         if pkg[0] == "#" or len(pkg)<2:
    #             continue
    #         install_with_pip(pkg, log=logger)
    #     finally:
    #         pass

    # Check MSD installation
    # is_msd_fftw3_installed(logger, namepkg="MSD_FFTW3")

    # Check for topology installation
    try:
        import topology
    except ImportError:
        m = "ERROR. Topology library is not correctly installed\n"
        m +="ERROR. Please check ./topology/install.log"
        print(m) if logger is None else logger.info(m)
        exit()

    # # Setup POLYANAGRO ===========================================
    from Cython.Build import cythonize
    import numpy

    # Extensions
    extensions = setup_external_extensions()
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t RUNNING SETUP FROM SETUPTOOLS {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    print(os.getcwd())
    setup(
        ext_modules=cythonize(extensions, compiler_directives={'language_level': sys.version_info[0]}),
        include_dirs=[numpy.get_include()],
        )

