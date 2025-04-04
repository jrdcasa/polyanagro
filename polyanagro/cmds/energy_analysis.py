import argparse
import os
import utils
import sys
import datetime
from collections import defaultdict
import pandas as pd
from subprocess import PIPE, Popen
from polyanagro.Energy import energy_analysis
from polyanagro.Custom_Plots import Custom_Plots


# ************************************
def call_subprocess(cmd, cwd=None):
    """
    Execute shell command `cmd`, using the environment `env`
    and wait for the results.
    """
    try:
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd)
        rs = p.communicate()
    except Exception as e:
        msg1 = "Subprocess fails with error: {}".format(e)
        msg2 = "Command: {}\n".format(cmd)
        print(msg1 + msg2)

    return rs

# =============================================================================
def parse_arguments():

    # ************************************
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    # ************************************
    def join_energy_files():

        # Check that all files have the same extensions
        if len(ext_set) > 1:
            print("\nERROR: All energy files must have the same extension:")
            for item in args.energy_list:
                print("\t"+item)
            exit()
        # For list of edrs we need the path to gromacs
        if mdpackage == "GROMACS" and args.joinpath is None:
            print("A list of energy files requires to kwown the PATH to GROMACS in your system.")
            print("Use the option --joinpath <PATH_TO_GROMACS_GMX>")
            print("Example: --joinpath /usr/bin/gmx")
            exit()
        else:
            workdir = os.getcwd()
            outTpr = os.path.join(workdir, 'united_edr.edr')
            if os.path.isfile(outTpr):
                cmd = ["rm", outTpr]
                call_subprocess(cmd, cwd=workdir)
            cmd = [args.joinpath, "eneconv", "-f"]
            for item in args.energy_list:
                cmd.append(item)
            cmd.append('-o')
            cmd.append(outTpr)
            call_subprocess(cmd, cwd=workdir)

        return outTpr

    # ****** MAIN ARG PARSER ******
    desc = """Energy analysis.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)
    subparser = parser.add_subparsers(dest='command', required=True)

    info = subparser.add_parser('info')
    calc = subparser.add_parser('calc')

    info.add_argument("-e", "--energy", dest="energy_list", nargs="+",
                        help="Energy file from MD package.\n "
                             "The package is detected by the extension of the file",
                        action="store", required=True, default=None)

    calc.add_argument("-e", "--energy", dest="energy_list", nargs="+",
                        help="Energy file from MD package.\n "
                             "The package is detected by the extension of the file",
                        action="store", required=True, default=None)

    info.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="energy_analysis_info.log")

    calc.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="energy_analysis.log")

    calc.add_argument("--tbegin", dest="tbegin", type=float,
                        help="Starting time in ps to perform the analysis.\n"
                             "Example: A value of 10, start the analysis at 10ps.",
                        action="store", required=False, default=None)

    calc.add_argument("--tend", dest="tend", type=float,
                        help="Ending time in ps to perform the analysis.\n"
                             "Example: A value of 30, end the analysis at 30ps.",
                        action="store", required=False, default=None)

    calc.add_argument("--joinpath", dest="joinpath", type=str,
                        help="Path to the program to join energy files.\n"
                             "Example: For GROMACS --> /usr/bin/gmx",
                        action="store", required=False, default=None)

    calc.add_argument("--groupterms", dest="groupterm",
                        help="Group energy terms (bond, nonbond, ...).",
                        action="store_true", required=False)

    calc.add_argument("--avg", dest="avg",
                        help="Calculate averages from <begin> to <end>.",
                        action="store_true", required=False)

    calc.add_argument("--acf", dest="acf_list", nargs="+",
                        help="Autocorrelation function (ACF) of the time series.\n "
                             "A list with the labels of the parameter to calculate the ACF",
                        action="store", required=False, default=None)


    args = parser.parse_args()

    # Check the extension of the energy files
    ext_list = []
    for iener in args.energy_list:
        ext_list.append(os.path.splitext(iener)[1])
        if not os.path.isfile(iener):
            print("\nERROR: File {} does not exist\n".format(iener))
            exit()
    ext_set = set(ext_list)
    if ext_list[0] == ".edr":
        mdpackage = "GROMACS"
    elif ext_list[0] == ".log":
        mdpackage = "LAMMPS"
    else:
        mdpackage = None
        m = "\nERROR: Extension '{}' for energy files is unkwown.\n".format(ext_list[0])
        print(m)
        exit()

    if len(args.energy_list) > 1:
        args.energy = join_energy_files()
    else:
        args.energy = args.energy_list[0]

    return args

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                                Energy Analysis 
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This utility is part of the polyanagro library. Polyanagro is an 
        open-source python library to analyze simulations of polymer systems.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger_log is None else logger_log.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    m += "\t\t\t         or\n"
    m += "\t\t\tenergy_analysis".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)

# =============================================================================
def energy_calc_info_function(e, args, log):

    m = "\t\t**** ENERGY INFO ****\n"
    m += "\t\t  Number of points in energy files   : {}\n".format(len(e._df))
    ll = e._df['Time'].to_list()
    start_time = ll[0]
    end_time = ll[-1]
    if args.tbegin is None:
        args.tbegin = start_time
    if args.tend is None:
        args.tend = end_time
    m += "\t\t  Original trajectory from {} (ps) to {} (ps).\n".format(start_time, end_time)
    m += "\t\t  Analyzed trajectory from {} (ps) to {} (ps).\n".format(args.tbegin, args.tend)
    m += "\t\t  Timestep {} (ps).\n".format(ll[1]-ll[0])
    m += "\t\t**** End ENERGY INFO ****\n"
    print(m) if log is None else log.info(m)

    return start_time, end_time

# =============================================================================
def energy_info_function(e, args, log):

    m = "\t\t**** ENERGY INFO ****\n"
    m += "\t\t  Number of points in energy files   : {}\n".format(len(e._df))
    ll = e._df['Time'].to_list()
    start_time = ll[0]
    end_time = ll[-1]
    m += "\t\t  Original trajectory from {} (ps) to {} (ps).\n".format(start_time, end_time)
    m += "\t\t  Timestep {} (ps).\n".format(ll[1]-ll[0])
    m += "\t\t**** End ENERGY INFO ****\n"
    print(m) if log is None else log.info(m)

    m = "\t\t**** TERMS PRESENT IN THE ENERGY FILE ****\n"
    for iseries in e._df:
        label = iseries
        m += "\t\t  "+label+"\n"
    m += "\t\t**** END TERMS PRESENT IN THE ENERGY FILE ****\n"
    print(m) if log is None else log.info(m)

    return start_time, end_time

# =============================================================================
def individual_energy_figures_function(e, args, log):

    # INDIVIDUAL ENERGY FIGURES
    workdir = os.path.join(os.getcwd(), "individual_energy")
    if os.path.isdir(workdir):
        cmd = ["rm", "-r", workdir]
        call_subprocess(cmd, cwd=workdir)
    os.makedirs("individual_energy")
    if args.tbegin is not None:
        e.plot_energy_time(path_to_save=workdir, begin=args.tbegin)
    else:
        e.plot_energy_time(path_to_save=workdir, begin=0)
    m = "\t\t Individual energy plots have been saved in:\n"
    m += "\t\t  {}\n".format(workdir)
    print(m) if log is None else log.info(m)

# =============================================================================
def grouped_energy_figures_function(e, args, log):

    workdir = os.path.join(os.getcwd(), "grouped_energy")
    if os.path.isdir(workdir):
        cmd = ["rm", "-r", workdir]
        call_subprocess(cmd, cwd=workdir)
    os.makedirs("grouped_energy")
    e.plot_energy_group(skip_data = args.tbegin, path_to_save=workdir)
    workdir = os.path.join(os.getcwd(), "grouped_density")
    if os.path.isdir(workdir):
        cmd = ["rm", "-r", workdir]
        call_subprocess(cmd, cwd=workdir)
    os.makedirs("grouped_density")
    e.plot_density_group(skip_data = args.tbegin, path_to_save=workdir)
    m = "\t\t Grouped density plots have been saved in:\n"
    m += "\t\t  {}\n".format(workdir)
    print(m) if log is None else log.info(m)

# =============================================================================
def avg_energy_function(e, args, log):

    pd.options.display.float_format = '{:.5f}'.format
    e.avg_data(tbegin=args.tbegin, tend=args.tend)
    m = "\t\t**** ENERGY AVG INFO ****\n"
    print(m) if log is None else log.info(m)
    print(e._dfavg) if log is None else log.info(e._dfavg)
    print("\n") if log is None else log.info("\n")
    m = "\t\t**** END ENERGY AVG INFO ****\n"
    print(m) if log is None else log.info(m)

# =============================================================================
def acf_energy_function(e, args, end_time, log):

    dict_acf = defaultdict()
    workdir = os.path.join(os.getcwd(), "individual_acf")
    if os.path.isdir(workdir):
        cmd = ["rm", "-r", workdir]
        call_subprocess(cmd, cwd=workdir)
    os.makedirs("individual_acf")
    for ilabel in args.acf_list:
        ydata_acf = e.autocorrelation_data(ilabel, tbegin=0, tend=end_time)
        xdata = [i for i in range(0, len(ydata_acf))]
        p = Custom_Plots(figtitle=ilabel + "_ACF")
        p.simple_line_plot(xdata, ydata_acf, xlabel="lag (frames)",
                           ylabel=ilabel + "_acf", path_to_save=workdir)
        dict_acf[ilabel] = ydata_acf
    acf_data = pd.DataFrame.from_dict(dict_acf)
    filename = os.path.join(workdir, "acf_data.dat")
    acf_data.to_csv(filename, sep=' ')


# =============================================================================
def main_app():

    import polyanagro as pag

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)

    # All data in edr file is represented in a png file. The files
    # are stored in ./test01_a_figures
    file_energy = os.path.abspath(args.energy)
    e = energy_analysis(file_energy, logger=log)
    e.read_energy()

    # Print info of the trajectory
    if args.command == "info":

        start_time, end_time = energy_info_function(e, args, log)

    elif args.command == "calc":

        start_time, end_time = energy_calc_info_function(e, args, log)

        if args.tbegin >= end_time:
            m = "\n"
            print(m) if log is None else log.info(m)
            m = "\nThe number of skipped points cannot be greater that the total number of points."
            print(m) if log is None else log.error(m)
            m = "\n"
            print(m) if log is None else log.info(m)

        else:
            # Save data and individual figures of energy
            individual_energy_figures_function(e, args, log)

            # Grouped energy terms
            if args.groupterm:
                grouped_energy_figures_function(e, args, log)

            # Averages calculation
            if args.avg:
                avg_energy_function(e, args, log)

            # Autocorrelation plots
            if args.acf_list:
                acf_energy_function(e, args, end_time, log)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

# =============================================================================
if __name__ == "__main__":

    main_app()
