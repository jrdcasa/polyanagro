import argparse
import glob
import os
import numpy as np
from scipy.optimize import least_squares
import matplotlib
import matplotlib.pyplot as plt
import utils
import sys
import datetime
import pandas as pd
from subprocess import PIPE, Popen
from polyanagro.Energy import energy_analysis


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

    # ****** MAIN ARG PARSER ******
    desc = """Tg fitting.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)
    subparser = parser.add_subparsers(dest='command', required=True)

    annealing = subparser.add_parser('annealing')
    stepwise = subparser.add_parser('stepwise')
    refit = subparser.add_parser('refit')

    # Annealing
    annealing.add_argument("-e", "--energy", dest="energy_list", nargs="+",
                         help="Energy file from MD package.\n "
                             "The package is detected by the extension of the file",
                         action="store", required=True, default=None)

    annealing.add_argument("--bin", dest="bin_temp", type=float,
                        help="Bin width for the temperature.\n ",
                        action="store", required=False, default=1)

    annealing.add_argument("--log", dest="log",
                         help="Name of the file to write logs from this command",
                         action="store", required=False, default="tg_analysis_anneal_info.log")

    annealing.add_argument("--range", dest="range", nargs=2,
                          help="Range to refit the data",
                          action="store", required=True)

    annealing.add_argument("--tgguess", dest="tgguess", type=float,
                          help="Guess value for Tg in Kelvin",
                          action="store", required=True)

    annealing.add_argument("--scaleparam", dest="scaleparam", type=float,
                          help="This rescales the temperatures so that the characteristic domain "
                               "of the hyperbola is in a box with unit sides. This is needed to ensure that "
                               "the fit is agnostic to the actual temperature range, thereby facilitating fitting",
                          action="store", required=True)

    # Stepwise
    stepwise.add_argument("-d", "--direnergy", dest="energy_list", nargs="+",
                         help="Energy file from MD package.\n "
                             "The package is detected by the extension of the file",
                         action="store", required=True, default=None)

    stepwise.add_argument("--log", dest="log",
                         help="Name of the file to write logs from this command",
                         action="store", required=False, default="tg_analysis_stepwise_info.log")

    stepwise.add_argument("--bin", dest="bin_temp", type=float,
                          help="Bin width for the temperature.\n ",
                          action="store", required=False, default=1)
    stepwise.add_argument("--range", dest="range", nargs=2,
                          help="Range to refit the data",
                          action="store", required=True)

    stepwise.add_argument("--tgguess", dest="tgguess", type=float,
                          help="Guess value for Tg in Kelvin",
                          action="store", required=True)
    stepwise.add_argument("--scaleparam", dest="scaleparam", type=float,
                          help="This rescales the temperatures so that the characteristic domain "
                               "of the hyperbola is in a box with unit sides. This is needed to ensure that "
                               "the fit is agnostic to the actual temperature range, thereby facilitating fitting",
                          action="store", required=True)

    # Refit
    refit.add_argument("--data", dest="data_refit",
                       help="Density versus Tg to refit",
                       action="store", required=True, default=None)
    refit.add_argument("--log", dest="log",
                       help="Name of the file to write logs from this command",
                       action="store", required=False, default="tg_analysis_refit_info.log")
    refit.add_argument("--range", dest="range", nargs=2,
                          help="Range to refit the data",
                          action="store", required=True)
    refit.add_argument("--tgguess", dest="tgguess", type=float,
                          help="Guess value for Tg in Kelvin",
                          action="store", required=True)
    refit.add_argument("--scaleparam", dest="scaleparam", type=float,
                          help="This rescales the temperatures so that the characteristic domain "
                               "of the hyperbola is in a box with unit sides. This is needed to ensure that "
                               "the fit is agnostic to the actual temperature range, thereby facilitating fitting",
                          action="store", required=True)


    args = parser.parse_args()

    ext_list = []
    list_edrs = []
    if args.command != "refit":
        # Check the extension of the energy files
        list_edrs = []
        for iener in args.energy_list:
            # Is a file?
            if os.path.isfile(iener):
                ext_list.append(os.path.splitext(iener)[1])
                list_edrs.append(iener)
            else:  # Is a dir?
                if os.path.isdir(iener):
                    ext_list.append(".edr")
                    path = os.path.join(iener, "*.edr")
                    list_edrs += sorted(glob.glob(path), reverse=True)
                else:
                    print("\nERROR: File or directory {} does not exist\n".format(iener))
                    exit()
        if ext_list[0] == ".edr":
            mdpackage = "GROMACS"
        elif ext_list[0] == ".log":
            mdpackage = "LAMMPS"
        else:
            mdpackage = None
            m = "\nERROR: Extension '{}' for energy files is unkwown.\n".format(ext_list[0])
            print(m)
            exit()

    if args.command == "refit":
        if not os.path.isfile(args.data_refit):
            print("\nERROR: File or directory {} does not exist\n".format(args.data_refit))
            exit()

    if float(args.range[1]) < float(args.range[0]):
            print("\nERROR: First value in the range must be lower than the second 1st: {} 2nd: {}\n".
                  format(args.range[0], args.range[1]))
            exit()

    return args, list_edrs

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                    Glass Transtion Temperature Analysis 
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
    m += "\t\t\t{}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)

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

    m = "\t\t**** TERMS DENSITY and TEMPERATURE PRESENT IN THE ENERGY FILE ****\n"
    try:
        e._df["Density"]
        m +="\t\tDensity\n"
    except KeyError:
        m1 = "ERROR!!!! Density is not present in the energy file. {}".format(args.energy)
        print(m1) if log is None else log.info(m1)
        exit()
    try:
        e._df["Temperature"]
        m += "\t\tTemperature\n"
    except KeyError:
        m1 = "ERROR!!!! Temperature is not present in the energy file. {}".format(args.energy)
        print(m1) if log is None else log.info(m1)
        exit()
    m += "\t\t**** END TERMS DENSITY and TEMPERATURE PRESENT IN THE ENERGY FILE ****\n"
    print(m) if log is None else log.info(m)

    return start_time, end_time

# =============================================================================
def extract_data_function(e, temp_group_value=1.0):

    lines = '#Time(ps) Temperature(K) Density(g/cm^3)\n'
    select_columns = e._df[['Time', 'Temperature', 'Density']]

    # Write the selected columns to a file with specific formatting
    with open('density_temperature_raw.csv', 'w') as fraw:
        # Write data rows
        for index, row in select_columns.iterrows():
            lines += '{0:8.1f},{1:6.1f},{2:8.4f}\n'.format(row['Time'], row['Temperature'], row['Density']/1000)
        fraw.writelines(lines)

    df = pd.DataFrame(e._df, columns=['Time', 'Temperature', 'Density'])

    df.head()
    # Round to the nearest interger
    df['Temperature'] = df['Temperature'].round(0)
    # Create bins for every temp_group_value degrees
    df['Temperature_bin'] = (df['Temperature'] // temp_group_value) * temp_group_value
    df['Temperature_bin'] +=  temp_group_value/2.0

    # Group
    group = df.groupby('Temperature_bin')['Density'].agg(['mean', 'std'])
    # Replace NaN values with 0
    group.fillna(0, inplace=True)

    # Write the selected columns to a file with specific formatting
    lines = "#Index Temp(K) Density(g/cm3) Std(g/cm3)\n"
    with open('density_temperature_{0:03d}K.csv'.format(int(temp_group_value)), 'w') as fraw:
        # Write data rows
        idx = 0
        for index, row in group.iterrows():
            lines += '{0:04d}, {1:8.1f},{2:8.4f},{3:8.4f}\n'.format(idx, index, row['mean']/1000, row['std']/1000)
            idx += 1
        fraw.writelines(lines)

    return group

# =============================================================================
def hyperbolafun_2(hyper_param, x):

    """
    Evaluates a parameterized hyperbola
    """

    x0 = hyper_param[0]
    y0 = hyper_param[1]
    alpha = hyper_param[2]
    gamma = hyper_param[3]
    delta = hyper_param[4]

    dx = x - x0
    hx = y0 - alpha * dx - gamma * (dx / 2 + np.sqrt(dx**2 / 4 + np.exp(-delta)))
    return hx

# =============================================================================
def hyperbola_angle_2(hyper_param, r):

    """
    Returns x at which the hyperbola slope is the fraction r and (1-r) of its two asymptotic values.
    The parameter 0<r<1 can be thought of as the percent convergence of the hyperbola to its asymptotes.
    Here we assume the hyperbola parametrization according to the function hyperbolafun_2.m
    hyper_param = [x0;y0;alpha;gamma;delta]
    """

    x0 = hyper_param[0]
    delta = hyper_param[4]

    dX = np.exp(-delta / 2) * (2 * r - 1) / np.sqrt(r * (1 - r))
    xLow = x0 - dX
    xHigh = x0 + dX

    xBnds = np.array([xLow, xHigh])
    return xBnds

# =============================================================================
def residuals(params, x, y, scalepower):
    if scalepower == 0:
        return hyperbolafun_2(params, x) - y
    else:
        return (hyperbolafun_2(params, x) / (x ** scalepower)) - (y / (x ** scalepower))

# =============================================================================
def tg_bootstrap(densities, temperatures, num_synthetic_sets, tgguess_0=300, scaleparam_0=300, rejection=0, showfigs=0, logger=None):

    """
    Performs a bootstrap analysis to evaluate the uncertainty in the glass transition temperature Tg.

    Parameters:
    densities (np.array): Vector of densities
    temperatures (np.array): Vector of temperatures corresponding to the densities
    numsyntheticsets (int): Number of synthetic datasets to generate
    rejection (int): Flag for running the rejection criterion (default is 0, disabled)
    showfigs (int): Flag for showing figures (default is 0, disabled)
    runinparallel (int): Flag for running synthetic dataset computations in parallel (default is 0, disabled)

    Returns:
    tuple: (meantg, withinvar, trialcentertgs, rejected)
    """

    import warnings
    warnings.filterwarnings("ignore")

    m = "\t\t**** GLASS TRANSITION FITTING INFO ****"
    print(m) if logger is None else logger.info(m)

    """
    % *** NOTES TO USER; these parameters are needed for stability of the
    % fitting process
    % Set tgguess to a rough estimate of where you think Tg is
    % I also find that scaleparam should be set to a value that is roughly
    % where Tg is (or perhaps slightly below).
    """
    # Guess value for Tg in Kelvin
    tgguess = tgguess_0
    # This rescales the temperatures so that the characteristic domain of
    # the hyperbola is in a box with unit sides.
    # This is needed to ensure that the fit is agnostic to the actual temperature range, thereby
    # facilitating fitting
    scaleparam = scaleparam_0
    # If greater than zero, terms in the least squares fit will be scaled by (i.e. divided by) temperature 
    # to this power. Setting this to anything besides zero performs a weighted least squares
    scalepower = 1.25
    # What should the temperature scaling of the noise model be?  Note, this parameter is different from
    # scalepower in that this is the scaling used to
    # generate synthetic noise, not weighted least
    # squares.  Note that residuals will be divided by temperature raised to the power tempscale
    tempscale = 1.25
    # Greater than zero, less than one
    # how converged do you require the data
    # to be to the hyperbola asympototes?
    Percent_converged = 0.90
    # Debugging
    # Allows you to dial down or up the noise residuals; 1 is default
    residscalefactor = 1
    # Set to 1 to see individual hyperbola fits for synthetic data
    showraw2 = 0

    # Re-organizing data and initializing a few variables
    numtemps = len(temperatures)
    testtemps = np.linspace(min(temperatures), max(temperatures), num=100)

    # Guess parameters for the(rescaled) hyperbola fit
    xo = [tgguess / scaleparam, 1, 1, 1, 1]

    # 
    result = least_squares(residuals, xo, args=(temperatures / scaleparam, densities, scalepower))
    popt = result.x

    # Hyperbola fit parameters
    savedouts = popt
    # Vector of residuals and the sum of errors squared

    rvec = densities - hyperbolafun_2(popt, temperatures / scaleparam)
    resnormsx = np.sum(rvec ** 2)

    m = "\t\tStarting analysis for current dataset.\n"
    m += "\t\tTg from initial fit (in K): {0:4.1f}\n".format(popt[0] * scaleparam)
    m += "\t\tNormed residual (smaller is better): {0:10.6f}".format(resnormsx)
    print(m) if logger is None else logger.info(m)

    # Rejection criterion
    # The next few lines will use the hyperbola fits to determine where the
    # asymptotic regimes are relative to the data.
    rejected = 0
    if rejection == 1:
        m = "\t\tConvergence parameter set to: {0:3.1f}%".format(100 * Percent_converged)
        print(m) if logger is None else logger.info(m)
        bnds = hyperbola_angle_2(popt, Percent_converged)
        lowerbound = bnds[0] * scaleparam
        upperbound = bnds[1] * scaleparam

        if showfigs == 1:
            plt.figure(1)
            tempouts = np.linspace(min([0.9 * lowerbound, 1.1 * lowerbound, min(temperatures)]), max([1.1 * upperbound, max(temperatures)]), 100)
            plt.plot(temperatures, densities, 'bx')
            plt.plot(tempouts, hyperbolafun_2(popt, tempouts / scaleparam), 'r')
            plt.xlabel('Temperatures (in K)')
            plt.ylabel('Density (g/cm^3)')
            plt.title('Best fit hyperbola using all simulated data')
            plt.plot([lowerbound, upperbound], [1, 1], 'g', linewidth=3)
            plt.plot(popt[0] * scaleparam, popt[1], 'bo')
            plt.plot([80, popt[0] * scaleparam], [popt[1] - popt[2] * (80 / scaleparam - popt[0]), popt[1]], 'k')
            plt.legend(['Simulated data', 'Best-fit hyperbola', 'Non-asymptotic regimes', 'Hyperbola Center', 'Asymptotes'])
            plt.plot([popt[0] * scaleparam, 500], [popt[1], popt[1] - (popt[2] + popt[3]) * (500 / scaleparam - popt[0])], 'k')
            plt.savefig('figure1.png')
            # plt.show()
        create_gnu_template(temperatures, densities, popt, tempouts, scaleparam, lowerbound, upperbound)

        if lowerbound < min(temperatures) or upperbound > max(temperatures):
            rejected = 1
            m = "\t\tData does not satisfy requirements for extracting Tg"
            print(m) if logger is None else logger.info(m)

        else:
            m = "\t\tDataset not rejected"
            print(m) if logger is None else logger.info(m)


    resids = (densities - hyperbolafun_2(popt, temperatures / scaleparam)) * residscalefactor

    # Noise modeling
    trialsigma = np.std(resids / (temperatures ** tempscale))
    trialnoise = np.random.normal(0, 1, (num_synthetic_sets, numtemps)) *\
                 trialsigma * (temperatures ** tempscale).reshape(1, -1)
    synthetic_densities = hyperbolafun_2(popt, temperatures / scaleparam) + trialnoise[0, :]

    if showfigs == 1:
        plt.figure(2)
        plt.plot(temperatures, densities, 'bx')
        plt.plot(testtemps, hyperbolafun_2(savedouts, testtemps / scaleparam), 'r')
        plt.xlabel('Temperatures (in K)')
        plt.ylabel('Density (g/cm^3)')
        plt.title('Best fit hyperbola using all simulated data')
        plt.savefig('figure2.png')
        # plt.show()

        plt.figure(3)
        plt.plot(temperatures, resids, 'x')
        plt.xlabel('Temperatures (in K)')
        plt.ylabel('Density residuals (g/cm^3)')
        plt.title('Residuals from best fit hyperbola')
        plt.savefig('figure3.png')
        # plt.show()

        plt.figure(4)
        plt.plot(temperatures, resids / (temperatures ** tempscale), 'x')
        plt.xlabel('Temperatures (in K)')
        plt.ylabel('Scaled residuals in units of Density/Temperature^{3/2} [g/(K^{3/2} \times cm^3)]')
        plt.title('Scaled residuals')
        plt.savefig('figure4.png')
        # plt.show()

        plt.figure(5)
        plt.plot(temperatures, trialnoise[0, :], 'x')
        plt.xlabel('Temperatures (in K)')
        plt.ylabel('Synthetic residuals based on noise model (g/cm^3)')
        plt.title('Example of synthetic residuals of best fit hyperbola')
        plt.savefig('figure5.png')
        # plt.show()

        plt.figure(6)
        plt.plot(temperatures, synthetic_densities, 'bx')
        plt.xlabel('Temperatures (in K)')
        plt.ylabel('Synthetic Density (g/cm^3)')
        plt.title('Example of Synthetic Density vs Temperature curve')
        plt.plot(testtemps, hyperbolafun_2(popt, testtemps / scaleparam), 'r')
        # plt.show()
        plt.savefig('figure6.png')

    trial_center_tgs = np.zeros(num_synthetic_sets)

    for j in range(num_synthetic_sets):
        synthetic_densities = hyperbolafun_2(savedouts, temperatures / scaleparam) + trialnoise[j, :]
        result = least_squares(residuals, xo, args=(temperatures / scaleparam, synthetic_densities, scalepower))
        popt = result.x
        trial_center_tgs[j] = popt[0] * scaleparam

    if showfigs == 1:
        plt.figure(7)
        plt.hist(trial_center_tgs)
        plt.xlabel('T_g (in K)')
        plt.ylabel('Number of synthetic sets')
        plt.title('T_g values from hyperbola center (entirely synthetic data)')
        plt.savefig('figure7.png')
        #plt.show()

    x = savedouts[0] * scaleparam
    y = trial_center_tgs
    ste = np.std(trial_center_tgs)

    m = "\t\tEstimated Tg (in K): {0:5.2f}\n".format(x)
    m += "\t\tStandard Error: {0:5.2f}\n".format(ste)
    m += "\t\tRejected: {}".format(rejected)
    print(m) if logger is None else logger.info(m)

    m = "\t\t**** END GLASS TRANSITION FITTING INFO ****"
    print(m) if logger is None else logger.info(m)

    return x, y, ste, rejected

# =============================================================================
def uncertainty_quantification_hyperbola(df, tgguess_0=300, scaleparam_0=300, logger=None):

    # How many synthetic sets to generate
    numsyntheticsets=240
    # Options; set to 1 to enable, 0 to disable
    testreject=1           # perform asymptotic convergence test to see if dataset would be rejected
    showfigs=1             # show figures
    runinparallel=0        # generate and analyze synthetic sets in parallel; requires parallel toolbox

    density_mean = df['mean'].to_numpy()/1000 # g/cm^3
    temperature = df.index.to_numpy()  # K

    tg_bootstrap(density_mean, temperature, numsyntheticsets, tgguess_0=tgguess_0, scaleparam_0=scaleparam_0,
                 rejection=testreject,
                 showfigs=showfigs, logger=logger)

    return None

# =============================================================================
def create_gnu_template(temperatures, densities, popt, tempouts, scaleparam, lowerbound, upperbound):

    # Write data to produce gnuplots
    with open("hyperbola_fit_a.dat", 'w') as f1:
        line = "#Temperature(K) Density(g/cm3)\n"
        for i in range(len(temperatures)):
            line += "{0:4.1f} {1:6.4f}\n".format(temperatures[i], densities[i])
        f1.writelines(line)

    # Write data to produce gnuplots
    with open("hyperbola_fit_b.dat", 'w') as f1:
        line = "#Tempouts(K) Hyperbolafun(g/cm3)\n"
        set2 = hyperbolafun_2(popt, tempouts / scaleparam)
        for i in range(len(tempouts)):
            line += "{0:4.1f} {1:6.4f}\n".format(tempouts[i], set2[i])
        f1.writelines(line)

    # Write data to produce gnuplots
    with open("hyperbola_fit_c.dat", 'w') as f1:
        line = "#X1 Y1 X2 Y2 X3 Y3\n"
        set1 = [80, popt[0] * scaleparam]
        set2 = [popt[1] - popt[2] * (80 / scaleparam - popt[0]), popt[1]]
        set3 = [popt[0] * scaleparam, 500]
        set4 = [popt[1], popt[1] - (popt[2] + popt[3]) * (500 / scaleparam - popt[0])]
        set5 = [lowerbound, upperbound]
        set6 = [1, 1]
        for i in range(2):
            line += "{0:4.1f} {1:6.4f} {2:4.1f} {3:6.4f} {4:f} {5:f}\n".format(set1[i], set2[i],
                                                                               set3[i], set4[i],
                                                                               set5[i], set6[i])
        f1.writelines(line)


    lines = 'reset\n'
    lines += '# ====================================================================\n'
    lines += 'set style line 1 lt 1 ps 1.2 lc rgb "black"   pt 6  lw 1.0\n'
    lines += 'set style line 2 lt 1 ps 1.2 lc rgb "black"   pt 6  lw 1.5\n'
    lines += 'set style line 3 lt 1 ps 1.2 lc rgb "blue"    pt 6  lw 1.5\n'
    lines += 'set style line 4 lt 1 ps 1.2 lc rgb "green"   pt 6  lw 2.5\n'
    lines += '# ====================================================================\n'
    lines += 'set term qt 1 enhanced dashed size 600,400 font "{/Arial:Bold}"\n'
    lines += 'set encoding iso_8859_1\n'
    lines += 'set multiplot layout 1,1\n'
    lines += '\n'
    lines += '# ============ PLOT 1 ==================\n'
    lines += 'set xlabel "Temperature (K)" font "Arial,16"\n'
    lines += '#set ylabel "{/:Bold R_g^2 ({\305}^2)}" font "Arial,16" offset -2\n'
    lines += 'set ylabel "Density (g/cm^3)" font "Arial,16" offset -2\n'
    lines += 'set format x "%4.0f"\n'
    lines += 'set format y "%4.3f"\n'
    lines += 'set xrange[50:750]\n'
    lines += 'set xtics 100 format "{/:Bold {/=12 %.0f}}"\n'
    lines += 'set ytics 0.05 format "{/:Bold {/=12 %.3f}}"\n'
    lines += 'set mxtics 2\n'
    lines += 'set mytics 2\n'
    lines += 'set lmargin 12\n'
    lines += 'set key font ",8"\n'
    lines += 'unset grid\n'
    lines += '\n'
    lines += 'set title ""\n'
    lines += 'set label 1 "" at {}, {} point pointtype 7 pointsize 1.4 \n'.format(popt[0] * scaleparam, popt[1])
    lines += 'p "./hyperbola_fit_a.dat" u 1:2 ls 1 title "Simulated data", \\\n'
    lines += '  "./hyperbola_fit_b.dat" u 1:2 ls 2 w l title "Best-fit hyperbola",\\\n'
    lines += '  "./hyperbola_fit_c.dat" u 1:2 ls 3 w l title "Asymptotes",\\\n'
    lines += '  "./hyperbola_fit_c.dat" u 3:4 ls 3 w l notitle ,\\\n'
    lines += '  "./hyperbola_fit_c.dat" u 5:6 ls 4 w l title "Non-asymptotic regimes"\n'
    lines += '  \n'
    lines += '\n'
    lines += 'unset multiplot\n'

    with open("hyperbola_fit.gnu", 'w') as fgnu:
        fgnu.writelines(lines)

# =============================================================================
def main_app():

    import polyanagro as pag

    # Parse arguments
    args, list_edrs = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)

    # # All data in edr file is represented in a png file. The files
    # # are stored in ./test01_a_figures
    if args.command == "refit":
        list_edrs = [args.data_refit]
    e = energy_analysis(list_edrs, logger=log)
    e.read_energy()
    start_time, end_time = energy_info_function(e, args, log)
    # Extract data
    if args.command == "refit":
        d_vs_T_df = e._df
    else:
        d_vs_T_df = extract_data_function(e, temp_group_value=args.bin_temp)
    # Perform fitting analysis
    # http://dx.doi.org/10.1016/j.polymer.2016.01.074
    min_value = float(args.range[0])
    max_value = float(args.range[1])
    filtered_df = d_vs_T_df[(d_vs_T_df.index >= min_value) & (d_vs_T_df.index <= max_value)]
    uncertainty_quantification_hyperbola(filtered_df, tgguess_0=args.tgguess,
                                         scaleparam_0=args.scaleparam, logger=log)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

# =============================================================================
if __name__ == "__main__":

    main_app()
