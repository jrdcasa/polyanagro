import os
from abc import ABCMeta, abstractmethod
from polyanagro.Custom_Plots import Custom_Plots
import pandas as pd
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import acf
from scipy.optimize import curve_fit



# =============================================================================
class Energy(metaclass=ABCMeta):

    # =========================================================================
    def __init__(self, energy_list_filenames, logger=None):

        """
        Class to handle energy files from MD packages.

        Args:
            energy_list_filenames (str): Name of the file to be analyzed.
        """

        self._energy_list_filenames = energy_list_filenames
        self._logger = logger
        self._mdpackage = None
        self._df = None
        self._dfavg = None
        self._properties = None

        # Get extension file
        ext = os.path.splitext(energy_list_filenames[0])[1]

        if ext == ".edr":
            self._mdpackage = "GROMACS"
        elif ext == ".log":
            self._mdpackage = "LAMMPS"
        elif ext == ".dat" or ext == ".csv":
            self._mdpackage = "REFIT"
        else:
            m = "\t\tMD package '{}' is not implemented\n".format(self._mdpackage)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

    # =========================================================================
    @abstractmethod
    def read_energy(self):
        pass

    # =========================================================================
    def assign_units(self, label):

        if label in self._energylabels:
            title = label + " energy"
            yunits = "kJ/mol"
        elif label in self._temperaturelabels:
            title = label
            yunits = "K"
        elif label in self._pressurelabels:
            title = label
            yunits = "bar"
        elif label in self._virialcomplabels:
            title = label
            yunits = "kJ/mol"
        elif label in self._pressurecomplabels:
            title = label
            yunits = "bar"
        elif label in self._dimensionlabels:
            title = label
            yunits = "nm"
        elif label in self._volumelabels:
            title = label
            yunits = "nm^3"
        elif label in self._densitylabels:
            title = label
            yunits = "kg/m^3"
        elif label in self._timelabels:
            title = label
            yunits = "ps"
        elif label in self._nounits:
            title = label
            yunits = ""
        else:
            m = "{} is not found in the energy file".format(label)
            print(m) if self._logger is None else self._logger.warning(m)
            title = label
            yunits = ""

        return title, yunits

    # =========================================================================
    def plot_energy_time(self, path_to_save=".", begin=None, label=None):

        if not os.path.isdir(path_to_save):
            m = "{} does not exist. Please create the directory first".format(path_to_save)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        cplt = Custom_Plots()

        xlabel = self._properties[0]
        xdata = self._df[xlabel]
        if label is None:
            for ylabel in self._properties[1:]:
                ydata = self._df[ylabel]
                title, yunits = self.assign_units(ylabel)
                cplt.simple_line_plot(xdata, ydata, title=title, xlabel=xlabel, ylabel=ylabel,
                                      begin=begin, xunits="ps", yunits=yunits, path_to_save=path_to_save)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        pd.set_option('display.colheader_justify', 'center')
        pd.set_option('display.precision', 6)
        filename_csv = os.path.join(path_to_save, 'energy_data.dat')
        self._df.to_csv(filename_csv, sep=' ', float_format="%.6f")

    # =========================================================================
    def avg_data(self, tbegin, tend):

        # =======================
        def linear_fit(x, a, b):
            return a * x + b
        # =======================
        import warnings
        warnings.simplefilter("ignore", Warning)
        # Collect statistics from data
        data = []
        for iseries in self._df:
            label = iseries
            s = self._df[iseries]
            N = len(s[tbegin:tend])
            mean = s[tbegin:tend].mean()
            std = s[tbegin:tend].std()
            maxv = s[tbegin:tend].max()
            minv = s[tbegin:tend].min()
            l , units = self.assign_units(label)
            # Calculate drift of the time series (similar to GROMACS)
            x = s[tbegin:tend].index
            y = s[tbegin:tend].values
            fit_y_fit_a, fit_y_fit_b = curve_fit(linear_fit, x, y)
            y1 = linear_fit(tbegin, fit_y_fit_a[0], fit_y_fit_a[1])
            y2 = linear_fit(tend, fit_y_fit_a[0], fit_y_fit_a[1])
            drift = y2 - y1
            data.append([label, N, mean, std, minv, maxv, drift, units ])

        self._dfavg = pd.DataFrame(data, columns=['Name', 'N', 'Mean', 'Std', 'Min', 'Max', 'Drift', 'Units'])

    # =========================================================================
    def autocorrelation_data(self, label, tbegin, tend):

        series = self._df[label][tbegin:tend]
        figure = plot_acf(series)
        acf_data = acf(series)

        return acf_data

    # =========================================================================
    @abstractmethod
    def plot_energy_group(self, path_to_save="."):
        pass


# =============================================================================
def energy_analysis(energy_filenames, logger=None) -> Energy:

    from polyanagro.EnergyGromacs import EnergyGromacs
    from polyanagro.EnergyLAMMPS import EnergyLammps
    from polyanagro.EnergyDat import EnergyDat

    if type(energy_filenames) is not list:
        energy_filenames = [energy_filenames]

    # Get extension file examining the first item in the list
    if os.path.isdir(energy_filenames[0]):
        ext = os.path.splitext(energy_filenames[0])[1]
    else:
        ext = os.path.splitext(energy_filenames[0])[1]

    if ext == ".edr":
        e = EnergyGromacs(energy_list_filenames=energy_filenames, logger=logger)
    elif ext == ".log":
        e = EnergyLammps(energy_list_filenames=energy_filenames, logger=logger)
    elif ext == ".dat" or ext == ".csv":
        e = EnergyDat(energy_list_filenames=energy_filenames, logger=logger)
    else:
        e = None
        m = "\t\tExtension '{}' for energy is unkwown.\n".format(ext)
        print(m) if logger is None else logger.error(m)
        exit()

    return e
