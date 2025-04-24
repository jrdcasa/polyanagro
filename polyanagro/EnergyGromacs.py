import panedr
from polyanagro.Energy import Energy
from polyanagro.Custom_Plots import Custom_Plots
import pandas as pd

# =============================================================================
class EnergyGromacs(Energy):

    # =============================================================================
    def __init__(self, energy_list_filenames, logger=None):

        Energy.__init__(self, energy_list_filenames=energy_list_filenames, logger=logger)

        self._energylabels = ["Bond", "Angle", "Ryckaert-Bell.", "LJ-14", "Coulomb-14", "LJ (SR)",
                              "Disper. corr.", "Coulomb (SR)", "Coul. recip.", "Potential", "Kinetic En.",
                              "Total Energy", "Conserved En.", "Proper-Dih.", "Improper-Dih.", "Improper Dih.",
                              "pV", "Enthalpy"]
        self._timelabels = ["Time"]
        self._temperaturelabels = ["Temperature", "T-System"]
        self._pressurelabels = ["Pres. DC", "Pressure", "#Surf*SurfTen"]
        self._virialcomplabels = ["Vir-XX", "Vir-XY", "Vir-XZ",
                                  "Vir-YX", "Vir-YY", "Vir-YZ",
                                  "Vir-ZX", "Vir-ZY", "Vir-ZZ"]
        self._pressurecomplabels = ["Pres-XX", "Pres-XY", "Pres-XZ",
                                    "Pres-YX", "Pres-YY", "Pres-YZ",
                                    "Pres-ZX", "Pres-ZY", "Pres-ZZ"]
        self._dimensionlabels = ['Box-X', 'Box-Y', 'Box-Z']
        self._volumelabels = ['Volume']
        self._densitylabels = ['Density']
        self._nounits = ['Constr. rmsd', "Lamb-System"]

    # =========================================================================
    def read_energy(self, skip_lines=10):

        """
        Class to handle energy files from GROMACS.

        """

        m = "\t\tReading energy data from: \n"
        for iedr in self._energy_list_filenames:
            m += "\t\t\t{}\n".format(iedr)
        print(m) if self._logger is None else self._logger.info(m)

        # Read the EDR file
        # Annealing
        if len(self._energy_list_filenames) == 1:
            self._df = panedr.edr_to_df(iedr)
            self._properties = self._df.columns.tolist()
        # Stepwise
        else:
            mean_values_list = []
            std_values_list = []
            for iedr in self._energy_list_filenames:
                tmp_df = panedr.edr_to_df(iedr)
                # Discard the first 10 points to calculate the average
                tmp_df_discart_first_lines = tmp_df.iloc[skip_lines:]
                mean_values_list.append(pd.DataFrame(tmp_df_discart_first_lines.mean()).transpose())
                std_values_list.append(pd.DataFrame(tmp_df_discart_first_lines.std()).transpose())
            self._df = pd.concat(mean_values_list, ignore_index=True)
            self._df_std = pd.concat(std_values_list, ignore_index=True)
            self._properties = self._df.columns.tolist()
            self._properties_std = self._df_std.columns.tolist()

    # =========================================================================
    def plot_energy_group(self, skip_data=0, path_to_save="."):

        cplt = Custom_Plots()

        energy_group = ["Bonded", "Non-Bonded", "1-4 terms",
                        "Dispersion correction", "Potential", "Total"]
        bonded = ["Bond", "Angle", "Ryckaert-Bell.", "Proper-Dih.", "Improper-Dih."]
        nonbonded = ["LJ (SR)", "Coulomb (SR)", "Coul. recip."]
        nonbonded14 = ["LJ-14", "Coulomb-14"]
        dispersion_corr = ["Disper. corr."]

        # Group energy plots
        nsubplots = [2, 3]
        nplots = 2*3

        xlabel = self._properties[0]
        xdata = self._df[xlabel]
        gromacs_df = pd.DataFrame(xdata)

        # Bonded energy Group
        sum_bonded = pd.DataFrame(0, index=self._df['Time'], columns=['Bonded Energy'])
        for ibond in bonded:
            try:
                if len(sum_bonded) == 0:
                    sum_bonded['Bonded Energy'] = self._df[ibond]
                else:
                    sum_bonded['Bonded Energy'] += self._df[ibond]
            except KeyError:
                continue
        # NonBonded energy Group
        sum_nonbonded = pd.DataFrame(0, index=self._df['Time'], columns=['Nonbonded Energy'])
        for inb in nonbonded:
            try:
                if len(sum_nonbonded) == 0:
                    sum_nonbonded['Nonbonded Energy'] = self._df[inb]
                else:
                    sum_nonbonded['Nonbonded Energy'] += self._df[inb]
            except KeyError:
                continue

        # 1-4 Energy Group
        sum_14 = pd.DataFrame(0, index=self._df['Time'], columns=['Nonbonded14 Energy'])
        for inb in nonbonded14:
            try:
                if len(sum_14) == 0:
                    sum_14['Nonbonded14 Energy'] = self._df[inb]
                else:
                    sum_14['Nonbonded14 Energy'] += self._df[inb]
            except KeyError:
                continue

        # Dispersion correction
        disp_corr = pd.DataFrame(0, index=self._df['Time'], columns=['Dispersioncorr'])
        for inb in dispersion_corr:
            try:
                if len(disp_corr) == 0:
                    disp_corr['Dispersioncorr'] = self._df[inb]
                else:
                    disp_corr['Dispersioncorr'] += self._df[inb]
            except KeyError:
                continue

        # Join dataframes
        if len(xdata) != len(sum_bonded):
            m = "Length of time ({}) and bonded ({}) series are different.".format(len(xdata), len(sum_bonded))
            print(m) if self._logger is None else self._logger.warning(m)
            exit()
        elif len(xdata) != len(sum_nonbonded):
            m = "Length of time ({}) and nonbonded ({}) series are different.".format(len(xdata), len(sum_nonbonded))
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif len(xdata) != len(sum_14):
            m = "Length of time ({}) and nonbonded-14 ({}) series are different.".format(len(xdata), len(sum_14))
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif len(xdata) != len(disp_corr):
            m = "Length of time ({}) and dispersion correction ({}) series are different.".\
                format(len(xdata), len(disp_corr))
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            gromacs_df = gromacs_df.assign(bonded=sum_bonded.values,
                                           nonbonded=sum_nonbonded.values,
                                           nonbonded14=sum_14.values,
                                           dispersioncorr=disp_corr.values,
                                           potential=self._df["Potential"],
                                           total_energy=self._df["Total Energy"])

            gromacs_df.columns=["Time", "Bond energy", "Non-Bonded energy", "1-4 nonbonded energy",
                                "Dispersion correction", "Potential energy", "Total energy"]

        yunits_list = ["kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol"]
        cplt.group_line_plot(nsubplots, gromacs_df, begin=skip_data, path_to_save=path_to_save,
                             xunits="ps", yunits=yunits_list, fig_size=(15, 10))

    # =========================================================================
    def plot_density_group(self, skip_data=0, path_to_save="."):

        cplt = Custom_Plots()

        # Group density plots
        nsubplots = [2, 4]

        xlabel = self._properties[0]
        xdata = self._df[xlabel]
        gromacs_df = pd.DataFrame(xdata)

        gromacs_df = gromacs_df.assign(density=self._df["Density"],
                                       volume=self._df["Volume"],
                                       pressure=self._df["Pressure"],
                                       temperature=self._df["Temperature"],
                                       kinetic=self._df["Kinetic En."],
                                       box_x=self._df["Box-X"],
                                       box_y=self._df["Box-Y"],
                                       box_z=self._df["Box-Z"])

        gromacs_df.columns=["Time", "Density", "Volume", "Pressure", "Temperature",
                            "Kin. Energy", "Box-X", "Box-Y", "Box-Z"]
        yunits_list = ["kg/cm^3", "nm^3", "bar", "K", "kJ/mol", "nm", "nm", "nm"]
        cplt.group_line_plot(nsubplots, gromacs_df, begin=skip_data, path_to_save=path_to_save,
                             xunits="ps", yunits=yunits_list, fig_size=(20, 10), outputname='density.jpg')



