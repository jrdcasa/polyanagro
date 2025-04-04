import panedr
from polyanagro.Energy import Energy
from polyanagro.Custom_Plots import Custom_Plots
from collections import defaultdict
import pandas as pd
import numpy as np
import re



# =============================================================================
class EnergyLammps(Energy):

    # =============================================================================
    def __init__(self, energy_list_filenames, logger=None):

        Energy.__init__(self, energy_list_filenames=energy_list_filenames, logger=logger)

        self._energylabels = ["E_bond", "E_angle", "E_dihed", "E_impro", "E_pair", "E_vdwl",
                              "E_coul", "E_long", "E_tail", "PotEng"]
        self._timelabels = ["Time"]
        self._temperaturelabels = ["Temp"]
        self._pressurelabels = ["Press"]
        self._dimensionlabels = ['Lx', 'Ly', 'Lz']
        self._densitylabels = ['Density']
        self._timesteplabel = ['timestep ']
        self._timestep = None  # in ps

        self._virialcomplabels = []
        self._pressurecomplabels = []
        self._volumelabels = []
        self._nounits = []


    # =========================================================================
    def read_energy(self):

        """
        Class to handle energy files from GROMACS.
        """

        m = "\t\tReading energy data from {}\n".format(self._energy_list_filenames)
        print(m) if self._logger is None else self._logger.info(m)

        # Read the LOG file from LAMMPS
        self._log_to_df()
        self._properties = self._df.columns.tolist()

    # =========================================================================
    def _log_to_df(self):

        """
        Read energy terms in a log LAMMPS file creating a
        DataFrame (Pandas) with the data (self._df).
        """

        kcal_to_kJ = 4.184
        angs_to_nano = 0.1
        atm_to_bar = 1.01325

        # Starting the dictionaty with the default labels
        dict_df = defaultdict()
        dict_idx_label = defaultdict()
        with open(self._energy_list_filenames, 'r') as flog:
            lines = flog.readlines()
            flog.seek(0)
            for num, iline in enumerate(flog, 1):
                # Get the timestep and transform in ps
                if re.match(r"(\s*)timestep (\s*)", iline):
                    self._timestep = float(iline.split()[1])/1000  # in ps
                elif re.match(r"^Step", iline):
                    label_list = iline.split()
                    label_list.insert(0, "Time")
                    idx_pos_startenergy = num
                    j = 0
                    for ilabel in label_list:
                        dict_df[ilabel] = []
                        dict_idx_label[j] = ilabel
                        j += 1

        # Locate the final position of the energies
        for idx in range(idx_pos_startenergy, len(lines)):
            if lines[idx].find("Loop time") != -1:
                idx_pos_endenergy = idx
                break

        # Get energy data
        newindex = []
        for idx in range(idx_pos_startenergy, idx_pos_endenergy):
            for num, item in enumerate(lines[idx].split(),1):
                label = dict_idx_label[num]
                dict_df[label].append(float(item))
                if num == 1:
                    dict_df['Time'].append(float(item)*self._timestep)
                    newindex.append(float(item)*self._timestep)

        # Change units
        for label, values in dict_df.items():
            if label in self._energylabels:
                m = [element * kcal_to_kJ for element in values]
                dict_df[label] = m         # kJ/mol
            elif label in self._dimensionlabels:
                m = [element * angs_to_nano for element in values]
                dict_df[label] = m         # nm
            elif label in self._densitylabels:
                m = [element * 1000 for element in values]
                dict_df[label] = m        # kg/m3
            elif label in self._pressurelabels:
                m = [element * atm_to_bar for element in values]
                dict_df[label] = m        # bar
            # else:
            #     print(label)

        # Reindex the dataframe with the values of Time
        s = pd.Series(newindex)
        self._df = pd.DataFrame.from_dict(dict_df)
        self._df = self._df.set_index([s])

    # =========================================================================
    def plot_energy_group(self, skip_data=0, path_to_save="."):

        cplt = Custom_Plots()


        energy_group = ["Bonded", "Non-Bonded",
                        "Dispersion correction", "Potential", "Total"]
        bonded = ["E_bond", "E_angle", "E_dihed", "E_impro"]
        nonbonded = ["E_pair"]
        dispersion_corr = ["E_tail"]

        # Group energy plots
        nsubplots = [2, 2]
        nplots = 2*2

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
        elif len(xdata) != len(disp_corr):
            m = "Length of time ({}) and dispersion correction ({}) series are different.".\
                format(len(xdata), len(disp_corr))
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            gromacs_df = gromacs_df.assign(bonded=sum_bonded.values,
                                           nonbonded=sum_nonbonded.values,
                                           dispersioncorr=disp_corr.values,
                                           potential=self._df["PotEng"],)

            gromacs_df.columns=["Time", "Bond energy", "Non-Bonded energy",
                                "Dispersion correction", "Potential energy"]

        yunits_list = ["kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol", "kJ/mol"]
        cplt.group_line_plot(nsubplots, gromacs_df, begin=skip_data, path_to_save=path_to_save,
                             xunits="ps", yunits=yunits_list, fig_size=(15, 10))

    # =========================================================================
    def plot_density_group(self, skip_data=0, path_to_save="."):

        cplt = Custom_Plots()

        # Group density plots
        nsubplots = [2, 3]

        xlabel = self._properties[0]
        xdata = self._df[xlabel]
        gromacs_df = pd.DataFrame(xdata)

        gromacs_df = gromacs_df.assign(density=self._df["Density"],
                                       pressure=self._df["Press"],
                                       temperature=self._df["Temp"],
                                       box_x=self._df["Lx"],
                                       box_y=self._df["Ly"],
                                       box_z=self._df["Lz"])

        gromacs_df.columns=["Time", "Density", "Pressure", "Temperature",
                            "Box-X", "Box-Y", "Box-Z"]
        yunits_list = ["kg/cm^3", "bar", "K", "nm", "nm", "nm"]
        cplt.group_line_plot(nsubplots, gromacs_df, begin=skip_data, path_to_save=path_to_save,
                             xunits="ps", yunits=yunits_list, fig_size=(20, 10), outputname='density.jpg')



