import glob
import os
import re
import csv
import datetime
import numpy as np
import pandas as pd
from collections import defaultdict
from subprocess import PIPE, Popen

# =============================================================================
class EnergyCohesive:

    __slots__ = ["_logger", "_filenametpr", "_gmxexepath", "_nkinds_molecules",
                 "_mols", "_atoms_molkind", "_ndxfile", "_filenamextc", "_composition",
                 "_dispersioncorrect", "_vdvcutoff", "_coulombcutoff", "_topology", "_dfPotential_isolated",
                 "_dfPotential_full", "_dfDensity_full", "_fraction_trj_avg"]

    # ==========================================================================
    def __init__(self, gmxexepath, filenametpr, filenamextc, filenametopo,
                 fraction_trj_avg, logger=None):

        # Atributtes
        self._gmxexepath = gmxexepath
        self._logger = logger
        self._filenametpr = filenametpr
        self._filenamextc = ""
        self._topology = filenametopo
        for item in filenamextc:
            self._filenamextc += item + " "
        self._nkinds_molecules = defaultdict()
        self._composition = list()
        self._mols = 0
        self._atoms_molkind = list()
        self._ndxfile = None
        self._dispersioncorrect = True
        self._vdvcutoff = 1.2
        self._coulombcutoff = 1.2
        self._dfPotential_isolated = None
        self._dfPotential_full = None
        self._dfDensity_full = None
        self._fraction_trj_avg = fraction_trj_avg

        # Get info of the system from the tpr
        self._getinfo_fromtpr()

        m = "\t\t Number of different kind of molecules: {}\n".format(self._nkinds_molecules)
        # for ikind in range(self._nkinds_molecules):
        #     m += "\t\t   Molecule kind {}: {} molecules with {} atoms each.\n".format(ikind,
        #                                                                               self._mols[ikind],
        #                                                                               self._atoms_molkind[ikind])

        print(m) if self._logger is None else self._logger.info(m)


    # ==========================================================================
    def _getinfo_fromtpr(self):

        # Get information from the full tpr
        workdir = os.getcwd()
        cmd = "{} dump -s {}".format(self._gmxexepath, self._filenametpr)
        rs = self.call_subprocess(cmd, cwd=workdir)

        m = "\t\t Reading the tpr file ({})\n".format(self._filenametpr)
        print(m) if self._logger is None else self._logger.info(m)

        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            output_dump_tpr = rs[0].decode("utf-8").split("\n")
            for num, line in enumerate(output_dump_tpr):

                if re.findall(r"moltype.*=", line):
                    name = line.split()[3]
                    next_line = output_dump_tpr[num+1]
                    nmols = next_line.split()[2]
                    self._composition.append([name, int(nmols)])
                    self._mols += int(nmols)

                if re.findall(r"moltype.*\(", line):
                    # i = line.split()[-1]
                    # i = i.replace("(","")
                    # i = i.replace("):","")
                    next_line = output_dump_tpr[num+1]
                    nnext_line = output_dump_tpr[num+3]
                    name = next_line.split("=")[-1]
                    natoms = nnext_line.split()[-1]
                    natoms = natoms.replace("(", "")
                    natoms = natoms.replace("):","")
                    self._nkinds_molecules[name] = int(natoms)


    # ==========================================================================
    def create_molindex_ndx(self, ndxfile):

        """
        Create a ndx file for molecules
        Args:
            ndxfile:

        Returns:

        """

        self._ndxfile = ndxfile

        m = "\t\t Creating index file ({})\n".format(ndxfile)
        print(m) if self._logger is None else self._logger.info(m)

        idx_global = 1
        imol_global = 1
        with open(ndxfile, "w") as fndx:
            for name, nmols in self._composition:
                for imol in range(nmols):
                    fndx.writelines("[mol{}]\n".format(imol_global))
                    for iatom in range(self._nkinds_molecules[name]):
                        fndx.writelines("{}\n".format(idx_global))
                        idx_global += 1
                    imol_global += 1

    # ==========================================================================
    def get_isolatedmol_energy(self, rvdw=None, rcoulomb=None, dispersioncorrect=None):

        """
        Get the intramolecular (isolated) molecules.
        Returns:

        """

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t Calculating intramolecular potential energy of the isolated chains. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        workdir = os.getcwd()
        # Parameters to calculate the energy
        if rvdw is not None:
            self._vdvcutoff = rvdw
        if rcoulomb is not None:
            self._coulombcutoff = rcoulomb
        if dispersioncorrect == 1:
            self._dispersioncorrect = True
        else:
            self._dispersioncorrect = False
        self._mdp_rerun()

        # Join trajectories if needed
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t\t Join trajectories to trj_full.xtc. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        cmd_join = "gmx trjcat -f {0:s} -o trj_full.xtc".format(self._filenamextc)
        rs = self.call_subprocess(cmd_join, cwd=workdir)
        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        imol_global = 0
        for name, nmols in self._composition:

            fullpath_sctopo = self._createtopology_for_sc(name)

            for imol in range(nmols):

                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\t Calculating energy for molecule {} of {}. ({})".format(imol_global+1, self._mols, now)
                print(m) if self._logger is None else self._logger.info(m)

                # Get TPR for the chain imol_global
                cmd_template_convert_tpr = "{0:s} convert-tpr -s {1:s} -o tmp_{2:04d}.tpr -n {3:s}". \
                    format(self._gmxexepath, self._filenametpr, imol_global, self._ndxfile)
                cmd01 = "echo {} |".format(imol_global) + cmd_template_convert_tpr
                rs = self.call_subprocess(cmd01, cwd=workdir)
                if len(rs[0]) == 0:
                    m = "\n\n" + rs[1].decode("utf-8")
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

                # Get XTC for the chain imol_global
                cmd_template_trjconv_01 = "{0:s} trjconv -n {1:s} -f {2:s} -o SC_traj_ch_{3:04d}.xtc". \
                    format(self._gmxexepath, self._ndxfile, "trj_full.xtc", imol_global)
                cmd02 = "echo {} |".format(imol_global) + cmd_template_trjconv_01
                rs = self.call_subprocess(cmd02, cwd=workdir)
                if len(rs[0]) == 0:
                    m = "\n\n" + rs[1].decode("utf-8")
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

                # Get GRO for the chain imol_global
                cmd_template_trjconv_02 = "{0:s} trjconv -n {1:s} -f {2:s} -s {3:s} " \
                                          "-skip 1000000000 -o SC_traj_ch_{4:04d}.gro". \
                    format(self._gmxexepath, self._ndxfile, "trj_full.xtc", self._filenametpr, imol_global)
                cmd03 = "echo {} |".format(imol_global) + cmd_template_trjconv_02
                rs = self.call_subprocess(cmd03, cwd=workdir)
                if len(rs[0]) == 0:
                    m = "\n\n" + rs[1].decode("utf-8")
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

                # Run grompp to mdrun rerun
                cmd04 = "gmx grompp -f {0:s} -c SC_traj_ch_{1:04d}.gro -p {2:s} -o SC_traj_ch_{1:04d}.tpr -maxwarn 5".\
                    format("rerun.mdp", imol_global, fullpath_sctopo)
                rs = self.call_subprocess(cmd04, cwd=workdir)
                if len(rs[0]) == 0:
                    m = "\n\n" + rs[1].decode("utf-8")
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

                # Run mdrun rerun
                cmd05 = "gmx mdrun -rerun SC_traj_ch_{0:04d}.xtc -s SC_traj_ch_{0:04d}.tpr " \
                        "-deffnm SC_traj_ch_{0:04d}_out".\
                    format(imol_global)
                rs = self.call_subprocess(cmd05, cwd=workdir)

                # Get potential energy
                cmd06 = "echo \"Potential\" | gmx energy -f SC_traj_ch_{0:04d}_out.edr -o SC_energy_{0:04d}.xvg"\
                    .format(imol_global)
                rs = self.call_subprocess(cmd06, cwd=workdir)

                #self._clean_files(imol_global)

                if imol_global == 0:
                    # Initialize the Dataframe
                    skip_rows = 0
                    cname = "Potential Mol {0:04d} (kJ/mol)".format(imol_global+1)
                    fn = "SC_energy_{0:04d}.xvg".format(imol_global)
                    with open(fn, 'r') as f:
                        for line in f:
                            if line.startswith('#') or line.startswith('@') :
                                skip_rows += 1
                            else:
                                break
                    self._dfPotential_isolated = pd.read_table(fn, skipinitialspace=True, delimiter=" ", skiprows=skip_rows,
                                                               names=["#timestep(ps)", cname])
                else:
                    cname = "Potential Mol {0:04d} (kJ/mol)".format(imol_global+1)
                    fn = "SC_energy_{0:04d}.xvg".format(imol_global)
                    tmpdf = pd.read_table(fn, skipinitialspace=True, delimiter=" ", skiprows=skip_rows,
                                          names=["#timestep(ps)", cname])

                    self._dfPotential_isolated = pd.concat([self._dfPotential_isolated, tmpdf[cname]], axis = 1)

                imol_global += 1

    # ==========================================================================
    def get_full_energy(self, rvdw=None, rcoulomb=None, dispersioncorrect=None):

        """
        Get the full energy of the system.
        Returns:

        """

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t Calculating potential energy of the full system. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        workdir = os.getcwd()
        # Parameters to calculate the energy
        if rvdw is not None:
            self._vdvcutoff = rvdw
        if rcoulomb is not None:
            self._coulombcutoff = rcoulomb
        if dispersioncorrect == 1:
            self._dispersioncorrect = True
        else:
            self._dispersioncorrect = False
        self._mdp_rerun()

        cmd01 = "echo 0 | gmx trjconv -f trj_full.xtc -s {0:s} -skip 1000000000 -o full.gro ".format(self._filenametpr)
        rs = self.call_subprocess(cmd01, cwd=workdir)
        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        cmd02 = "gmx grompp -f rerun.mdp  -c full.gro -p {0:s} -o TMPnewtopo.tpr".format(self._topology)
        rs = self.call_subprocess(cmd02, cwd=workdir)
        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        cmd03 = "gmx mdrun -rerun trj_full.xtc -s TMPnewtopo.tpr -deffnm Full_traj_out"
        rs = self.call_subprocess(cmd03, cwd=workdir)

        # Get potential energy
        cmd04 = "echo \"Potential\" | gmx energy -f Full_traj_out.edr -o SC_energy_full.xvg"
        rs = self.call_subprocess(cmd04, cwd=workdir)
        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Get density
        cmd05 = "echo \"Density\" | gmx energy -f Full_traj_out.edr -o SC_density_full.xvg"
        rs = self.call_subprocess(cmd05, cwd=workdir)
        if len(rs[0]) == 0:
            m = "\n\n" + rs[1].decode("utf-8")
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Initialize the Dataframe
        skip_rows = 0
        cname = "Full_Potential (kJ/mol)"
        fn = "SC_energy_full.xvg"
        with open(fn, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('@'):
                    skip_rows += 1
                else:
                    break
        self._dfPotential_full = pd.read_table(fn, skipinitialspace=True, delimiter=" ", skiprows=skip_rows,
                                               names=["#timestep(ps)", cname])
        fn = "SC_density_full.xvg"
        cname = "Density (kg/cm3)"
        self._dfDensity_full = pd.read_table(fn, skipinitialspace=True, delimiter=" ", skiprows=skip_rows,
                                               names=["#timestep(ps)", cname])

        self._clean_files(None)

    # ==========================================================================
    def calc_energycoh(self, totalmass, nmols):

        self._dfPotential_isolated["Sum isolated (kJ/mol)"] = self._dfPotential_isolated.iloc[:,1:].sum(axis=1)
        self._dfPotential_full["Sum_Isolated_Potential (kJ/mol)"] = self._dfPotential_isolated["Sum isolated (kJ/mol)"]
        self._dfPotential_full["Energy_cohesive (kJ/mol)"] = self._dfPotential_isolated["Sum isolated (kJ/mol)"] - \
                                                             self._dfPotential_full["Full_Potential (kJ/mol)"]
        self._dfPotential_full["Molar_volume (cm3/mol)"] = totalmass / (self._dfDensity_full.iloc[:,1]/1000)


        self._dfPotential_full["Sum_Isolated_Potential (kJ/mol)"] /= nmols
        self._dfPotential_full["Full_Potential (kJ/mol)"] /= nmols
        self._dfPotential_full["Energy_cohesive (kJ/mol)"] = self._dfPotential_full["Sum_Isolated_Potential (kJ/mol)"] - \
                                                             self._dfPotential_full["Full_Potential (kJ/mol)"]
        self._dfPotential_full["Molar_volume (cm3/mol)"] = (totalmass/nmols) / (self._dfDensity_full.iloc[:,1]/1000)


        self._dfPotential_full["CED (J/cm3)"] = \
                  1000*self._dfPotential_full["Energy_cohesive (kJ/mol)"] \
                  /    self._dfPotential_full["Molar_volume (cm3/mol)"]

        self._dfPotential_full["Solubility_parameter (J/cm3)^0.5"] = \
                  np.sqrt(1000*self._dfPotential_full["Energy_cohesive (kJ/mol)"] /
                               self._dfPotential_full["Molar_volume (cm3/mol)"])

        self._statistics()

        format_mapping = {'#timestep(ps)': '{:10.2f}',
                          'Full_Potential (kJ/mol)': '{:10.1f}',
                          'Sum_Isolated_Potential (kJ/mol)': '{:10.1f}',
                          'Energy_cohesive (kJ/mol)': '{:10.1f}',
                          'Molar_volume (cm3/mol)': '{:10.1f}',
                          'CED (J/cm3)': '{:10.1f}',
                          'Solubility_parameter (J/cm3)^0.5': '{:10.2f}'}

        for key, value in format_mapping.items():
             self._dfPotential_full[key] = self._dfPotential_full[key].apply(value.format)

        self._dfPotential_full.to_csv("CED_solubility.dat", sep=' ', index=False,
                                      columns=format_mapping, quoting=None, quotechar=' ')



    # #########################################################################
    def _statistics(self):

        nrows = self._dfPotential_full.shape[0]

        nrows_to_avg = int(nrows * self._fraction_trj_avg)

        m = "\n\t\t Number of rows to average: from {0:d} to {1:d}\n".format(nrows_to_avg, nrows)
        tstart = self._dfPotential_full.iloc[nrows_to_avg, 0]
        tend = self._dfPotential_full.iloc[nrows-1, 0]
        m += "\t\t From {0:.1f} to {1:.1f} ps".format(tstart, tend)
        print(m) if self._logger is None else self._logger.info(m)

        if nrows_to_avg >= 0:
            Pot_avg = self._dfPotential_full.iloc[nrows_to_avg:,1].mean()
            Pot_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 1].std()
            IsoPot_avg = self._dfPotential_full.iloc[nrows_to_avg:,2].mean()
            IsoPot_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 2].std()
            Ecoh_avg = self._dfPotential_full.iloc[nrows_to_avg:,3].mean()
            Ecoh_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 3].std()
            Molvol_avg = self._dfPotential_full.iloc[nrows_to_avg:,4].mean()
            Molvol_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 4].std()
            CED_avg = self._dfPotential_full.iloc[nrows_to_avg:,5].mean()
            CED_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 5].std()
            Solpar_avg = self._dfPotential_full.iloc[nrows_to_avg:,6].mean()
            Solpar_avgstd = self._dfPotential_full.iloc[nrows_to_avg:, 6].std()
            Density_avg = (self._dfDensity_full.iloc[:,1]/1000).mean()
            Density_avgstd = (self._dfDensity_full.iloc[:,1]/1000).std()
        else:
            Pot_avg = 999999.9
            Pot_avgstd = 0.0
            IsoPot_avg = 999999.9
            IsoPot_avgstd = 0.0
            Ecoh_avg = 999999.9
            Ecoh_avgstd = 0.0
            Molvol_avg = 999999.9
            Molvol_avgstd = 0.0
            CED_avg = 999999.9
            CED_avgstd = 0.0
            Solpar_avg = 999999.9
            Solpar_avgstd = 0.0
            Density_avg = 999999.9
            Density_avgstd = 0.0

        m = "\n\t*** Average values per chain***\n"
        m += "\t\t Full Potential Energy      (kJ/mol)     = {0:10.2f} +- {1:10.2f} (kJ/mol)\n".\
            format(Pot_avg, Pot_avgstd)
        m += "\t\t Isolated Potential Energy  (kJ/mol)     = {0:10.2f} +- {1:10.2f} (kJ/mol)\n".\
            format(IsoPot_avg, IsoPot_avgstd)
        m += "\t\t Energy cohesive            (kJ/mol)     = {0:10.2f} +- {1:10.2f} (kJ/mol)\n".\
            format(Ecoh_avg, Ecoh_avgstd)
        m += "\t\t Molar volume               (cm^3/mol)   = {0:10.2f} +- {1:10.2f} (cm^3/mol)\n".\
            format(Molvol_avg, Molvol_avgstd)
        m += "\t\t CED                        (J/cm^3)     = {0:10.2f} +- {1:10.2f} (J/cm^3)\n".\
            format(CED_avg, CED_avgstd)
        m += "\t\t Solubility parameter       (J/cm^3)^0.5 = {0:10.2f} +- {1:10.2f} (J/cm^3)^0.5\n".\
            format(Solpar_avg, Solpar_avgstd)
        m += "\t\t\t 1 (J/cm^3)^0.5 = 1 (MPa)^0.5\n"
        m += "\t\t Density                    (g/cm^3)     = {0:10.3f} +- {1:10.3f} (g/cm^3)\n".\
            format(Density_avg, Density_avgstd)
        print(m) if self._logger is None else self._logger.info(m)

        # Template plots dimensions
        diroutput = "./"
        dict_avg = defaultdict()
        dict_avg["Epot"] = [Pot_avg, Pot_avgstd]
        dict_avg["Ecoh"] = [Ecoh_avg, Ecoh_avgstd]
        dict_avg["Molvol"] = [Molvol_avg, Molvol_avgstd]
        dict_avg["CED"] = [CED_avg, CED_avgstd]
        dict_avg["Solparm"] = [Solpar_avg, Solpar_avgstd]
        fnamegnu = os.path.join(diroutput, "gnuplot_ecoh.gnu")
        self._gnuplot_template_dimensions(fnamegnu, dict_avg)
        #
        # # Template plots distributions
        # fnamegnu = os.path.join(diroutput, "gnuplot_distributions.gnu")
        # self._gnuplot_template_distributions(fnamegnu, dict_avg)
        #
        # # Template plots characteristic ratio
        # fnamegnu = os.path.join(diroutput, "gnuplot_charratio.gnu")
        # self._gnuplot_template_charratio(fnamegnu, dict_avg)




    # =============================================================================
    @staticmethod
    def call_subprocess(cmd, cwd=None):
        """
        Execute a shell command `cmd`, using the environment `env`
        and wait for the results.
        """
        try:
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd)
            rs = p.communicate()
        except Exception as e:
            rs = None
            msg1 = "Subprocess fails with error: {}\n".format(e)
            msg2 = "Command: {}\n".format(cmd)
            print(msg1 + msg2)

        # It is a tuple (stdout_data, stderr_data). The data will
        # be strings if streams were opened in text mode; otherwise, bytes
        return rs


    # ==================================================================================
    def _mdp_rerun(self, dt=0.001):

        line = "; RUN CONTROL PARAMETERS\n"
        line += "integrator               = md\n"
        line += "dt                       = {}\n".format(dt)
        line += "nsteps                   = 0\n"
        line += "tinit                    = 0\n"
        line += "init-step                = 0\n"
        line += "comm-mode                = Linear\n"
        line += "nstcomm                  = 100\n"
        line += "comm-grps                = system\n"
        line += "\n"
        line += "; OUTPUT CONTROL OPTIONS\n"
        line += "nstxout                  = 100\n"
        line += "nstvout                  = 0\n"
        line += "nstfout                  = 0\n"
        line += "nstlog                   = 0\n"
        line += "nstcalcenergy            = 100\n"
        line += "nstenergy                = 100\n"
        line += "nstxout-compressed       = 1\n"
        line += "compressed-x-precision   = 0\n"
        line += "\n"
        line += "; NEIGHBORSEARCHING PARAMETERS\n"
        line += "nstlist                  = 5\n"
        line += "pbc                      = xyz\n"
        line += "rlist                    = 1.2\n"
        line += "\n"
        line += "; OPTIONS FOR ELECTROSTATICS AND VDW\n"
        line += "cutoff-scheme            = verlet\n"
        line += "\n"
        line += "coulombtype              = PME\n"
        line += "rcoulomb                 = {0:.1f}\n".format(self._coulombcutoff)
        line += "pme_order                = 4\n"
        line += "fourierspacing           = 0.16\n"
        line += "\n"
        line += "vdw_type                 = Cut-off\n"
        line += "rvdw                     = {0:.1f}\n".format(self._vdvcutoff)
        line += "rvdw_switch              = {0:.1f}\n".format(self._vdvcutoff)
        line += "\n"
        if self._dispersioncorrect:
            line += "DispCorr = EnerPres\n"
        line += "\n"
        line += "; OPTIONS FOR TEMPERATURE COUPLING ALGORITHMS\n"
        line += "tcoupl                   = V-rescale\n"
        line += "tc-grps                  = System\n"
        line += "tau_t                    = 0.5\n"
        line += "ref_t                    = 400\n"
        line += "\n"
        line += "; OPTIONS FOR PRESSURE COUPLING ALGORITHMS\n"
        line += ";Pcoupl                  = berendsen\n"
        line += ";Pcoupl                  = Parrinello-Rahman\n"
        line += "Pcoupl                   = C-rescale\n"
        line += "Pcoupltype               = isotropic\n"
        line += "tau_p                    = 5.0\n"
        line += "compressibility          = 4.5e-5\n"
        line += "ref_p                    = 1.0\n"
        line += "\n"
        line += "; GENERATE VELOCITIES FOR STARTUP RUN\n"
        line += ";gen_vel                  = yes\n"
        line += ";gen_temp                 = 473\n"
        line += ";gen_seed                 = 173529\n"
        line += "continuation = yes\n"
        if dt > 0.001:
            line += '; OPTIONS FOR BONDS\n'
            line += '; -----------------------------------------------------------------------------------------\n'
            line += 'constraints              = h-bonds    ; [1] water is done via shake\n'
            line += '; Type of constraint algorithm\n'
            line += 'constraint-algorithm     = Lincs      ; default\n'
            line += 'Shake-SOR                = no         ; default\n'
            line += 'shake-tol                = 0.0001     ; default also same as in [1]\n'
            line += 'lincs-order              = 4          ; default\n'
            line += 'lincs-iter               = 1          ; default\n'
            line += 'lincs-warnangle          = 90         ; default\n'
            line += 'morse                    = no         ; default\n'


        with open("rerun.mdp", 'w') as fmdp:
            fmdp.writelines(line)

    # ==================================================================================
    def _createtopology_for_sc(self, name):

        with open(self._topology, 'r') as ftopo:
            topo_lines = ftopo.readlines()

        ffpath = os.path.realpath(self._topology)
        ffpath = os.path.split(ffpath)[0]

        for idx, iline in enumerate(topo_lines):
            if re.findall(r"\[ *molecules *\]", iline):
                idx_cut = idx

        fullpath_sctopo = os.path.join(ffpath,"polymer_SC.top")
        with open(fullpath_sctopo, "w") as fsctopo:
            fsctopo.writelines(topo_lines[:idx_cut-1])
            fsctopo.writelines("\n[ molecules ]\n")
            name = name.replace("\"","")
            fsctopo.writelines(name + " 1\n")

        return fullpath_sctopo

    # ==================================================================================
    @staticmethod
    def _clean_files(imol):

        # Clean isolated files
        if imol is not None:

            os.remove("SC_traj_ch_{0:04d}_out.trr".format(imol))
            os.remove("SC_traj_ch_{0:04d}_out.log".format(imol))
            os.remove("SC_traj_ch_{0:04d}_out.edr".format(imol))
            os.remove("SC_traj_ch_{0:04d}.gro".format(imol))
            os.remove("SC_traj_ch_{0:04d}.tpr".format(imol))
            os.remove("SC_traj_ch_{0:04d}.xtc".format(imol))
            os.remove("tmp_{0:04d}.tpr".format(imol))

        # Clean xvg files
        else:
            ll = glob.glob("*.xvg")
            for item in ll:
                os.remove(item)
            os.remove("Full_traj_out.edr".format(imol))
            os.remove("Full_traj_out.log".format(imol))
            os.remove("Full_traj_out.trr".format(imol))
            os.remove("full.gro".format(imol))
            os.remove("mdout.mdp".format(imol))
            os.remove("trj_full.xtc".format(imol))
            os.remove("rerun.mdp".format(imol))
            os.remove("TMPnewtopo.tpr".format(imol))

            ll = glob.glob("#*")
            for item in ll:
                os.remove(item)

    # #########################################################################
    @staticmethod
    def _gnuplot_template_dimensions(filenamegnu, dict_avg):

        basedir, file = os.path.split(filenamegnu)
        if basedir == "":
            basedir = "./"

        # Defaults ===========================
        ps=[1.4, 1.0]
        lw=[2.0, 1.0]
        colors=["black", "blue", "red", "orange"]
        dt=[1, 2, 3, 4]
        pt_empty=[4, 6, 8, 12]
        pt_full = [5, 7, 9, 13]
        # Defaults ===========================
        d = defaultdict()
        d["Ecoh"] = {"fname":os.path.join(basedir,"CED_solubility.dat"), "cols":[1, 3],
                   "labels": ["t (ps)", "<Ecoh> (kJ/mol)"]}
        d["Molvol"] = {"fname":os.path.join(basedir,"CED_solubility.dat"), "cols":[1, 4],
                   "labels": ["t (ps)", "<Molvol> (cm^3/mol)"]}
        d["CED"] = {"fname":os.path.join(basedir,"CED_solubility.dat"), "cols":[1, 5],
                        "labels": ["t (ps)", "<CED> (J/cm^3)"]}
        d["Solparm"] = {"fname":os.path.join(basedir,"CED_solubility.dat"), "cols":[1, 6],
                        "labels": ["t (ps)", "Solubility Parameter (J/cm^3)^0.5"]}

        nfiles = len(d)
        if nfiles == 1:
            dim_plot = [600, 400, 1, 1]
        elif nfiles == 2:
            dim_plot = [1000, 600, 1, 2]
        else:
            dim_plot = [1200, 1200, 2, 2]

        line = 'reset\n'
        line += 'set term wxt 1 enhanced dashed size {},{} font "Arial,10"\n'.format(dim_plot[0], dim_plot[1])
        line += 'set multiplot layout {},{}\n'.format(dim_plot[2], dim_plot[3])

        idx = 1
        for ikey, value in d.items():

            line += 'set xlabel "{}"\n'.format(d[ikey]["labels"][0])
            line += 'set ylabel "{}"\n'.format(d[ikey]["labels"][1])
            line += 'set grid\n'
            line += 'set style fill transparent solid 0.5 noborder\n\n'
            if ikey == "Ecoh":
                title = "Cohesive energy"
            elif ikey == "Molvol":
                title = "Molar volume"
            elif ikey == "CED":
                title = "Cohesive Energy Density"
            elif ikey == "Solparm":
                title = "Solubility Parameter"
            line += 'set title "{}"\n'.format(title)
            line += 'p "{}" u {}:{} w l notitle lc "{}" lw {} dt {},\\\n'.\
                format(d[ikey]["fname"], d[ikey]["cols"][0], d[ikey]["cols"][1],
                       colors[idx-1], lw[0], dt[0] )
            line += '  {} lc "{}" lw 3 dt 2 notitle,\\\n'.format(dict_avg[ikey][0], colors[idx-1])
            line += '  {} with filledcurves y1={} lt 1 lc "grey" notitle, ' \
                    '{} with filledcurves y1={} lt 1 lc "grey" notitle\n\n'.\
                format(dict_avg[ikey][0]-dict_avg[ikey][1], dict_avg[ikey][0],
                       dict_avg[ikey][0]+dict_avg[ikey][1], dict_avg[ikey][0])
            idx += 1
        line += "\nunset multiplot"

        with open(filenamegnu, "w") as f:
            f.writelines(line)
