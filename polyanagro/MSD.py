import math
import datetime
import os
import subprocess
import numpy as np
import polyanagro as pag
import psutil

class MSD(object):

    r"""Class to calculate the mean square displacement"""

    # =============================================================================================
    def __init__(self, trj, isnotjump, list_atoms = None,
                 shift_atoms_from_center= None, method="msd_fftw3_cython",
                 startingfactor=None, com=0,
                 outputname="msd_atom.dat", start=0, end=-1, step=1, logger=None):

        r"""Constructor of the msd object


        Parameters
        ----------
        universe : Universe object from MDAnalysis
        list_atoms: list
            list of atom to take into account in the calculation
        shift_atoms_from_center: integer
            Offset of the atoms from and to the center atom in a polymer chain.
            Example: shift_atoms_from_center = 4 (center atom is calculated, taken into account only the
                     backbone atoms)
                       (36) -- (37) -- (38) -- (39) -- [40] -- (41) -- (42) -- (43) -- (44)

        The mean square displacement (msd) on the selected atoms is calculated. The selection of atoms is
        performed using the parameter **list_atoms** or **shift_atoms_from_center**, otherwise all atoms
        are ued in the calculations. IN the case of the both parameters are not **None**, the **list_atoms**
        is prevalent.

        """

        implemented_methods = ["classic_py", "classic_opt_cython",
                               "classic_multitau", "msd_fftw3_fast", "msd_fftw3_cython"]

        self._trajectory = trj
        self._isnotjump = isnotjump
        self._logger = logger
        self._filename_msd = outputname
        self._start = start
        self._freq = 100
        self._startingfactor = startingfactor
        if end == -1:
            self._end = len(self._trajectory.universe.trajectory)
        else:
            self._end = end
        self._step = step
        self._positions_signal = None
        self._dt = self._trajectory.universe.trajectory.dt * self._step  # in ps
        if list_atoms is not None:
#            print("list_atoms")
            self._s_atoms = list_atoms
        elif shift_atoms_from_center is not None:
#            print("shift_atoms_center")
            self._s_atoms = None
        else:
#            print("All atoms are considered")
            self._s_atoms = self._trajectory.universe.select_atoms("all")

        if not self._isnotjump:
            m = "\t\t The trajectory must be unwrapped and without jumps."
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        self._n_atoms = len(self._s_atoms)
        self._n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
        self._n_dim = 3
        self._n_mols = len(self._trajectory.topology._nmols)
        self._com = com
        self._nmols_array, self._l_neigh_array = self._trajectory.topology.get_array_mols_neigh()

        if method == "classic_py":

            self.msd_classic_py()

        elif method == "classic_opt_cython":

            self.msd_classic_opt_cython()

        elif method == "classic_multitau":

            self.msd_classic_multitau()

        elif method== "msd_fftw3_cython":

            self.msd_fftw3_cython()

        elif method== "msd_fftw3_fast":

            self.msd_fftw3_fast()

        else:
            print ("MSD cannot be calculated.")
            print ("Method: {} is not implemented".format(method))
            print ("Use one of the following methods: {}".format(implemented_methods))
            exit()


    # =============================================================================================
    def msd_internal_multiplewindow_python(self):

        r"""Calculate mean-square-displacement using the multiple window method

        The method of blocks based on multiple windows is described in [Dubbeldam2009]_
        References
        ----------
        .. [Dubbeldam2009] Dubbeldam D., Ford D.C., Ellis D.E. and Snurra R.Q. 2009
           A new perspective on the order-n algorithm for computing correlation functions
           *Molecular Simulation* 35(12-13) 1084-1097.


        """

        # Sampling MSD along the trajectory
        for ts in self._trajectory.universe.trajectory:

            iframe = ts.frame

            # Determine current number of blocks to be updated
            number_of_blocks = 1
            i = int(iframe/self.number_of_block_elements)
            while i != 0:
                number_of_blocks += 1
                i = int(i / self.number_of_block_elements)

            number_of_blocks = min(number_of_blocks, self.max_number_of_blocks)

            #loop over all the blocks to test
            #which blocks need sampling
            for current_block in range(0, number_of_blocks):

                #Update current block info?
                m = iframe % math.pow(self.number_of_block_elements, current_block)
                if m == 0:

                    # Increase the current block length
                    self.block_length[current_block] += 1
                    # Compute the current length of the block, limited to size ’NumberOfBlockElements’
                    currentblocklength=min(self.block_length[current_block],self.number_of_block_elements);

                    # Loop over the selected particles, atoms, molecules ...
                    for k in self.s_atoms:
                        index_atom = k.ix
                        # shift to the left, set last index to the correlation value
                        for ielem in range(1, self.number_of_block_elements):
                            self.block_data[current_block, index_atom, ielem-1] = \
                                self.block_data[current_block, index_atom, ielem]

                        self.block_data[current_block, index_atom, self.number_of_block_elements-1] = ts.positions[index_atom]

                        # get the origin, take into account that blocks can be partially filled
                        index_origin = self.number_of_block_elements - currentblocklength
                        origin = self.block_data[current_block, index_atom, index_origin]

                        for ielem in range(0, currentblocklength):

                            a0 = self.block_data[current_block, index_atom, index_origin+ielem, 0] - origin[0]
                            self.msdav_x[current_block][ielem] += a0 * a0
                            a1 = self.block_data[current_block, index_atom, index_origin+ielem, 1] - origin[1]
                            self.msdav_y[current_block][ielem] += a1 * a1
                            a2 = self.block_data[current_block, index_atom, index_origin+ielem, 2] - origin[2]
                            self.msdav_z[current_block][ielem] += a2 * a2

                            self.msdcount[current_block][ielem] += 1.0
                            self.msdav[current_block][ielem] += a0 * a0 + \
                                                                a1 * a1 + \
                                                                a2 * a2


        # Print
        with open(self._filename_msd, 'w') as f:

            f.writelines("# column 1: time [ps]\n")
            f.writelines("# column 2: msd xyz [A^2]\n")
            f.writelines("# column 3: msd x [A^2]\n")
            f.writelines("# column 4: msd y [A^2]\n")
            f.writelines("# column 5: msd z [A^2]\n")
            f.writelines("# column 6: number of samples [-]\n")

            for current_block in range(0,number_of_blocks):
                currentblocklength = min(self.block_length[current_block], self.number_of_block_elements)

                if current_block == 0:
                    line = "{0:.3f}  {1:.6f}  {2:.1f}\n".format(0.0,
                                                                0.0,
                                                                0.0)
                    f.writelines(line)


                for ielem in range(1, currentblocklength):

                    itime = ielem*self.dt*math.pow(self.number_of_block_elements, current_block)

                    # Write time-index
                    if self.msdcount[current_block,ielem] > 0.0:
                        fac = 1.0 / self.msdcount[current_block][ielem]
                        line = "{0:.3f}  {1:.6f}  {2:.6f}  {3:.6f}  {4:.6f}  {5:.1f}\n".format(itime,
                                                               fac*self.msdav[current_block,ielem],
                                                               fac*self.msdav_x[current_block,ielem],
                                                               fac*self.msdav_y[current_block,ielem],
                                                               fac*self.msdav_z[current_block,ielem],
                                                               self.msdcount[current_block,ielem])
                        f.writelines(line)

    # =============================================================================================
    def msd_classic_py(self):

        n_atoms = len(self._s_atoms)
        n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
        positions = np.zeros((n_frames, n_atoms, 3))
        msd = np.zeros(n_frames)
        count = np.zeros(n_frames)

        m = "\t\t Calculating MSD for {} atoms/particles in {} frames.\n".format(n_atoms, n_frames)
        m1= "\t\t "+"*"*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
        for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
            positions[i] = self._s_atoms.positions

        indx_iframe = 0
        start_time = datetime.datetime.now()
        for i in range(n_frames):

            if indx_iframe % 100 == 0:
                msg = "\tFrame {0:9d} of {1:9d}".format(indx_iframe, n_frames)
                mid_time = datetime.datetime.now()
                elapsed_time = mid_time - start_time
                msg += "\ttime: {0:s} seconds".format(str(elapsed_time.total_seconds()))
                print(msg) if self._logger is None else self._logger.info(msg)

            for j in range(i, n_frames):
                disp = positions[j] - positions[i]
                sq_disp = np.sum(disp ** 2, axis=1)
                msd[j - i] += np.sum(sq_disp)
                count[j - i] += n_atoms
            indx_iframe += 1

        msd /= count
        times = np.arange(n_frames) * self._dt

       # Print
        with open(self._filename_msd, 'w') as f:
            for i, imsd in enumerate(msd):
                line = "{0:.3f}  {1:.6f}\n".format(i*self._dt, imsd)
                f.writelines(line)

        return times, msd

    # =============================================================================================
    def msd_classic_opt_cython(self):

        n_atoms = len(self._s_atoms)
        n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
        positions = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
        msd = np.zeros(n_frames)
        count = np.zeros(n_frames)

        m = "\t\t Calculating MSD for {} atoms/particles in {} frames.\n".format(n_atoms, n_frames)
        m1= "\t\t "+"*"*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
        for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
            positions[i] = self._s_atoms.positions

        filename = self._filename_msd.encode('utf-8')  # Convert str to bytes
        pag.msd_atomistic_cython(n_frames, n_atoms, 3, self._dt, positions, filename)

        filename = "parallel.dat".encode('utf-8')  # Convert str to bytes
        pag.msd_atomistic_opt_cython(n_frames, n_atoms, 3, self._dt, positions, filename)

        #return times, msd

    # =============================================================================================
    def msd_classic_multitau(self):

        n_atoms = len(self._s_atoms)
        n_mol = len(self._trajectory.topology._nmols)
        n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
        positions = np.zeros((n_frames, n_atoms, 3), dtype=np.float64)
        # Input parameter from command line
        weedf = self._startingfactor  # number of starting points of independent time series

        # Test if there is memory available in the system
        self._check_memory_multitau(weedf)

        m = "\t\t Calculating MSD for {} atoms/particles in {} frames.\n".format(n_atoms, n_frames)
        m1= "\t\t "+"*"*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
        for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
            positions[i] = self._s_atoms.positions

        pag.msd_all_c(positions, n_mol, self._dt, weedf, self._filename_msd)

        # COM MSD
        if self._com:

            m = "\t\t Calculating COM-MSD for {} particles in {} frames.\n".format(n_mol, n_frames)
            m1= "\t\t "+"*"*len(m)+"\n"
            print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)

            positions_com = np.zeros((n_frames, n_mol, 3), dtype=np.float64)
            mass = np.array(self._trajectory.topology.mass, dtype=np.float64)

            pag.msd_com_c(self._nmols_array, mass, positions, positions_com)

                # Write COM trj
            with open("com_trajectory.xyz", "w") as fxyz:
                iframe = 0
                indx_iframe = 0
                start_time = datetime.datetime.now()
                for iframe in range(0, n_frames):

                    fxyz.writelines("{}\n".format(n_mol))
                    fxyz.writelines("Frame {}\n".format(iframe))
                    for ich in range(0, n_mol):
                        line = "C {} {} {}\n".format(positions_com[iframe, ich, 0],
                                                   positions_com[iframe, ich, 1],
                                                   positions_com[iframe, ich, 2] )
                        fxyz.writelines(line)

            self._filename_msd = "msd_com.dat"
            pag.msd_all_c(positions_com, n_mol, self._dt, weedf, self._filename_msd)

    # =============================================================================================
    def msd_fftw3_cython(self):

        n_dim = 3
        # self._n_atoms = len(self._s_atoms)
        # self._n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
        self._positions_signal = np.zeros((self._n_atoms, n_dim, self._n_frames), dtype=np.float64)

        m = "\t Calculating MSD for {} atoms/particles in {} frames (FFTW3 algorithm).\n".\
            format(self._n_atoms, self._n_frames)
        m1= "\t "+"="*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
        s = datetime.datetime.now()
        self._freq = int(self._n_frames / (self._n_frames*0.001))

        # Test if there is memory available in the system
        self._check_memory_msdfftw3()

        for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
            self._positions_signal[:,:, i] = self._s_atoms.positions
            # Write info
            if i%self._freq == 0:
                elapsed_time = datetime.datetime.now() - s
                m = "\tIFRAME: {0:d} of {1:d} in {2:s} seconds".format \
                    (i, self._n_frames, str(elapsed_time.total_seconds()))
                print(m) if self._logger is None else self._logger.info(m)

        # Toy example:
        # self._positions_signal= np.zeros((3, 3, 4), dtype=np.float64)
        # self._n_atoms = 3
        # self._n_dim = 3
        # self._n_frames = 4
        # for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._n_frames:self._step]):
        #     for iatom in range(0, self._n_atoms):
        #         for idim in range(0, n_dim):
        #             self._positions_signal[iatom, idim, i] = self._s_atoms.positions[iatom, idim]

        msd = pag.msd_fftw3_cython(self._positions_signal, self._dt)

        # Print
        with open(self._filename_msd, 'w') as f:
            line = "#Timestep(ps) msd_x(A^2) msd_y(A^2) msd_z(A^2) msd(A^2)\n"
            f.writelines(line)
            for i in range(0, self._n_frames):
                line = "{0:.3f} ".format(i*self._dt)
                sum_msd_xyz = 0
                for k in range(0, n_dim):
                    tmp = msd[i+self._n_frames*k]
                    sum_msd_xyz += tmp
                    line += "{0:.6f} ".format(tmp)
                line += "{0:.6f}\n".format(sum_msd_xyz)
                f.writelines(line)

    # =============================================================================================
    def msd_fftw3_fast(self):

        n_atoms = len(self._s_atoms)
        n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])

        msd_fftw3_program = os.path.expanduser("~/.local/bin/msd")

        if not os.path.isfile(msd_fftw3_program):

            msg = "\t\t MSD_FFTW3 program is not installed\n"
            msg += "\t\t Please install from: https://github.com/jrdcasa/MeanSquareDisplacement.git\n"
            msg += "\t\t Please install from: https://github.com/RaulPPelaez/MeanSquareDisplacement\n"
            print(msg) if self._logger is None else self._logger.info(msg)
            return None

        m = "\t\t Calculating MSD for {} atoms/particles in {} frames (msd_fft3_fast).\n".format(n_atoms, n_frames)
        m1= "\t\t "+"*"*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)

        if not os.path.isfile("pos.dat"):
            indx_iframe = 0
            start_time = datetime.datetime.now()
            msg = "\t\t Writing trajectory to format needed in the msd_fftw3_fast function."
            print(msg) if self._logger is None else self._logger.info(msg)
            with open("pos.dat", "w") as f:
                for ts in self._trajectory.universe.trajectory:

                    if indx_iframe % 100 == 0:
                        msg = "\tFrame {0:9d} of {1:9d}".format(indx_iframe, n_frames)
                        mid_time = datetime.datetime.now()
                        elapsed_time = mid_time - start_time
                        msg += "\ttime: {0:s} seconds".format(str(elapsed_time.total_seconds()))
                        print(msg) if self._logger is None else self._logger.info(msg)

                    f.write("#\n")  # Frame separator
                    for atom in self._s_atoms:
                        f.write(f"{atom.position[0]} {atom.position[1]} {atom.position[2]}\n")

                    indx_iframe += 1
        else:
            msg = "\t\t pos.dat file is already present in the working directory."
            print(msg) if self._logger is None else self._logger.info(msg)

        # Command with parameters
        try:
            m = "\t\t Running MSD for {} atoms/particles in {} frames (msd_fft3_fast).\n".format(n_atoms, n_frames)
            start_time = datetime.datetime.now()
            print(m) if self._logger is None else self._logger.info(m)

            command = [os.path.expanduser("~/.local/bin/msd"), "-N", str(n_atoms), "-Nsteps", str(n_frames), "-dimensions", "3"]

            with open("pos.dat", "r") as infile, open("msd_fftw3.dat", "w") as outfile:
                subprocess.run(command, stdin=infile, stdout=outfile, text=True)

            mid_time = datetime.datetime.now()
            elapsed_time = mid_time - start_time
            msg = "\t\t MSD fftw3 time: {0:.2f} seconds".format(elapsed_time.total_seconds())
            print(msg) if self._logger is None else self._logger.info(msg)
        except subprocess.CalledProcessError as e:
            print(f"Error running MSD fftw3: {e.stderr}")

    # =============================================================================================
    def _check_memory_msdfftw3(self):

            msg = "\t\t Estimation of the memory neccesary \t\t\n"
            m1= "\t "+"*"*len(msg)+"\n"
            print(m1+msg+m1) if self._logger is None else self._logger.info(m1+msg+m1)

            process = psutil.Process(os.getpid())
            python_program_Gb = process.memory_info().rss / (1024 ** 3)
            positions_signal_bytes = self._positions_signal.nbytes / (1024 ** 3)
            S2c_signal_in_bytes = 8 * 2 * ((self._n_frames * 2) / 2 + 1) * self._n_atoms * 3 / (1024 ** 3)
            S2c_res_in_bytes = 8 * self._n_frames * self._n_atoms / (1024 ** 3)
            S2c_s2Averaged_in_bytes = 8 * self._n_frames * 3 / (1024 ** 3)
            S2c_sumAux_in_bytes = 8 * self._n_atoms * 3 / (1024 ** 3)
            S1c_s2Averaged_in_bytes = 8 * self._n_frames * 3 / (1024 ** 3)
            S1c_sumAux_in_bytes = 8 * self._n_frames * 3 / (1024 ** 3)
            ss = python_program_Gb + positions_signal_bytes + S2c_signal_in_bytes + S2c_res_in_bytes
            msg = "\t\t Number of atoms                              : {0:d} atoms\n".format(self._n_atoms)
            msg += "\t\t Number of frames                             : {0:d} frames\n".format(self._n_frames)
            msg += "\t\t Sample time                                  : {0:7.1f} ps\n".format(self._dt)
            msg += "\t\t Memory used by python program                : {0:7.3f} GB.\n".format(python_program_Gb)
            msg += "\t\t Memory used by signal (all atoms and frames) : {0:7.3f} GB.\n".format(positions_signal_bytes)
            msg += "\t\t   ======= Compute S2 =======\n"
            msg += "\t\t Memory used by c-FFTW3 (S2, signal, autocor) : {0:7.3f} GB.\n".format(S2c_signal_in_bytes)
            msg += "\t\t Memory used by c-FFTW3 (S2, res, autocor)    : {0:7.3f} GB.\n".format(S2c_res_in_bytes)
            msg += "\t\t Memory used by c-FFTW3 (S2Averaged )         : {0:7.3f} GB.\n".format(S2c_s2Averaged_in_bytes)
            msg += "\t\t Memory used by c-FFTW3 (sumAuxS2 )           : {0:7.3f} GB.\n".format(S2c_sumAux_in_bytes)
            msg += "\t\t Memory used by c-FFTW3 (S1Averaged )         : {0:7.3f} GB.\n".format(S1c_s2Averaged_in_bytes)
            msg += "\t\t Memory used by c-FFTW3 (sumAuxS1 )           : {0:7.3f} GB.".format(S1c_sumAux_in_bytes)

            print(msg) if self._logger is None else self._logger.info(msg)
            msg = "\t\t Total Memory estimation for FFTW3 MSD        : {0:7.3f} GB.".format(ss)
            m1= "\t "+"-"*len(msg)+"\n"
            print(m1+msg) if self._logger is None else self._logger.info(m1+msg)
            mem = psutil.virtual_memory()
            # Available memory in bytes
            free_bytes = mem.available
            free_gb = free_bytes / (1024 ** 3)

            msg = "\t\t Free memory                                  : {0:7.3f} GB.\n".format(free_gb)
            print(msg) if self._logger is None else self._logger.info(msg)

            if ss > free_gb:
                msg = "\n\t\t There is not free memory available. \n\t\t" \
                    " Increase memory and/or decrease the number of frames/atoms in the trajectory."
                print(msg) if self._logger is None else self._logger.error(msg)
                exit()

    # =============================================================================================
    def _check_memory_multitau(self, weedf):

            msg = "\t\t Estimation of the memory neccesary \t\t\n"
            m1= "\t "+"*"*len(msg)+"\n"
            print(m1+msg+m1) if self._logger is None else self._logger.info(m1+msg+m1)

            process = psutil.Process(os.getpid())
            python_program_Gb = process.memory_info().rss / (1024 ** 3)
            time_avg_bytes = self._n_frames * 9 / (1024 ** 3)
            times_bytes = weedf / (1024 ** 3)
            dxdydz_bytes = (self._n_atoms + self._n_mols) * self._n_frames * 8 / (1024 ** 3)
            ss = time_avg_bytes + times_bytes + dxdydz_bytes
            msg = "\t\t Number of atoms                              : {0:d} atoms\n".format(self._n_atoms)
            msg += "\t\t Number of molecules                          : {0:d} atoms\n".format(self._n_mols)
            msg += "\t\t Number of frames                             : {0:d} frames\n".format(self._n_frames)
            msg += "\t\t Sample time                                  : {0:7.1f} ps\n".format(self._dt)
            msg += "\t\t Memory used by python program                : {0:7.3f} GB.\n".format(python_program_Gb)

            print(msg) if self._logger is None else self._logger.info(msg)
            msg = "\t\t Total Memory estimation for FFTW3 MSD        : {0:7.3f} GB.".format(ss)
            m1= "\t "+"-"*len(msg)+"\n"
            print(m1+msg) if self._logger is None else self._logger.info(m1+msg)
            mem = psutil.virtual_memory()
            # Available memory in bytes
            free_bytes = mem.available
            free_gb = free_bytes / (1024 ** 3)

            msg = "\t\t Free memory                                  : {0:7.3f} GB.\n".format(free_gb)
            print(msg) if self._logger is None else self._logger.info(msg)

            if ss > free_gb:
                msg = "\n\t\t There is not free memory available. \n\t\t" \
                    " Increase memory and/or decrease the number of frames/atoms in the trajectory."
                print(msg) if self._logger is None else self._logger.error(msg)
                exit()
