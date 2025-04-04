import math
import datetime
import os
import subprocess
import numpy as np
import polyanagro as pag

class MSD(pag.Calculations):

    r"""Class to calculate the mean square displacement"""

    # =============================================================================================
    def __init__(self, trj, isnotjump, list_atoms = None,
                 shift_atoms_from_center= None, method="multiple_window",
                 number_of_block_elements = 10,
                 number_of_blocks = 10, outputname="msd_atom.dat", start=0, end=-1, step=1, logger=None):

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

        implemented_methods = ["multiple_window", "classic_py", "classic_opt_cython",
                               "msd_fftw3", "classic_multitau", "msd_fftw3_fast"]

        self._trajectory = trj
        self._dt = self._trajectory.universe.trajectory.dt #in ps
        self._isnotjump = isnotjump
        self._logger = logger
        self._filename_msd = outputname
        self._start = start
        if end == -1:
            self._end = len(self._trajectory.universe.trajectory)
        else:
            self._end = end
        self._step = step

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

        # ====================== MULTIPLE WINDOW METHOD ===================
        if method == "multiple_window":

            self.max_number_of_blocks = number_of_blocks
            self.number_of_block_elements = number_of_block_elements
            self.number_particles = len(self.s_atoms)

            self.block_data = np.zeros([self.max_number_of_blocks,
                                        self.number_particles,
                                        self.number_of_block_elements, 3], dtype=float)

            self.block_length = np.zeros (self.max_number_of_blocks, dtype=int)
            self.msdcount = np.zeros([self.max_number_of_blocks, self.number_of_block_elements], dtype = float)
            self.msdav = np.zeros([self.max_number_of_blocks, self.number_of_block_elements], dtype = float)

            self.msdav_x = np.zeros([self.max_number_of_blocks, self.number_of_block_elements], dtype = float)
            self.msdav_y = np.zeros([self.max_number_of_blocks, self.number_of_block_elements], dtype = float)
            self.msdav_z = np.zeros([self.max_number_of_blocks, self.number_of_block_elements], dtype = float)

            self.msd_internal_multiplewindow_python()

        elif method == "classic_py":

            self.msd_classic_py()

        elif method == "classic_opt_cython":

            self.msd_classic_opt_cython()

        elif method == "classic_multitau":

            self.msd_classic_multitau()

        elif method== "msd_fftw3":

            self.msd_fftw3_cython()

        elif method== "msd_fftw3_fast":

            self.msd_fftw3_fast()

        else:
            print ("MSD cannot be calculated.")
            print ("Method: {} is not implemented".format(method))
            print ("Use one of the following methods: {}".format(implemented_methods))
            exit()

        pass

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
        weedf = 2  # number of starting points of independent time series

        m = "\t\t Calculating MSD for {} atoms/particles in {} frames.\n".format(n_atoms, n_frames)
        m1= "\t\t "+"*"*len(m)+"\n"
        print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
        for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
            positions[i] = self._s_atoms.positions

        pag.msd_all_c(positions, n_mol, self._dt, weedf)


    # # =============================================================================================
    # def msd_fftw3_cython(self):
    #
    #     n_atoms = len(self._s_atoms)
    #     n_frames = len(self._trajectory.universe.trajectory[self._start:self._end:self._step])
    #     positions = np.zeros((n_frames, n_atoms, 3), dtype=np.float64)
    #
    #     m = "\t\t Calculating MSD for {} atoms/particles in {} frames (FFTW3 algorithm).\n".format(n_atoms, n_frames)
    #     m1= "\t\t "+"*"*len(m)+"\n"
    #     print(m1+m+m1) if self._logger is None else self._logger.info(m1+m+m1)
    #     for i, ts in enumerate(self._trajectory.universe.trajectory[self._start:self._end:self._step]):
    #         positions[i] = self._s_atoms.positions
    #
    #     msd_avg = pag.msd_fftw3_ext(positions)
    #
    #     # Print
    #     filename = "msd_fftw3.dat"
    #     with open(filename, 'w') as f:
    #         for i in range(n_frames):
    #             line = "{0:.3f}  {1:.6f}\n".format(i*self._dt, msd_avg[i])
    #             f.writelines(line)

    # print(msd_avg)

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
