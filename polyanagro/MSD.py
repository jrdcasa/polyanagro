import math
import numpy as np

class MSD(object):

    r"""Class to calculate the mean square displacement"""

    # =============================================================================================
    def __init__(self, universe, list_atoms = None,
                 shift_atoms_from_center= None, method="multiple_window",
                 number_of_block_elements = 10,
                 number_of_blocks = 10):

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

        implemented_methods = ["multiple_window", "classical"]

        self.universe = universe
        self.dt = universe.trajectory.dt #in ps

        if list_atoms is not None:
            print("list_atoms")
            self.s_atoms = list_atoms
        elif shift_atoms_from_center is not None:
            print("shift_atoms_center")
            self.s_atoms = None
        else:
            print("All atoms are considered")
            self.s_atoms = universe.select_atoms("index 0")
        print("msd internal ****")

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
        for ts in self.universe.trajectory:

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
        with open("msd_atom.dat", 'w') as f:

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
