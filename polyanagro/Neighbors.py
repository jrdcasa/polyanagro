import numpy as np
import os
from collections import defaultdict
from MDAnalysis import transformations as trans
import MDAnalysis as mda
from topology.internal_coordinates import distance_array

class Neighbors(object):

    __slots__ =  ["_com_list", "_trajectory", "_nmols_array", "_l_neigh_array", "_log", "_nmols",
                  "_com_filename_pdb", "_neighbor_sphere"]

    # #######################################################################
    def __init__(self, trj, calc_com=True, writetcl=True, log=None):

        self._com_list = []
        self._trajectory = trj
        self._nmols_array, self._l_neigh_array = self._trajectory.topology.get_array_mols_neigh()
        self._nmols = len(self._nmols_array)
        self._log = log
        self._com_filename_pdb = None
        self._neighbor_sphere = defaultdict(list)

        # Unwrap trajectory:
        ag = self._trajectory.universe.atoms
        transform = trans.unwrap(ag)
        self._trajectory.universe.trajectory.add_transformations(transform)

        # COM must be done with the coordinates unwrapped
        if calc_com:
            self._com()

        # Write a tmp pdb and redefine the universe in the the tajectory.
        # This is neccesary because MDAnalysis does not allow to make
        # two transformation to a given universe.
        ag.write("wrapped.pdb")
        self._trajectory.universe = mda.Universe("wrapped.pdb")
        ag = self._trajectory.universe.atoms
        transform = trans.wrap(ag)
        self._trajectory.universe.trajectory.add_transformations(transform)
        ag.write("wrapped.pdb")


    # #######################################################################
    def _com(self):

        for imol in self._nmols_array:
            seltxt = "index "
            for iat in imol:
                seltxt += "{0:d} ".format(iat)
            g = self._trajectory.universe.select_atoms(seltxt)
            self._com_list.append(g.center_of_mass())

        #Wrap COM
        tmp_com = []
        box_x = self._trajectory.universe.dimensions[0]
        box_y = self._trajectory.universe.dimensions[1]
        box_z = self._trajectory.universe.dimensions[2]
        for icom in self._com_list:
            x = icom[0] - np.floor(icom[0]/box_x)*box_x
            y = icom[1] - np.floor(icom[1]/box_y)*box_y
            z = icom[2] - np.floor(icom[2]/box_z)*box_z
            tmp_com.append([x, y, z])
        self._com_list = tmp_com


    # #######################################################################
    def _write_com_pdb(self, filename_com_pdb="coords_com.pdb", wrap=False):

        """
        Write a pdb file to check the center of mass.
        Adapted from MDAnalysis software (https://www.mdanalysis.org/)
        """

        self._com_filename_pdb = filename_com_pdb

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        a = self._trajectory.universe.dimensions[0]
        b = self._trajectory.universe.dimensions[1]
        c = self._trajectory.universe.dimensions[2]
        alpha = self._trajectory.universe.dimensions[3]
        beta = self._trajectory.universe.dimensions[4]
        gamma = self._trajectory.universe.dimensions[5]

        atom_kind_molecule_label = None
        spacegroup = "P -1"
        zvalue = 1
        radtodeg = 180/np.pi
        with open(filename_com_pdb, 'w') as fpdb:

            fpdb.write(fmt['REMARK'].format('Created with polyanogro(J.Ramos)'))
            fpdb.write(fmt['CRYST1'].format(a, b, c, alpha, beta, gamma,
                                            spacegroup, zvalue))

            for idx in range(self._nmols):

                resname = "{0:03d}".format(idx)

                if atom_kind_molecule_label is None:
                    fpdb.write(fmt['HETATM'].format(
                        serial=idx+1,
                        name='C',
                        altLoc=" ",
                        resName=resname,
                        chainID=" ",
                        resSeq=idx,
                        iCode=" ",
                        pos=[i for i in self._com_list[idx]],
                        occupancy=1.0,
                        tempFactor=1.0,
                        segID="    ",
                        element='C'
                    ))


            fpdb.write('END\n')

        self._write_tcl_atoms_com()

    # # #######################################################################
    def find_sphere_com(self, radius=None):

        """
        Find other com's inside a sphere of this radius
        Args:
            radius: in angstroms

        Returns:

        """

        a = self._trajectory.universe.dimensions[0]
        b = self._trajectory.universe.dimensions[1]
        c = self._trajectory.universe.dimensions[2]

        if radius is None:
            radius = min([0.5*a, 0.5*b, 0.5*c])

        if 2.0*radius > a or 2.0 * radius > b or 2.0 * radius > c:
            m = "Radius {} cannot be greater than half the box ({},{},{})\n".format(radius, a, b, c)
            print(m) if self._log is None else self._log.error(m)
            exit()

        ncoms = len(self._com_list)

        for ires, icom in enumerate(self._com_list):
            ref = np.zeros((1,3), dtype=np.float64)
            ref[0, :] = np.array(icom, dtype=np.float64)
            points = np.zeros((ncoms, 3), dtype=np.float64)
            idx = 0
            for jres, jcom in enumerate(self._com_list):
                points[idx, :] = np.array(jcom)
                idx += 1
            d = distance_array(ref, points, openmp=True)
            for jres in range(ncoms):
                if jres == ires:
                    continue
                if d[0][0][jres] < radius:
                    self._neighbor_sphere[ires].append(jres)

        print (self._neighbor_sphere)

    # #######################################################################
    def _write_tcl_atoms_com(self, filename_tcl="vmd_com.tcl"):

        lines = "proc newRep { sel type color rep imol} {\n"
        lines += "    mol selection $sel\n"
        lines += "    mol representation $type\n"
        lines += "    mol addrep $imol\n"
        lines += "    mol showrep $imol $rep on\n"
        lines += "    mol modcolor $rep $imol $color\n"
        lines += "}\n"
        lines += "\n"
        lines += "set dir \"{}\"\n".format(os.path.split(filename_tcl)[0])
        lines += "\n"
        lines += "display projection orthographic\n"
        lines += "axes location off\n"
        lines += "color Display Background white\n"
        lines += "display depthcue off\n"
        lines += "\n"

        lines += "\n"
        lines += "mol new {} type pdb\n".format(self._trajectory.trjpath[0])
        lines += "set imol1 [molinfo top]\n"
        lines += 'mol delrep 0 $imol1\n'
        lines += 'set rep1 0\n'
        lines += 'newRep "all" "CPK" "Name" $rep1 $imol1\n'
        lines += '\n'
        lines += "mol new {} type pdb\n".format(self._com_filename_pdb)
        lines += "set imol2 [molinfo top]\n"
        lines += 'mol delrep 0 $imol2\n'
        lines += 'set rep2 0\n'
        lines += 'newRep "all" "vdw 0.6" "Index" $rep2 $imol2\n'
        lines += '\n'
        lines += "pbc box\n"
        lines += '\n'

        with open(filename_tcl, 'w') as ftcl:
            ftcl.writelines(lines)
