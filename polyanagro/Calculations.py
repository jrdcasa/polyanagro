import numpy as np
import polyanagro as pag
import topology
import os
from collections import defaultdict

# ===============================================================
class Calculations(object):

    __slots__ = ["_trajectory", "_dt", "_stride", "_nmols_array",
                 "_l_neigh_array", "_freq", "_coords_unwrap", "_iframe",
                 "_logger"]
    # #########################################################################
    def __init__(self, trj, dt=1, stride=1, freq=50, logger=None):

        """
        Calculations constructor

        ``Parameters``:
            * ``trj`` (Trajectory object) : Trajectory object containing Universe, Topology and others
            * ``dt`` (float) : Time step (default in ps)
            * ``stride`` (int) : Stride between frames
            * ``freq`` (int) : Frequency

        ``Attributes``:
            * ``self._trajectory`` (Trajectory object) : Trajectory object to perform the calculations
            * ``self._dt`` (float) : Time step (default in ps)
            * ``self._stride`` (int) : Stride between frames
            * ``self._freq`` (int) : Frequency
            * ``self._nmols_array`` (ndarray-int32 [nmols, maxnumberofatomsinmol]) : Index of atoms for each molecule (chain)
                [ [ 0 1 2 ...] [24 25 26 ...] [48 49 50 ...]]
            * ``self._l_neigh_array`` (ndarray-int32 [natoms, 3]) : Number of neighbours for each atom. The value -1 represents not neighbour in that position
            * ``self._coords_unwrap`` (ndarray-float32 [natoms, 3]) : Unwrapped coordinates for the frame iframe
            * ``self._iframe`` (int) : Frame in which the trajectory is.
            * ``self._logger`` (int) : Frame in which the trajectory is.

        """

        self._trajectory = trj
        self._dt = dt  # Timestep in picoseconds
        self._stride = stride
        self._nmols_array, self._l_neigh_array = self._trajectory.topology.get_array_mols_neigh()
        self._freq = freq
        self._coords_unwrap = None
        self._iframe = None
        self._logger = logger

    # #########################################################################
    def _unwrap_coordinates(self, iframe):

        """
        This method unwraps the PBC coordinates for the frame ``iframe``

        ``Parameters``:
            * ``iframe`` (int) : Frame

        ``Return``:
            * ``None``

        """

        coords_t0_wrapped = self._trajectory.universe.trajectory[iframe].positions
        box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]

        self._coords_unwrap = topology.unwrap(coords_t0_wrapped, self._nmols_array, self._l_neigh_array,
                                         box_dimensions, iframe=iframe)

        self._iframe = iframe

        return None

    # #########################################################################
    @staticmethod
    def _calc_avg(data, fraction_trj):

        first_frame = int(data[0][0])
        last_frame = int(data[-1][0])
        total_frames = int(last_frame - first_frame + 1)
        delta = int(total_frames / data.shape[0])

        iframe_begin_avg = (1 - fraction_trj) * total_frames
        frames_avg = (total_frames - iframe_begin_avg + 1)

        avg = np.mean(data[iframe_begin_avg:total_frames:delta,1])
        std_avg = np.std(data[iframe_begin_avg:total_frames:delta, 1])

        std = 0.0
        N = 0
        for i in range(iframe_begin_avg, total_frames, delta):
            std += data[i,2]**2
            N += 1
        std = np.sqrt(std/N)

        return avg, std, std_avg

    # #########################################################################
    @staticmethod
    def _extract_data_for_avg(filename, cols):

        icol1 = cols[0]
        icol2 = cols[1]
        icol3 = cols[2]

        with open(filename, "r") as f:
            lines = []
            while True:
                l = f.readline()
                if not l:
                    break
                if l.find("#") != -1:
                    continue
                idx = int(l.split()[icol1])
                value = float(l.split()[icol2])
                if icol3 is not None:
                    std = float(l.split()[icol3])
                else:
                    std = 0.0
                lines.append([idx, value, std])
            data_array = np.array(lines)

        return data_array

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
        d["Rg"] = {"fname":os.path.join(basedir,"Rg.dat"), "cols":[2, 3, 4],
                   "labels": ["t (ps)", "<Rg^2> (angstroms)"]}
        d["Ree"] = {"fname":os.path.join(basedir,"Ree.dat"), "cols":[2, 3, 4],
                   "labels": ["t (ps)", "<Ree^2> (angstroms)"]}
        d["Ree2Rg2"] = {"fname":os.path.join(basedir,"Ree2Rg2.dat"), "cols":[2, 3, 4],
                        "labels": ["t (ps)", "<Ree^2>/<Rg^2>"]}
        d["Cn"] = {"fname":os.path.join(basedir,"Cn.dat"), "cols":[2, 3, 4],
                        "labels": ["t (ps)", "Cn"]}

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
            if ikey == "Rg":
                title = "Radius of gyration"
            elif ikey == "Ree":
                title = "End-to-End distance"
            elif ikey == "Ree2Rg2":
                title = "<R_{ee}^2>/<R_g^2>"
            elif ikey == "Cn":
                title = "Characteristic ratio"
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

