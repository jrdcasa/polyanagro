import numpy as np
import polyanagro as pag
import topology

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




