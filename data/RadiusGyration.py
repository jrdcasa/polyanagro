import numpy as np
from polyanagro.Calculations import unwrap
from ext_libc.c_unwrap_openmp import calc_rg_openmp, calc_rg_openmp_massweigth

class Calculations:

    # #########################################################################
    def __init__(self, trj, top):

        self._rg2_frame = {}

        self._trj = trj
        self._top = top

        self._nmols_array, self._l_neigh_array = self._top.get_array_mols_neigh()


    # #########################################################################
    def calculate_frame(self, itrj, iframe):

        coords = itrj[iframe].positions
        box_dimensions = itrj[iframe].dimensions
        coords_unwrap = unwrap(coords, self._nmols_array,
                               self._l_neigh_array, box_dimensions,
                               iframe=iframe, write_xyz=False)

        rg2_avg = self._single_rg(coords_unwrap)

        self._rg2_frame[iframe] = rg2_avg

    # #########################################################################
    def _single_rg(self, cu, ismassweight=False):

        mass = np.array(self._top._mass, dtype=np.float32)
        nchains = len(self._nmols_array)
        rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)

        if ismassweight:
            calc_rg_openmp_massweigth(self._nmols_array, mass, cu, rgsq_ich_iframe)
        else:
            calc_rg_openmp(self._nmols_array, mass, cu, rgsq_ich_iframe)


        rg2_avg = 0.0
        for ich in range(nchains):
            rg2_avg += rgsq_ich_iframe[ich,3]

        return rg2_avg/nchains

    # #########################################################################
    def write_rg_tofile(self, iframe):



        pass

