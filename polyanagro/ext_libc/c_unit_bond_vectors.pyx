# import both numpy and the Cython declarations for numpy
import cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from copy import copy
cdef bint USED_OPENMP

cdef extern from "calc_bondvectors.c":
    float c_unit_bond_vectors(int natoms, int nchains, int nbonds,
                             int maxnbondsperch, int* bonds,
                             float* coords, int* iatch, float* uux,
                             float* uuy, float* uuz)
    float c_bond_bond_orientation(int natoms, int nchains, int nbonds,
                                  int maxnbondsperch, int* bonds,
                                  float* coords, int* iatch, float* cb)
    float c_odf_intra(int iframe, int natoms, int nchains, int nbonds,
                      int maxnbondsperch, int* bonds,
                      float* coords, int* iatch, float* uux,
                      float* uuy, float* uuz)
    void c_setup_odf_intra(int maxnbondsperch)
    void c_avg_write_odf_intra(int nframes, int maxnbondsperch, char* filename)

OPENMP_ENABLED = True if USED_OPENMP else False
# ========================================================================================
def unit_bond_vectors(nchains,
                      np.ndarray[int, ndim=2, mode="c"] all_bonds,
                      np.ndarray[float, ndim=2, mode="c"] coords,
                      np.ndarray[int, ndim=1, mode="c"] iatch,
                      np.ndarray[float, ndim=2, mode="c"] uux,
                      np.ndarray[float, ndim=2, mode="c"] uuy,
                      np.ndarray[float, ndim=2, mode="c"] uuz):

    nbonds = all_bonds.shape[0]
    natoms = coords.shape[0]
    nchains = uux.shape[0]
    maxnbondsperch = uux.shape[1]

    if iatch.shape[0] != natoms:
        print("ERROR ---> iatch does not have the correct dimension")
        print("ERROR ---> Expected: {} Given: {} ".format(natoms,iatch.shape[0]))
        return None

    cn = c_unit_bond_vectors(natoms, nchains, nbonds, maxnbondsperch,
                             &all_bonds[0,0], &coords[0,0], &iatch[0],
                             &uux[0,0], &uuy[0,0], &uuz[0,0])

    return cn

# ========================================================================================
def bond_bond_correlation(nchains,
                          maxnbondsperch,
                          np.ndarray[int, ndim=2, mode="c"] all_bonds,
                          np.ndarray[float, ndim=2, mode="c"] coords,
                          np.ndarray[int, ndim=1, mode="c"] iatch,
                          np.ndarray[float, ndim=1, mode="c"] cbb):

    nbonds = all_bonds.shape[0]
    natoms = coords.shape[0]

    if iatch.shape[0] != natoms:
        print("ERROR ---> iatch does not have the correct dimension")
        print("ERROR ---> Expected: {} Given: {} ".format(natoms,iatch.shape[0]))
        return None

    c_bond_bond_orientation(natoms, nchains, nbonds, maxnbondsperch,
                            &all_bonds[0,0], &coords[0,0], &iatch[0], &cbb[0])


# ========================================================================================
def setup_odf_intra(int maxnbondsperch):

    c_setup_odf_intra(maxnbondsperch)

# ========================================================================================
def odf_intra(iframe, nchains,
              np.ndarray[int, ndim=2, mode="c"] all_bonds,
              np.ndarray[float, ndim=2, mode="c"] coords,
              np.ndarray[int, ndim=1, mode="c"] iatch,
              np.ndarray[float, ndim=2, mode="c"] uux,
              np.ndarray[float, ndim=2, mode="c"] uuy,
              np.ndarray[float, ndim=2, mode="c"] uuz):

    """
        Calculate intra chain orientation correlation (for persistence length)

        Calculate orientational correlation functions
        for intra-chain vectors as function of chemical distance.
        It also calculates components and 4th moments.
        The input is supposed to contain a multiple of nvec vectors.
    """

    nbonds = all_bonds.shape[0]
    natoms = coords.shape[0]
    nchains = uux.shape[0]
    maxnbondsperch = uux.shape[1]

    if iatch.shape[0] != natoms:
        print("ERROR ---> iatch does not have the correct dimension")
        print("ERROR ---> Expected: {} Given: {} ".format(natoms,iatch.shape[0]))
        return None

    cn = c_odf_intra(iframe, natoms, nchains, nbonds, maxnbondsperch,
                             &all_bonds[0,0], &coords[0,0], &iatch[0],
                             &uux[0,0], &uuy[0,0], &uuz[0,0])

    return cn

# ========================================================================================
def avg_write_odf_intra(nframes, maxnbondsperch, filename):

    ftmp = filename.encode('utf-8')
    cdef char* fname = ftmp
    c_avg_write_odf_intra(nframes, maxnbondsperch, fname)