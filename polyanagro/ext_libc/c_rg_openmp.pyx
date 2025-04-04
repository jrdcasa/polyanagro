
# import both numpy and the Cython declarations for numpy
import cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from copy import copy


# declare the interface to the C code
cdef extern from "calc_rg.c":
    void c_rg_chain(int* mols, float* mass, float* coords_unwrap,
                    int nchains, int maxatomsch, float* rgsq_ich_iframe)

    void c_rg_chain_massweigth(int* mols, float* mass, float* coords_unwrap,
                               int nchains, int maxatomsch, float* rgsq_ich_iframe)


#OPENMP_ENABLED = True if USED_OPENMP else False


# ========================================================================================
def calc_rg_openmp(np.ndarray[int, ndim=2, mode="c"] mols,
                   np.ndarray[float, ndim=1, mode="c"] mass,
                   np.ndarray[float, ndim=2, mode="c"] coords_unwrap,
                   np.ndarray[float, ndim=2, mode="c"] rgsq_ich_iframe):

    cdef int nchains = mols.shape[0]
    cdef int maxatomsch = mols.shape[1]

    c_rg_chain(&mols[0,0], &mass[0], &coords_unwrap[0,0], nchains, maxatomsch, &rgsq_ich_iframe[0,0])

    return None

# ========================================================================================
def calc_rg_openmp_massweigth(np.ndarray[int, ndim=2, mode="c"] mols,
                              np.ndarray[float, ndim=1, mode="c"] mass,
                              np.ndarray[float, ndim=2, mode="c"] coords_unwrap,
                              np.ndarray[float, ndim=2, mode="c"] rgsq_ich_iframe):

    cdef int nchains = mols.shape[0]
    cdef int maxatomsch = mols.shape[1]

    c_rg_chain_massweigth(&mols[0,0], &mass[0], &coords_unwrap[0,0], nchains, maxatomsch, &rgsq_ich_iframe[0,0])

    return None

