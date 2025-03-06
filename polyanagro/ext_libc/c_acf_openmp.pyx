import cython
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

# declare the interface to the C code
cdef extern from "calc_acf.c":
    cdef bint USED_OPENMP
    void c_acf_e2e (int ndim, int nchains, int ndump,
                    float* uree, float* rEE_acf)
    void c_acf2_e2e (int ndim, int nchains, int ndump,
                     float* ree, float* rEE_acf)

OPENMP_ENABLED = True if USED_OPENMP else False

# ========================================================================================
def calc_acf_ete(np.ndarray[np.float32_t, ndim=3, mode="c"] uree,
                 np.ndarray[np.float32_t, ndim=1, mode="c"] rEE_acf):

    # print("idump:0, ich:0",uree[:,0,0])
    # print("idump:0, ich:1",uree[:,1,0])
    # print("idump:0, ich:2",uree[:,2,0])
    #
    # print("idump:1, ich:0",uree[:,0,1])
    # print("idump:1, ich:1",uree[:,1,1])
    # print("idump:1, ich:2",uree[:,2,1])
    #
    # print("idump:2, ich:0",uree[:,0,2])
    # print("idump:2, ich:1",uree[:,1,2])
    # print("idump:2, ich:2",uree[:,2,2])

    cdef int ndim = uree.shape[0]
    cdef int nchains = uree.shape[1]
    cdef int ndumps = uree.shape[2]

    c_acf_e2e(ndim, nchains, ndumps, &uree[0,0,0], &rEE_acf[0])

    return None

# ========================================================================================
def calc_acf2_ete(np.ndarray[np.float32_t, ndim=3, mode="c"] ree,
                  np.ndarray[np.float32_t, ndim=1, mode="c"] rEE2_acf):

    # print("idump:0, ich:0",ree[:,0,0])
    # print("idump:0, ich:1",ree[:,1,0])
    # print("idump:0, ich:2",ree[:,2,0])
    #
    # print("idump:1, ich:0",ree[:,0,1])
    # print("idump:1, ich:1",ree[:,1,1])
    # print("idump:1, ich:2",ree[:,2,1])
    #
    # print("idump:2, ich:0",ree[:,0,2])
    # print("idump:2, ich:1",ree[:,1,2])
    # print("idump:2, ich:2",ree[:,2,2])

    cdef int ndim = ree.shape[0]
    cdef int nchains = ree.shape[1]
    cdef int ndumps = ree.shape[2]

    c_acf2_e2e(ndim, nchains, ndumps, &ree[0,0,0], &rEE2_acf[0])
    return None
