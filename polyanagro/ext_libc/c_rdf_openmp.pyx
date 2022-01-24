import cython
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

# declare the interface to the C code
cdef extern from "calc_rdf.c":
    void c_setup_rdf(float, float )
    void c_calc_rdf_equal(int, long*, float*, float*)


# ========================================================================================
def setup_rdf_init(max_box, delta_r = 0.02):

    """

    :param max_box: Maximum box distance in angstroms
    :param delta_r: Bin width in angstroms (default = 0.02A)

    """

    cdef float c_max_box = max_box
    cdef float c_delta_r = delta_r

    c_setup_rdf(c_max_box, c_delta_r)

# ========================================================================================
def rdf_calc(nat_A, nat_B,
             np.ndarray[long  , ndim=1, mode="c"] atindex_A,
             np.ndarray[long  , ndim=1, mode="c"] atindex_B,
             np.ndarray[float, ndim=2, mode="c"] coords_A,
             np.ndarray[float, ndim=2, mode="c"] coords_B,
             np.ndarray[float, ndim=1, mode="c"] box):

    """

    :param nat_A: Number of atoms of the set A
    :param nat_B: Number of atoms of the set B
    :param coords_A: Wrapped coordinates of the set A [nat_A,3] (in angstroms)
    :param coords_B: Wrapped coordinates of the set B [nat_B,3] (in angstroms)
    :param box: Simulation box (in angstroms)

    :return:
    """

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B
    cdef int dim1 = coords_A.shape[0]
    cdef int dim2 = coords_B.shape[0]

    if dim1 != nat_A or dim2 !=nat_B:
        return False
    if dim1 != atindex_A.shape[0] or dim2 != atindex_B.shape[0]:
        return False

    if np.array_equal(atindex_A, atindex_B):
        c_calc_rdf_equal(c_nat_A, &atindex_A[0],
                         &coords_A[0,0], &box[0])
    else:
        print("J")
    #
    # c_calc_rdf_equal(c_nat_A, c_nat_B, &atindex_A[0], &atindex_B[0],
    #                  &coords_A[0,0], &coords_B[0,0], &box[0])

    return True