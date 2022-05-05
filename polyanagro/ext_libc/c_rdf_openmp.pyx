import cython
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

# declare the interface to the C code
cdef extern from "calc_rdf.c":
    void c_setup_rdf(int, float* )
    void c_rdf_hist(int, int, int, int, float, int*, int*, float*, float*, int*)
    void c_rdf_gr()


# ========================================================================================
def setup_rdf_init(nbins,
                   np.ndarray[float, ndim=1, mode="c"] total_rdf):

    """

    :param nbins: Number of bins in the histogram
    :param total_rdf: Total RDF historgram (dim: nbins)

    """

    cdef int c_nbins = nbins

    c_setup_rdf(c_nbins, &total_rdf[0])

# ========================================================================================
def rdf_hist(nat_A, nat_B, nbins, delta_r,
             np.ndarray[int, ndim=1, mode="c"] atindex_A,
             np.ndarray[int, ndim=1, mode="c"] atindex_B,
             np.ndarray[float, ndim=2, mode="c"] wrapped_coords,
             np.ndarray[float, ndim=1, mode="c"] box,
             np.ndarray[int, ndim=1, mode="c"] hist_rdf):

    """

    :param nat_A: Number of atoms of the set A
    :param nat_B: Number of atoms of the set B
    :param nbins
    :param delta_r
    :param atindex_A: Atom index of the first set
    :param atindex_B: Atom index of the second set
    :param wrapped_coords: Wrapped coordinates of all atoms [nat,3] (in angstroms)
    :param box: Simulation box (in angstroms)
    :param total_rdf: Total RDF historgram (dim: nbins)

    :return:
    """

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B
    cdef int c_natoms = wrapped_coords.shape[0]

    c_rdf_hist(c_natoms, c_nat_A, c_nat_B, nbins, delta_r,
               &atindex_A[0], &atindex_B[0],
               &wrapped_coords[0,0], &box[0], &hist_rdf[0])

    return True

# ========================================================================================
def rdf_gr():

    print("Pass")