import cython
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

# declare the interface to the C code
cdef extern from "calc_rdf.c":
    cdef bint USED_OPENMP

    void c_setup_rdf(int, float* )
    void c_rdf_hist(int, int, int, int, float, float, int*, int*,
                    float*, float*, int*, int*, int*, int*, int*,
                    int*, int*, int*)
    void c_rdf_hist_openmp(int, int, int, int, float, float, int*, int*,
                    float*, float*, int*, int*, int*, int*, int*,
                    int*, int*, int*)
    void c_rdf_hist_excl_openmp(int, int, int, int, float, float, int*, int*, int*,
                                float*, float*, int*, int*, int*, int*, int*,
                                int*, int*, int*)
    void c_rdf_gr(int, int, int, int, float, float,
                  int*,
                  int*,
                  int*,
                  int*,
                  float*, int,
                  float*, int,
                  float*, int)
OPENMP_ENABLED = True if USED_OPENMP else False


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
def rdf_hist(nat_A, nat_B, nbins, delta_r, cutoff,
             np.ndarray[int, ndim=1, mode="c"] atindex_A,
             np.ndarray[int, ndim=1, mode="c"] atindex_B,
             np.ndarray[float, ndim=2, mode="c"] wrapped_coords,
             np.ndarray[float, ndim=1, mode="c"] box,
             np.ndarray[int, ndim=1, mode="c"] iatch,
             np.ndarray[int, ndim=1, mode="c"] hist_total_rdf,
             np.ndarray[int, ndim=1, mode="c"] hist_intra_rdf,
             np.ndarray[int, ndim=1, mode="c"] hist_inter_rdf,
             np.ndarray[int, ndim=1, mode="c"] hist_self_rdf):

    """

    :param nat_A: Number of atoms of the set A
    :param nat_B: Number of atoms of the set B
    :param nbins
    :param delta_r
    :param cutoff
    :param atindex_A: Atom index of the first set
    :param atindex_B: Atom index of the second set
    :param wrapped_coords: Wrapped coordinates of all atoms [nat,3] (in angstroms)
    :param box: Simulation box (in angstroms)
    :param iatch
    :param hist_total_rdf: Total RDF historgram (dim: nbins)
    :param hist_intra_rdf: Total RDF historgram (dim: nbins)
    :param hist_inter_rdf: Total RDF historgram (dim: nbins)
    :param hist_self_rdf: Total RDF historgram (dim: nbins)

    :return:
    """

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B
    cdef int c_natoms = wrapped_coords.shape[0]

    cdef int c_npairs_rdf_total[1]
    c_npairs_rdf_total[:] = [0]
    cdef int c_npairs_rdf_intra[1]
    c_npairs_rdf_intra[:] = [0]
    cdef int c_npairs_rdf_inter[1]
    c_npairs_rdf_inter[:] = [0]

    c_rdf_hist(c_natoms, c_nat_A, c_nat_B, nbins, delta_r, cutoff,
               &atindex_A[0], &atindex_B[0],
               &wrapped_coords[0,0], &box[0], &iatch[0],
               &hist_total_rdf[0], &hist_intra_rdf[0],
               &hist_inter_rdf[0], &hist_self_rdf[0],
               &c_npairs_rdf_total[0], &c_npairs_rdf_intra[0], &c_npairs_rdf_inter[0])

    return c_npairs_rdf_total[0], c_npairs_rdf_intra[0], c_npairs_rdf_inter[0]

# ========================================================================================
def rdf_hist_openmp(nat_A, nat_B, nbins, delta_r, cutoff,
                    np.ndarray[int, ndim=1, mode="c"] atindex_A,
                    np.ndarray[int, ndim=1, mode="c"] atindex_B,
                    np.ndarray[float, ndim=2, mode="c"] wrapped_coords,
                    np.ndarray[float, ndim=1, mode="c"] box,
                    np.ndarray[int, ndim=1, mode="c"] iatch,
                    np.ndarray[int, ndim=1, mode="c"] hist_total_rdf,
                    np.ndarray[int, ndim=1, mode="c"] hist_intra_rdf,
                    np.ndarray[int, ndim=1, mode="c"] hist_inter_rdf,
                    np.ndarray[int, ndim=1, mode="c"] hist_self_rdf):

    """

    :param nat_A: Number of atoms of the set A
    :param nat_B: Number of atoms of the set B
    :param nbins
    :param delta_r
    :param cutoff
    :param atindex_A: Atom index of the first set
    :param atindex_B: Atom index of the second set
    :param wrapped_coords: Wrapped coordinates of all atoms [nat,3] (in angstroms)
    :param box: Simulation box (in angstroms)
    :param iatch
    :param hist_total_rdf: Total RDF historgram (dim: nbins)
    :param hist_intra_rdf: Total RDF historgram (dim: nbins)
    :param hist_inter_rdf: Total RDF historgram (dim: nbins)
    :param hist_self_rdf: Total RDF historgram (dim: nbins)

    :return:
    """

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B
    cdef int c_natoms = wrapped_coords.shape[0]

    cdef int c_npairs_rdf_total[1]
    c_npairs_rdf_total[:] = [0]
    cdef int c_npairs_rdf_intra[1]
    c_npairs_rdf_intra[:] = [0]
    cdef int c_npairs_rdf_inter[1]
    c_npairs_rdf_inter[:] = [0]

    c_rdf_hist_openmp(c_natoms, c_nat_A, c_nat_B, nbins, delta_r, cutoff,
               &atindex_A[0], &atindex_B[0],
               &wrapped_coords[0,0], &box[0], &iatch[0],
               &hist_total_rdf[0], &hist_intra_rdf[0],
               &hist_inter_rdf[0], &hist_self_rdf[0],
               &c_npairs_rdf_total[0], &c_npairs_rdf_intra[0], &c_npairs_rdf_inter[0])

    return c_npairs_rdf_total[0], c_npairs_rdf_intra[0], c_npairs_rdf_inter[0]

# ========================================================================================
def rdf_hist_excl_openmp(nat_A, nat_B, nbins, delta_r, cutoff,
                         np.ndarray[int, ndim=1, mode="c"] atindex_A,
                         np.ndarray[int, ndim=1, mode="c"] atindex_B,
                         np.ndarray[int, ndim=2, mode="c"] excl_array,
                         np.ndarray[float, ndim=2, mode="c"] wrapped_coords,
                         np.ndarray[float, ndim=1, mode="c"] box,
                         np.ndarray[int, ndim=1, mode="c"] iatch,
                         np.ndarray[int, ndim=1, mode="c"] hist_total_rdf,
                         np.ndarray[int, ndim=1, mode="c"] hist_intra_rdf,
                         np.ndarray[int, ndim=1, mode="c"] hist_inter_rdf,
                         np.ndarray[int, ndim=1, mode="c"] hist_self_rdf):

    """

    :param nat_A: Number of atoms of the set A
    :param nat_B: Number of atoms of the set B
    :param nbins
    :param delta_r
    :param cutoff
    :param excl_array
    :param atindex_A: Atom index of the first set
    :param atindex_B: Atom index of the second set
    :param wrapped_coords: Wrapped coordinates of all atoms [nat,3] (in angstroms)
    :param box: Simulation box (in angstroms)
    :param iatch
    :param hist_total_rdf: Total RDF historgram (dim: nbins)
    :param hist_intra_rdf: Total RDF historgram (dim: nbins)
    :param hist_inter_rdf: Total RDF historgram (dim: nbins)
    :param hist_self_rdf: Total RDF historgram (dim: nbins)

    :return:
    """

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B
    cdef int c_natoms = wrapped_coords.shape[0]

    cdef int c_npairs_rdf_total[1]
    c_npairs_rdf_total[:] = [0]
    cdef int c_npairs_rdf_intra[1]
    c_npairs_rdf_intra[:] = [0]
    cdef int c_npairs_rdf_inter[1]
    c_npairs_rdf_inter[:] = [0]

    c_rdf_hist_excl_openmp(c_natoms, c_nat_A, c_nat_B, nbins, delta_r, cutoff,
                           &atindex_A[0], &atindex_B[0], &excl_array[0,0],
                           &wrapped_coords[0,0], &box[0], &iatch[0],
                           &hist_total_rdf[0], &hist_intra_rdf[0],
                           &hist_inter_rdf[0], &hist_self_rdf[0],
                           &c_npairs_rdf_total[0], &c_npairs_rdf_intra[0], &c_npairs_rdf_inter[0])

    return c_npairs_rdf_total[0], c_npairs_rdf_intra[0], c_npairs_rdf_inter[0]




# ========================================================================================
def rdf_gr(nframes, nat_A, nat_B, nbins, cutoff, volume_avg,
           np.ndarray[int, ndim=1, mode="c"] hist_total_rdf,
           np.ndarray[int, ndim=1, mode="c"] hist_intra_rdf,
           np.ndarray[int, ndim=1, mode="c"] hist_inter_rdf,
           np.ndarray[int, ndim=1, mode="c"] hist_self_rdf,
           np.ndarray[float, ndim=1, mode="c"] total_rdf, npairs_rdf_total,
           np.ndarray[float, ndim=1, mode="c"] intra_rdf, npairs_rdf_intra,
           np.ndarray[float, ndim=1, mode="c"] inter_rdf, npairs_rdf_inter):

    cdef int c_nat_A = nat_A
    cdef int c_nat_B = nat_B

    c_rdf_gr(nframes, c_nat_A, c_nat_B, nbins, cutoff, volume_avg,
             &hist_total_rdf[0], &hist_intra_rdf[0],
             &hist_inter_rdf[0], &hist_self_rdf[0],
             &total_rdf[0], npairs_rdf_total,
             &intra_rdf[0], npairs_rdf_intra,
             &inter_rdf[0], npairs_rdf_inter)
