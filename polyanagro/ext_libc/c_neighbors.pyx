# import both numpy and the Cython declarations for numpy
import cython
cimport cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from copy import copy
cdef bint USED_OPENMP

# cdef extern from "calc_internal_distances.c":
#     void c_setup_internal_distances(int natoms)
#     void c_internal_distances_iframe(int nchains, int* head_array, int* ichatbb, float* coords)
#     int c_internal_distances_return(float* rsq_intdist, int* rsqcount_intdist, float* rsq_avgdist)
#
# OPENMP_ENABLED = True if USED_OPENMP else False
#
# # ========================================================================================
# def setup_internal_distances(int maxnbondsperch):
#
#     c_setup_internal_distances(maxnbondsperch)
#
# # ========================================================================================
# def internal_distances_iframe(nchains,
#                               np.ndarray[int, ndim=1, mode="c"] head_array,
#                               np.ndarray[int, ndim=1, mode="c"] ichatbb,
#                               np.ndarray[float, ndim=2, mode="c"] coords):
#
#     c_internal_distances_iframe(nchains, &head_array[0], &ichatbb[0], &coords[0,0])
#
# # ========================================================================================
# def internal_distances_return(np.ndarray[float, ndim=1, mode="c"] rsq_intdist,
#                               np.ndarray[int, ndim=1, mode="c"] rsqcount_intdist,
#                               np.ndarray[float, ndim=1, mode="c"] rsq_avgdist):
#
#     maxnbondsperchain = c_internal_distances_return(&rsq_intdist[0],
#                                                     &rsqcount_intdist[0],
#                                                     &rsq_avgdist[0])
#
#     return maxnbondsperchain
#
#
# # ========================================================================================
# @cython.boundscheck(False)
# @cython.wraparound(False)
# def insternalchaindist_cython(int nchains, float[:,:] coords_unwrap, int[:] head_array, int[:] ichatbb,
#                               double[:] rsq_intdist,
#                               int[:] rsqcount_intdist):
#
#
#         cdef int ich, iatom_o, iatom_l, ilength
#         cdef double r0, r1, r2
#
#         #with cython.nogil:
#         for ich in range(0, nchains):
#             iatom_o = head_array[ich]
#             iatom_l = ichatbb[iatom_o]
#             while iatom_l>=0:
#                 ilength = 1
#                 while iatom_l >= 0:
#                     r0 = coords_unwrap[iatom_l, 0] - coords_unwrap[iatom_o, 0]
#                     r1 = coords_unwrap[iatom_l, 1] - coords_unwrap[iatom_o, 1]
#                     r2 = coords_unwrap[iatom_l, 2] - coords_unwrap[iatom_o, 2]
#                     rsq_intdist[ilength] += (r0*r0 + r1*r1 + r2*r2)
#                     rsqcount_intdist[ilength] += 1
#                     # with open("ord.txt",'a') as f:
#                     #     line = "{0:d}  {1:d}   {2:d}\n".format(iatom_o, iatom_l, ilength)
#                     #     f.write(line)
#                     # if ilength == 1:
#                     #     with open("bb.txt",'a') as f:
#                     #         line = "{0:d} {1:f} {2:f} {3:d} {4:d} {5:20.10f} {6:d}\n".format(ich, r0*r0 + r1*r1 + r2*r2, np.sqrt(r0*r0 + r1*r1 + r2*r2), iatom_l+1, iatom_o+1, rsq_intdist[ilength],  rsqcount_intdist[ilength])
#                     #         f.write(line)
#                     iatom_l = ichatbb[iatom_l]
#                     ilength += 1
#                 iatom_o = ichatbb[iatom_o]
#                 iatom_l = ichatbb[iatom_o]
#
