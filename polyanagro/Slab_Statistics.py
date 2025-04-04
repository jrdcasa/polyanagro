import numpy as np
import datetime
import polyanagro as pag


# ===============================================================
class Slab_Statistics(pag.Calculations):

    # #########################################################################
    def __init__(self, trj, normal, dt=1, stride=1, log=None):

        """
        Initialize the instance to make calculations using slabs along the
        direction normal
        Args:
            trj:
            normal (int) : A number indicating the normal direction. 0: x, 1: y and 2: z
            dt:
            stride:
            log:
        """

        super().__init__(trj, dt=dt, stride=stride, logger=log)

        self._nslabs = None    # Number of slabs along the normal axis
        self._ds = None        # Thickness of the slab along the normal axis
        self._normal = normal  # Normal direction. 0:x, 1:y and 2:z


    # #########################################################################
    def normal_coord_max_min(self, begin=0):

        """
        Get the maximum and minimum coordinates in the normal direction, for the full
        trajectory

        Returns:

        """

        nchains = len(self._nmols_array)
        nframes = self._trajectory.get_numframes()
        if self._stride == 1:
            nframes_analysed = nframes
        else:
            nframes_analysed = int(((nframes - begin) / self._stride)) + 1
        m = "\t Num of frames to analyse: {}".format(nframes_analysed)
        print(m) if self._logger is None else self._logger.info(m)

        normal_coordinates_maxmin = np.zeros([nframes, 2], dtype=int)
        x_min_max = np.zeros([nframes, 2], dtype=np.float32)
        y_min_max = np.zeros([nframes, 2], dtype=np.float32)
        z_min_max = np.zeros([nframes, 2], dtype=np.float32)
        boxx = np.zeros([nframes], dtype=np.float32)
        boxy = np.zeros([nframes], dtype=np.float32)
        boxz = np.zeros([nframes], dtype=np.float32)

        # Main loop of frames ====================================================
        s = datetime.datetime.now()
        idx_f = 0
        for iframe in range(begin, nframes, self._stride):

            # Estimated time (Use the 100 first frames to estimate the time)
            if idx_f == 0:
                f = datetime.datetime.now()
            if idx_f == 100:
                elapsed_time = datetime.datetime.now() - f
                k = int(((nframes - begin)/self._stride))/100
                estimated_time = k*elapsed_time.total_seconds()
                m = "\tESTIMATED TIME using {0:d} frames ({1:s} seconds): {2:.2f} seconds".format \
                    (100, str(elapsed_time.total_seconds()), estimated_time)
                print(m) if self._logger is None else self._logger.info(m)

            # Write info
            if iframe%self._freq == 0:
                elapsed_time = datetime.datetime.now() - s
                m = "\tIFRAME: {1:d} of {0:d} in {2:s} seconds".format \
                    (nframes_analysed, iframe, str(elapsed_time.total_seconds()))
                print(m) if self._logger is None else self._logger.info(m)

            # Check if the trajectory is unwrapped or not
            self._coords_wrap = self._trajectory.universe.trajectory[iframe].positions
            boxx[iframe] = self._trajectory.universe.trajectory[iframe].dimensions[0]
            boxy[iframe] = self._trajectory.universe.trajectory[iframe].dimensions[1]
            boxz[iframe] = self._trajectory.universe.trajectory[iframe].dimensions[2]

            idxx1 = idxx2 = None
            idxy1 = idxy2 = None
            idxz1 = idxz2 = None

            if (self._coords_wrap[:,0] <0.0).any() < 0.0 or (self._coords_wrap[:,0] >boxx[iframe]).any():
                iswrapped = False
                idxx1 = np.where(self._coords_wrap[:,0]<0.0)
                idxx2 = np.where(self._coords_wrap[:,0]>boxx[iframe])
            elif (self._coords_wrap[:,1] <0.0).any() < 0.0 or (self._coords_wrap[:,1] >boxy[iframe]).any():
                iswrapped = False
                idxy1 = np.where(self._coords_wrap[:,1]<0.0)
                idxy2 = np.where(self._coords_wrap[:,1]>boxy[iframe])
            elif (self._coords_wrap[:,2] <0.0).any() < 0.0 or (self._coords_wrap[:,2] >boxz[iframe]).any():
                iswrapped = False
                idxz1 = np.where(self._coords_wrap[:,2]<0.0)
                idxz2 = np.where(self._coords_wrap[:,2]>boxz[iframe])
            else:
                iswrapped = True

            if  not iswrapped:
                m = "\t\t ERROR!!!. The trajectory should be wrapped.\n"
                m += "\t\t          Please wrap the trajectory before to use the program.\n"
                m += "\t\t          Frame {} of the trajectory".format(iframe)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            x_min_max[iframe] = [np.min(self._coords_wrap[:,0]), np.max(self._coords_wrap[:,0])]
            y_min_max[iframe] = [np.min(self._coords_wrap[:,1]), np.max(self._coords_wrap[:,1])]
            z_min_max[iframe] = [np.min(self._coords_wrap[:,2]), np.max(self._coords_wrap[:,2])]

            if self._normal == 0:
                normal_coordinates_maxmin[iframe,:] = [np.min(self._coords_wrap[:,0]), np.max(self._coords_wrap[:,0])]
            elif self._normal == 1:
                normal_coordinates_maxmin[iframe,:] = [np.min(self._coords_wrap[:,1]), np.max(self._coords_wrap[:,1])]
            elif self._normal == 2:
                normal_coordinates_maxmin[iframe,:] = [np.min(self._coords_wrap[:,2]), np.max(self._coords_wrap[:,2])]

        # Write averages
        x_min_max_avg = [np.mean(x_min_max[:,0]/10.0), np.mean(x_min_max[:,1]/10.0)]
        y_min_max_avg = [np.mean(y_min_max[:,0]/10.0), np.mean(y_min_max[:,1]/10.0)]
        z_min_max_avg = [np.mean(z_min_max[:,0]/10.0), np.mean(z_min_max[:,1]/10.0)]
        x_min_max_std = [np.std(x_min_max[:,0]/10.0), np.std(x_min_max[:,1]/10.0)]
        y_min_max_std = [np.std(y_min_max[:,0]/10.0), np.std(y_min_max[:,1]/10.0)]
        z_min_max_std = [np.std(z_min_max[:,0]/10.0), np.std(z_min_max[:,1]/10.0)]
        boxx_avg = np.mean(boxx)/10.0
        boxy_avg = np.mean(boxy)/10.0
        boxz_avg = np.mean(boxz)/10.0
        boxx_std = np.std(boxx)/10.0
        boxy_std = np.std(boxy)/10.0
        boxz_std = np.std(boxz)/10.0

        m = "\t ======== AVG MINIMUN and MAXIMUM atomistic coordinates ========\n"
        m += "\t\t Wrapped trajectory = {}\n".format(iswrapped)
        m += "\t\t Box dimensions = [{0:.2f} +- {1:.2f}," \
                                     "{2:.2f} +- {3:.2f}," \
                                     "{4:.2f} +- {5:.2f}] (nm)\n".format(boxx_avg, boxx_std,
                                                                        boxy_avg, boxy_std,
                                                                        boxz_avg, boxz_std )
        m += "\t\t x min = {0:.2f} +- {1:.2f} x max = {2:.2f} +- {2:.2f} (nm)\n".\
            format(x_min_max_avg[0], x_min_max_std[0], x_min_max_avg[1], x_min_max_std[1])
        m += "\t\t y min = {0:.2f} +- {1:.2f} y max = {2:.2f} +- {2:.2f} (nm)\n".\
            format(y_min_max_avg[0], y_min_max_std[0], y_min_max_avg[1], y_min_max_std[1])
        m += "\t\t z min = {0:.2f} +- {1:.2f} z max = {2:.2f} +- {2:.2f} (nm)\n".\
            format(z_min_max_avg[0], z_min_max_std[0], z_min_max_avg[1], z_min_max_std[1])

        if self._normal == 0:
            m += "\t\t Atomistic thickness (X normal) = {0:.2f} +- {1:.2f} (nm)\n".\
                format(x_min_max_avg[1] - x_min_max_avg[0], x_min_max_std[1] - x_min_max_std[0])
        elif self._normal == 1:
            m += "\t\t Atomistic thickness (Y normal) = {0:.2f} +- {1:.2f} (nm)\n".\
                format(y_min_max_avg[1] - y_min_max_avg[0], y_min_max_std[1] - y_min_max_std[0])
        elif self._normal == 2:
            m += "\t\t Atomistic thickness (Z normal) = {0:.2f} +- {1:.2f} (nm)\n".\
                format(z_min_max_avg[1] - z_min_max_avg[0], z_min_max_std[1] - z_min_max_std[0])
        m += "\t ======== AVG MINIMUN and MAXIMUM atomistic coordinates ========"

        print(m) if self._logger is None else self._logger.info(m)

        return iswrapped, normal_coordinates_maxmin
