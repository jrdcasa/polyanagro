import polyanagro as pag
import datetime
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from collections import defaultdict


class Map_Density(pag.Calculations):

    __slots__ = ['_phi_name', '_psi_name', '_pairs', '_dihedral_phi_2DArray', '_dihedral_psi_2DArray',
                 '_dihvalues1DArray', '_dihlabels1DArray', '_dihedralgrid', '_dyadsDict']

    # #######################################################################
    def __init__(self, trj, pairdihlist, dt=1, stride=1, log=None):

        """
        Constructor of the Map_Density object

        Args:
            trj:
            pairdihlist:
            dt:
            stride:
            log:
        """
        super().__init__(trj, dt=dt, stride=stride, logger=log)

        self._phi_name = pairdihlist[0]
        self._psi_name = pairdihlist[1]

        self._dihvalues1DArray = None
        self._dihlabels1DArray = None
        self._dihedral_phi_2DArray = None
        self._dihedral_psi_2DArray = None
        self._dihedralgrid = []
        self._dyadsDict = defaultdict()


    # #######################################################################
    def _generatepairs(self, ndx_filename):

        with open(ndx_filename, 'r') as fndx:
            contents = fndx.readlines()
            fndx.seek(0)
            for num, line in enumerate(fndx, 1):
                if re.match(r'.*' + self._phi_name + ' ]$', line):
                    idx_phi = num
            fndx.seek(0)
            for num, line in enumerate(fndx, 1):
                if re.match(r'.*' + self._psi_name + ' ]$', line):
                    idx_psi = num
        # Extract indices of phi angles
        dihedral_phi_2DArray_tmp = []
        if idx_phi != -1:
            while True:
                try:
                    tmp = [int(i) - 1 for i in contents[idx_phi].split()]
                    dihedral_phi_2DArray_tmp.append(tmp)
                    idx_phi += 1
                except ValueError:
                    break
                except IndexError:
                    m = "\t\t{} not label in {} index file".format(self._phi_name, ndx_filename)
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()
        # Extract indices of psi angles
        dihedral_psi_2DArray_tmp = []
        if idx_psi != -1:
            while True:
                try:
                    tmp = [int(i) - 1 for i in contents[idx_psi].split()]
                    dihedral_psi_2DArray_tmp.append(tmp)
                    idx_psi += 1
                except ValueError:
                    break
                except IndexError:
                    m = "\t\t{} not label in {} index file".format(self._psi_name, ndx_filename)
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

        # Order the arrays
        dihedral_phi_2DArray_tmp = np.array(dihedral_phi_2DArray_tmp, dtype=np.int32)
        dihedral_phi_2DArray_tmp = dihedral_phi_2DArray_tmp[dihedral_phi_2DArray_tmp[:, 0].argsort()]
        dihedral_psi_2DArray_tmp = np.array(dihedral_psi_2DArray_tmp, dtype=np.int32)
        dihedral_psi_2DArray_tmp = dihedral_psi_2DArray_tmp[dihedral_psi_2DArray_tmp[:, 0].argsort()]

        self._dihedral_phi_2DArray = []
        self._dihedral_psi_2DArray = []
        for iphix, iphi in enumerate(dihedral_phi_2DArray_tmp):
            try:
                idx_right_dihedral = np.where((dihedral_psi_2DArray_tmp[:, 0] == iphi[1]) &
                                              (dihedral_psi_2DArray_tmp[:, 1] == iphi[2]) &
                                              (dihedral_psi_2DArray_tmp[:, 2] == iphi[3]))[0]
                if len(idx_right_dihedral) == 0:
                    continue
                self._dihedral_phi_2DArray.append(dihedral_phi_2DArray_tmp[iphix,:])
                self._dihedral_psi_2DArray.append(dihedral_psi_2DArray_tmp[int(idx_right_dihedral), :])
            except TypeError:
                pass

        self._dihedral_phi_2DArray = np.array(self._dihedral_phi_2DArray, dtype=np.int32)
        self._dihedral_psi_2DArray = np.array(self._dihedral_psi_2DArray, dtype=np.int32)

    # #######################################################################
    def calculate(self, begin=0, unwrap_pbc=True, ndx_filename=None, half_dihedral=False):

        """
         Calculate a 2D Map of distributions

        Returns:

        """

        nframes = self._trajectory.get_numframes()
        if self._stride == 1:
            nframes_analysed = nframes
        else:
            nframes_analysed = int(((nframes - begin) / self._stride)) + 1
        m = "\t Num of frames to analyse: {}".format(nframes_analysed)
        print(m) if self._logger is None else self._logger.info(m)
        natoms = self._trajectory.natoms

        # Start calculations for each frame
        #nframes = self._trajectory.nframes
        ini = begin
        m = "\tUnwrap PBC coordinates: {}".format(unwrap_pbc)
        print(m) if self._logger is None else self._logger.info(m)

        # Main loop of frames ====================================================
        s = datetime.datetime.now()
        idx_f = 0
        X1 = np.zeros([natoms], float)
        Y1 = np.zeros([natoms], float)
        Z1 = np.zeros([natoms], float)

        self._generatepairs(ndx_filename)

        for iframe in range(ini, nframes, self._stride):

            # Estimated time (Use the 100 first frames to estimate the time)
            if idx_f == 0:
                f = datetime.datetime.now()
            if idx_f == 100:
                elapsed_time = datetime.datetime.now() - f
                k = int(((nframes - ini)/self._stride))/100
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

            # If pbc is false it is assumed that the trajectory is unwrapped
            if unwrap_pbc:
                self._unwrap_coordinates(iframe)
            else:
                self._coords_unwrap = self._trajectory.universe.trajectory[iframe].positions

            X1[:] = self._coords_unwrap[:, 0]
            Y1[:] = self._coords_unwrap[:, 1]
            Z1[:] = self._coords_unwrap[:, 2]

            self._dihvalues1DArray = np.zeros(len(self._dihedral_phi_2DArray), dtype=np.float64)
            self._dihlabels1DArray = np.zeros(len(self._dihedral_phi_2DArray), dtype=np.int32)
            iserror3 = self.dihDistNeigh(self._dihedral_phi_2DArray, X1, Y1, Z1,
                                         self._dihvalues1DArray, self._dihlabels1DArray)
            dihvaluespsi1DArray = np.zeros(len(self._dihedral_psi_2DArray), dtype=np.float64)
            dihlabelspsi1DArray = np.zeros(len(self._dihedral_psi_2DArray), dtype=np.int32)
            iserror3 = self.dihDistNeigh(self._dihedral_psi_2DArray, X1, Y1, Z1,
                                         dihvaluespsi1DArray, dihlabelspsi1DArray)

            isize = self._dihlabels1DArray.shape[0]
            for i in range(0, isize):
                self._dihedralgrid.append([self._dihvalues1DArray[i], dihvaluespsi1DArray[i]])

        self._writeDyadStatistics()
        # ================== Half dihedral ==================
        # Range [-180..180]

        if half_dihedral:
            xse = [-180.0, 180.0]
            yse = [-180.0, 180.0]
            ddist_range_tmp = list()
            for idih in self._dihedralgrid:
                l_tmp = []
                for j in range(0, 2):
                    if idih[j] >= 180.:
                        l_tmp.append(idih[j] - 360.)
                    else:
                        l_tmp.append(idih[j])
                ddist_range_tmp.append(l_tmp)
            self._dihedralgrid = ddist_range_tmp
            stat, density, x_edge, y_edge, xrandom_elements, yrandom_elements = \
                self._write2DHist(xse=xse, yse=yse)
            self._writeMultiHist(stat, x_edge, y_edge, xse=xse, yse=yse)
        else:
            xse = [0.0, 360.0]
            yse = [0.0, 360.0]
            stat, density, x_edge, y_edge, xrandom_elements, yrandom_elements = \
                self._write2DHist(xse=xse, yse=yse)
            self._writeMultiHist(stat, x_edge, y_edge, xse=xse, yse=yse)

        self._writeGNUPlot(stat, density, x_edge, y_edge, xrandom_elements, yrandom_elements, xse=xse, yse=yse)

    ########################################################################
    @staticmethod
    def dihDistNeigh(dih2DArray, X1, Y1, Z1, dihvalues1DArray, dihlabels1DAray):

        iserror = 1
        if len(dih2DArray) != 0:
            iserror = pag.dihDistCNeigh(dih2DArray, X1, Y1, Z1, dihvalues1DArray, dihlabels1DAray)
        return iserror

    ########################################################################
    def _write2DHist(self, nxplots=2, nyplots=3, xse=[0.0, 360.0], yse=[0.0, 360.0]):

        starttime = datetime.datetime.now()

        # Real data from trajectory
        x = [i[0] for i in self._dihedralgrid]
        y = [i[1] for i in self._dihedralgrid]
        # This is a dirty trick to include the limits in the map
        x.append(xse[0])
        x.append(xse[1])
        y.append(yse[0])
        y.append(yse[1])

        fig, axs = plt.subplots(nxplots, nyplots, figsize=(26, 16))

        # Create the histogram
        stat, x_edge, y_edge, binnumber = binned_statistic_2d(x, y, values=None, statistic='count', bins=(360, 360))

        # Calculate the bin centers
        x_center = (x_edge[:-1] + x_edge[1:]) / 2
        y_center = (y_edge[:-1] + y_edge[1:]) / 2
        # Calculate the total number of points
        total_points = len(x)

        # Normalize the counts to get the density
        density = stat / total_points

        # Add labels and titles
        xticks_labels = np.linspace(xse[0], xse[1], 7)
        yticks_labels = np.linspace(yse[0], yse[1], 7)
        for i in range(0, nxplots):
            for j in range(0, nyplots):
                axs[i, j].set_xlabel(r'$\phi$ (º)')
                axs[i, j].set_ylabel(r'$\psi$ (º)')
                axs[i, j].set_xlim(xse[0], xse[1])
                axs[i, j].set_xticks(xticks_labels)
                axs[i, j].set_yticks(yticks_labels)
                axs[i, j].set_ylim(yse[0], yse[1])

        # Plots
        npoints = 15000
        step_size = len(x) // npoints
        npoints_real = len(list(range(0, len(x), step_size)))
        if len(x) <= npoints_real:
            c0 = axs[0,0].scatter(x, y, 2, marker='o', color='black')
            axs[0, 0].set_title('All points ({})'.format(npoints_real))
            xrandom_elements = x
            yrandom_elements = y
        else:
            # xrandom_elements = random.sample(x, 15000)
            # random_indices = [x.index(element) for element in xrandom_elements]
            # yrandom_elements = [y[idx] for idx in random_indices]
            indices = list(range(0, len(x), step_size))
            xrandom_elements = []
            yrandom_elements = []
            for i in indices:
                xrandom_elements.append(x[i])
                yrandom_elements.append(y[i])
            c0 = axs[0,0].scatter(xrandom_elements, yrandom_elements, 2, marker='o', color='black')
            axs[0, 0].set_title('Only {} points randomly sampled'.format(npoints_real))
        c1 = axs[0,1].pcolormesh(x_center, y_center, stat.T, cmap='twilight', shading='auto')
        c2 = axs[0,2].pcolormesh(x_center, y_center, density.T, cmap='twilight', shading='auto')
        levels = np.linspace(0, np.max(stat.T), 10)
        c3 = axs[1,0].contour(x_center, y_center, stat.T, cmap='twilight', linewidths=1, levels=levels)
        # levels = np.linspace(0, np.max(density.T), 4)
        # c3 = axs[1,0].contour(x_center, y_center, density.T, cmap='twilight', linewidths=1, levels=levels)
        axs[1,0].clabel(c3, c3.levels, inline=True, fontsize=6)

        axs[0, 1].set_title('Frequency Histogram 2D')
        axs[0, 2].set_title('Frequency Histogram 2D')
        axs[1, 0].set_title('Contour plot')
        axs[1,1].set_aspect('equal')

        # Add color bar for each subplot
        cb1 = plt.colorbar(c1, ax=axs[0,1])
        cb1.set_label('Count')
        cb2 = plt.colorbar(c2, ax=axs[0,2])
        cb2.set_label('Density')
        cb3 = plt.colorbar(c3, ax=axs[1,0])
        cb3.set_label('Count')

        # Add 3D bar graph
        hist, xedges, yedges, binnumber = binned_statistic_2d(x, y, values=None, statistic='count', bins=(45, 45))
        # Calculate the total number of points
        total_points = len(x)
        # Normalize the counts to get the density
        density_b = hist / total_points

        # Construct arrays for the anchor positions of the bars.
        xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25, indexing="ij")
        xpos = xpos.ravel()
        ypos = ypos.ravel()
        zpos = 0

        # Construct arrays with the dimensions for the bars.
        dx = dy = 0.5 * np.ones_like(zpos)
        dz = density_b.ravel()
        axs[1,1].remove()
        axs[1,1] = fig.add_subplot(235, projection="3d")
        axs[1, 1].set_xlabel(r'$\phi$ (º)')
        axs[1, 1].set_ylabel(r'$\psi$ (º)')
        axs[1, 1].set_ylabel(r'Density')
        axs[1, 1].set_xlim(xse[0], xse[1])
        axs[1, 1].set_xticks(xticks_labels)
        axs[1, 1].set_yticks(xticks_labels)
        axs[1, 1].set_ylim(yse[0], yse[1])

        axs[1,1].bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
        axs[1,2].remove()

        # Save the graph
        filename = '{}_vs_{}.png'.format(self._phi_name, self._psi_name)
        plt.savefig(filename)

        endtime = datetime.datetime.now()
        m = "\tHistograms written in {}".format(filename)
        print(m) if self._logger is None else self._logger.info(m)
        m = "\t\tTime: {0:.2f} seconds".format((endtime-starttime).total_seconds())
        print(m) if self._logger is None else self._logger.info(m)

        return stat, density, x_edge, y_edge, xrandom_elements, yrandom_elements

    ########################################################################
    def _writeMultiHist(self, stat, x_edge, y_edge, xse=[0.0, 360.0], yse=[0.0, 360.0]):

        starttime = datetime.datetime.now()

        # Real data from trajectory
        x = [i[0] for i in self._dihedralgrid]
        y = [i[1] for i in self._dihedralgrid]
        # This is a dirty trick to include the limits in the map
        x.append(xse[0])
        x.append(xse[1])
        y.append(yse[0])
        y.append(yse[1])

        # Start with a square Figure.
        fig = plt.figure(figsize=(8, 8))
        # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
        # the size of the marginal axes and the main axes in both directions.
        # Also adjust the subplot parameters for a square plot.
        heigths = [1, 4]
        widths = [4, 1]

        gs = fig.add_gridspec(2, 2, width_ratios=widths, height_ratios=heigths,
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.20, hspace=0.20)

        # Create the Axes.
        ax = fig.add_subplot(gs[1, 0], xlim=(xse[0], xse[1]), ylim=(yse[0], yse[1]),
                             xlabel=r'$\phi$ (º)', ylabel=r'$\psi$ (º)')
        ax.pcolormesh(x_edge, y_edge, stat.T, cmap='twilight', shading='auto')

        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax, ylabel=r'PDF')
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax, xlabel="PDF")

        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # now determine nice limits by hand:
        binwidth = 1.0
        xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
        lim = (int(xymax / binwidth) + 1) * binwidth

        bins = np.arange(-lim, lim + binwidth, binwidth)

        histx, binx = np.histogram(x, bins=bins, density=True)
        histy, biny = np.histogram(y, bins=bins, density=True)
        ax_histx.bar(binx[:-1], histx, width=(binx[1] - binx[0]), color="black")
        ax_histy.barh(biny[:-1], histy, height=(bins[1] - bins[0]), color="black")

#        ax_histx.hist(x, bins=bins, , color="black", )
#        ax_histy.hist(y, bins=bins, orientation='horizontal', density=True, color="black")

        filename = '{}_vs_{}_b.png'.format(self._phi_name, self._psi_name)
        plt.savefig(filename)

        endtime = datetime.datetime.now()
        m = "\tHistograms written in {}".format(filename)
        print(m) if self._logger is None else self._logger.info(m)
        m = "\t\tTime: {0:.2f} seconds".format((endtime-starttime).total_seconds())
        print(m) if self._logger is None else self._logger.info(m)

    ########################################################################
    def _writeDyadStatistics(self):

        label_dyads = ["ug", "tg", "ut", "uu", "gg", "gt", "tt", "tu", "gu"]

        # Classify the angles
        for item in self._dihedralgrid:

            label = ""
            for iangles in item:
                if 0 < iangles < 120:
                    label += "g"
                elif 120 < iangles < 240:
                    label += "t"
                else:
                    label += "u"
            if label in self._dyadsDict.keys():
                self._dyadsDict[label] += 1
            else:
                self._dyadsDict[label] = 1

        for ilabel in label_dyads:
            if ilabel not in self._dyadsDict.keys():
                self._dyadsDict[ilabel] = 0

        # Write the data
        totaldiah = sum([item for key, item in self._dyadsDict.items()])
        filename = '{}_vs_{}.dat'.format(self._phi_name, self._psi_name)
        with open(filename, "w") as f:
            line = "0 tt {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["tt"], float(self._dyadsDict["tt"] / float(totaldiah)),
                       float(self._dyadsDict["tt"] / float(totaldiah))*100)
            line += "1 gg {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["gg"], float(self._dyadsDict["gg"] / float(totaldiah)),
                       float(self._dyadsDict["gg"] / float(totaldiah))*100)
            line += "2 uu {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["uu"], float(self._dyadsDict["uu"] / float(totaldiah)),
                       float(self._dyadsDict["uu"] / float(totaldiah))*100)
            line += "3 tg {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["tg"], float(self._dyadsDict["tg"] / float(totaldiah)),
                       float(self._dyadsDict["tg"] / float(totaldiah))*100)
            line += "4 gt {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["gt"], float(self._dyadsDict["gt"] / float(totaldiah)),
                       float(self._dyadsDict["gt"] / float(totaldiah))*100)
            line += "5 tu {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["tu"], float(self._dyadsDict["tu"] / float(totaldiah)),
                       float(self._dyadsDict["tu"] / float(totaldiah))*100)
            line += "6 ut {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["ut"], float(self._dyadsDict["ut"] / float(totaldiah)),
                       float(self._dyadsDict["ut"] / float(totaldiah))*100)
            line += "7 gu {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["gu"], float(self._dyadsDict["gu"] / float(totaldiah)),
                       float(self._dyadsDict["gu"] / float(totaldiah))*100)
            line += "8 ug {0:<10d} {1:<7.5f} {2:>4.1f}\n". \
                format(self._dyadsDict["ug"], float(self._dyadsDict["ug"] / float(totaldiah)),
                       float(self._dyadsDict["ug"] / float(totaldiah))*100)
            f.writelines(line)

    ########################################################################
    def _writeGNUPlot(self, stat, density, xedge, yedge, xr, yr, xse=[0.0, 360.0], yse=[0.0, 360.0]):

        # Real data from trajectory
        x = [i[0] for i in self._dihedralgrid]
        y = [i[1] for i in self._dihedralgrid]

        # Write data of the random points generated =================
        filename_random = 'gnuplot_{}_vs_{}_random.dat'.format(self._phi_name, self._psi_name)
        with open(filename_random, 'w') as fname:
            for idx in range(0, len(xr)):
                line = "{0:1f} {1:1f}\n".format(xr[idx], yr[idx])
                fname.writelines(line)

        # Write data 2D histogram frequency map =====================
        filename = 'gnuplot_{}_vs_{}_xedges.dat'.format(self._phi_name, self._psi_name)
        with open(filename, "w") as fname:
            for idx in range(0, len(xedge)):
                line = "{0:1f}\n".format(xedge[idx])
                fname.writelines(line)

        filename = 'gnuplot_{}_vs_{}_yedges.dat'.format(self._phi_name, self._psi_name)
        with open(filename, "w") as fname:
            for idx in range(0, len(yedge)):
                line = "{0:1f}\n".format(yedge[idx])
                fname.writelines(line)

        filename_binned = 'gnuplot_{}_vs_{}_binned.dat'.format(self._phi_name, self._psi_name)
        np.savetxt(filename_binned, np.transpose(stat))

        filename_density = 'gnuplot_{}_vs_{}_density.dat'.format(self._phi_name, self._psi_name)
        np.savetxt(filename_density, np.transpose(density))

        # Data for contour ==================================
        xlen, ylen = np.shape(stat)
        filename_data = 'gnuplot_{}_vs_{}_data.dat'.format(self._phi_name, self._psi_name)
        with open(filename_data, "w") as ff:
            for ix in range(0, xlen):
                for iy in range(0, ylen):
                    if stat[ix,iy] == 0: continue
                    line = "{} {} {}\n".format(xedge[ix], yedge[iy], stat[ix, iy])
                    ff.writelines(line)
            ff.writelines("\n")

         # Write gnuplot template
        line = 'reset\n'
        line += 'set style line 1 lt 1 ps 0.4 lc rgb "black"  pt 6 lw 0.2\n'
        line += 'set style line 2 lt 1 ps 0.4 lc rgb "red"    pt 4 lw 2.0\n'
        line += 'set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0\n'
        line += 'set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0\n'
        line += 'set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0\n'
        line += 'set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0\n'
        line += '\n'
        line += '########################################################\n'
        line += 'set encoding utf8\n'
        line += 'set term qt 1 enhanced dashed size 600,600 font "Arial,14"\n'
        line += '\n'
        line += 'set xlabel "{/Symbol f} (º)"\n'
        line += 'set ylabel "{/Symbol y} (º)"\n'
        line += 'set format x "%.0f"\n'
        line += 'set format y "%.0f"\n'
        line += 'set xrange[{0:d}:{1:d}]\n'.format(int(xse[0]), int(xse[1]))
        line += 'set yrange[{0:d}:{1:d}]\n'.format(int(yse[0]), int(yse[1]))
        line += 'set xtics 60\n'
        line += 'set mxtics 6\n'
        line += 'set ytics 60\n'
        line += 'set mytics 6\n'
        line += 'set grid\n'
        line += '\n'
        line += 'set title "{}"\n'.format("Points randomly sampled")
        line += 'p "{}" u 1:2 w p ls 1 notitle\n'.format(filename_random)
        line += '\n'
        line += 'set term qt 2 enhanced dashed size 600,600 font "Arial,14"\n'
        line += 'reset\n'
        line += 'stats "gnuplot_phi_vs_psi_binned.dat" matrix name "STATS"\n'
        line += 'max_frequency=STATS_max\n'
        line += 'set view map\n'
        line += 'set palette rgbformulae 22,13,-31\n'
        line += 'set size ratio -1\n'
        line += 'set xlabel "{/Symbol f} (º)"\n'
        line += 'set ylabel "{/Symbol y} (º)"\n'
        line += 'set format x "%.0f"\n'
        line += 'set format y "%.0f"\n'
        line += 'set xrange[{0:d}:{1:d}]\n'.format(int(xse[0]), int(xse[1]))
        line += 'set yrange[{0:d}:{1:d}]\n'.format(int(yse[0]), int(yse[1]))
        line += 'set xtics 60\n'
        line += 'set mxtics 6\n'
        line += 'set ytics 60\n'
        line += 'set mytics 6\n'
        line += 'set cbrange[1:max_frequency]\n'
        line += 'set title "{}"\n'.format("Heat map plot")
        line += 'splot "{}" using 1:2:3 w p pt 5 ps 0.5 palette linewidth 0 notitle\n'.format(filename_data)
        line += '\n'
        line += 'set term qt 3 enhanced dashed size 600,600 font "Arial,14"\n'
        line += 'unset view\n'
        line += 'unset title \n'
        line += 'set hidden3d\n'
        line += 'set zlabel "Frequency" rotate by 90\n'
        line += 'set dgrid3d 50,50 exp\n'
        line += 'splot "{}" with lines notitle\n'.format(filename_data)
        line += 'pause -1\n'

        with open("phi_vs_psi.gnu", 'w') as fgnu:
            fgnu.writelines(line)
