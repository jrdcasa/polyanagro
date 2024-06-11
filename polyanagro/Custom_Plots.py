import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


class Custom_Plots(Figure):

    # =========================================================================
    def __init__(self, figtitle=None):

        super(Custom_Plots, self).__init__()

        self._fig = None
        self._ax = None

    # =========================================================================
    def simple_line_plot(self, xdata, ydata, title=None, xlabel="x", ylabel="y",
                         xunits=None, yunits=None,
                         begin=0, errdata=None, path_to_save="./"):

        self._fig, self._ax = plt.subplots(1)

        if xunits is not None:
            xlabel += " ({})".format(xunits)
        if xunits is not None:
            ylabel += " ({})".format(yunits)

        if not os.path.isdir(path_to_save):
            m = "{} directory does not exist. Please create the directory first"
            return m

        if errdata is None:
            self._ax.set_title(title)
            self._ax.set_xlabel(xlabel)
            self._ax.set_ylabel(ylabel)
            ymin = min(ydata[begin:])
            ymin -= np.abs(ymin*0.001)
            ymax = max(ydata[begin:])
            ymax += np.abs(ymax*0.001)
            if ymin == ymax:
                return None
            self._ax.set_ylim((ymin, ymax))
            xdata_np = xdata.to_numpy()
            ydata_np = ydata.to_numpy()
            #idx_begin = np.int(np.where(xdata_np == begin)[0])
            idx_begin = int(np.where(xdata_np == begin)[0])
            self._ax.plot(xdata_np[idx_begin:], ydata_np[idx_begin:])
            filename = os.path.join(path_to_save, '{}.png'.format(ylabel.split()[0]))
            self._fig.savefig(filename, dpi=self._fig.dpi)
            plt.close()

        return None

    # =========================================================================
    def group_line_plot(self, nsubplots, data_df, xunits=None, yunits=None,
                        begin=0, path_to_save="./", fontsize=12,
                        fig_size=(15, 10), outputname='energy.jpg'):

        if not os.path.isdir(path_to_save):
            m = "{} directory does not exist. Please create the directory first"
            return m

        if isinstance(nsubplots, list):
            nx_subplots = nsubplots[0]
            ny_subplots = nsubplots[1]
        else:
            nx_subplots = 1
            ny_subplots = 1

        self._fig, self._ax = plt.subplots(nx_subplots, ny_subplots, figsize=fig_size)

        column_names = data_df.columns.tolist()

        if xunits is not None:
            xlabel = column_names[0] + " ({})".format(xunits)
        else:
            xlabel = column_names[0]

        for irow in range(nx_subplots):
            for icol in range(ny_subplots):
                igroup = icol + ny_subplots*irow + 1
                if yunits is not None:
                    ylabel = column_names[igroup] + " ({})".format(yunits[igroup-1])
                else:
                    ylabel = column_names[igroup]
                ymin = min(data_df[column_names[igroup]][begin:])
                ymin -= np.abs(ymin * 0.001)
                ymax = max(data_df[column_names[igroup]][begin:])
                ymax += np.abs(ymax * 0.001)
                if ymin == ymax:
                    continue
                try:
                    self._ax[irow, icol].set_ylim((ymin, ymax))
                    self._ax[irow, icol].set_title(column_names[igroup], fontsize=fontsize)
                    self._ax[irow, icol].set_xlabel(xlabel, fontsize=fontsize)
                    self._ax[irow, icol].set_ylabel(ylabel, fontsize=fontsize)
                    self._ax[irow, icol].tick_params(axis='both', which='major', labelsize=fontsize-4)
                    self._ax[irow, icol].plot(np.array(data_df[column_names[0]][begin:]),
                                              np.array(data_df[column_names[igroup]][begin:]), 'k-')
                    self._ax[irow, icol].grid()
                except IndexError:
                    self._ax[irow].set_ylim((ymin, ymax))
                    self._ax[irow].set_title(column_names[igroup], fontsize=fontsize)
                    self._ax[irow].set_xlabel(xlabel, fontsize=fontsize)
                    self._ax[irow].set_ylabel(ylabel, fontsize=fontsize)
                    self._ax[irow].tick_params(axis='both', which='major', labelsize=fontsize-4)
                    self._ax[irow].plot(np.array(data_df[column_names[0]][begin:]),
                                        np.array(data_df[column_names[igroup]][begin:]), 'k-')
                    self._ax[irow].grid()

        self._fig.tight_layout()

        base = os.path.splitext(outputname)[0]
        filename_jpg = os.path.join(path_to_save, '{}'.format(outputname))
        filename_csv = os.path.join(path_to_save, '{}.csv'.format(base))
        self._fig.savefig(filename_jpg)
        data_df.to_csv(filename_csv, sep=" ")


    # =========================================================================
    def simple_line_plot_improve(self, xdata, ydata, filename, type, title=None, xlabel="x", ylabel="y",
                                 xunits=None, yunits=None,
                                 begin=0, errdata=None, path_to_save="./", xmaxmin = None,
                                 ymaxmin = None):

        if xmaxmin is None:
            xmax = max(xdata[begin:])
            xmin = min(xdata[begin:])
            delta = (xmax - xmin) * 0.05
            xmin -= np.abs(xmin*delta)
            xmax += np.abs(xmax*delta)
            if xmin == xmax:
                return None
        elif not isinstance(xmaxmin, list):
            print ("LIST")
            exit()

        if ymaxmin is None:
            ymax = max(ydata[begin:])
            ymin = min(ydata[begin:])
            delta = (ymax - ymin) * 0.05
            ymin -= np.abs(ymin*delta)
            ymax += np.abs(ymax*delta)
            if ymin == ymax:
                return None
        elif not isinstance(xmaxmin, list):
            print ("LIST")
            exit()

        self._fig, self._ax = plt.subplots(1)

        if xunits is not None:
            xlabel += " ({})".format(xunits)
        if xunits is not None:
            ylabel += " ({})".format(yunits)

        if not os.path.isdir(path_to_save):
            m = "{} directory does not exist. Please create the directory first"
            return m

        if errdata is None:
            self._ax.set_title(title)
            self._ax.set_xlabel(xlabel)
            self._ax.set_ylabel(ylabel)

            self._ax.set_ylim((ymin, ymax))
            self._ax.plot(xdata[begin:], ydata[begin:])
            filename = os.path.join(path_to_save, '{}'.format(filename))
            plt.grid()
            self._fig.savefig(filename, dpi=self._fig.dpi)
            plt.close()

        return None




