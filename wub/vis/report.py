#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Report:

    # Maybe it would be a good idea to convert these utility methods to object
    # oriented matplotlib style.

    def __init__(self, pdf):
        """Class for plotting utilities on the top of matplotlib. Plots are saved in the specified file through the PDF backend.

        :param self: object.
        :param pdf: Output pdf.
        :returns: The report object.
        :rtype: Report

        """
        self.pdf = pdf
        self.plt = plt
        self.pages = PdfPages(pdf)

    def _set_properties_and_close(self, fig, title, xlab, ylab):
        """Utility method to set title, axis labels and close the figure.

        :param self: object.
        :param fig: The current figure.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :retuns: None
        :rtype: object
        """
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.close(fig)

    def plot_arrays(self, data_map, title="", xlab="", ylab="", marker='.', legend_loc='best', legend=True):
        """Plot multiple pairs of data arrays.

        :param self: object.
        :param data_map: A dictionary with labels as keys and tupples of data arrays (x,y) as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param marker: Marker passed to the plot function.
        :param legend_loc: Location of legend.
        :param legend: Plot legend if True
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, data_arrays in data_map.iteritems():
            plt.plot(data_arrays[0], data_arrays[1], marker, label=label)

        if legend:
            plt.legend(loc=legend_loc)
        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_bars_simple(self, data_map, title="", xlab="", ylab="", alpha=0.6, xticks_rotation=0, auto_limit=False):
        """Plot simple bar chart from input dictionary.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param alpha: Alpha value.
        :param xticks_rotation: Rotation value for x tick labels.
        :param auto_limit: Set y axis limits automatically.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        labels = data_map.keys()
        data = data_map.values()
        positions = np.arange(len(labels))
        plt.bar(positions, data, align='center', alpha=alpha)
        plt.xticks(positions, labels, rotation=xticks_rotation)

        if auto_limit:
            low, high = min(data), max(data)
            plt.ylim([(low - 0.5 * (high - low)), (high + 0.5 * (high - low))])

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_heatmap(self, data_matrix, title="", xlab="", ylab="", colormap=plt.cm.jet):
        """Plot heatmap of data matrix.

        :param self: object.
        :param data_matrix: 2D array to be plotted.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param colormap: matplotlib color map.
        :retuns: None
        :rtype: object
        """
        """
        """
        fig = plt.figure()

        p = plt.contourf(data_matrix)
        plt.colorbar(p, orientation='vertical', cmap=colormap)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_pcolor(self, data, title="", xlab="", ylab="", xticks=None, yticks=None, invert_yaxis=False, colormap=plt.cm.Blues, tick_size=5, tick_rotation=90):
        """Plot square heatmap of data matrix.

        :param self: object.
        :param data: 2D array to be plotted.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param xticks: X axis tick labels..
        :param yticks: Y axis tick labels..
        :param invert_yaxis: Invert Y axis if true.
        :param colormap: matplotlib color map.
        :param tick_size: Font size on tick labels.
        :param tick_rotation: Rotation of tick labels.
        :retuns: None
        :rtype: object
        """
        """
        """

        fig, ax = plt.subplots()
        hm = plt.pcolor(data, cmap=colormap)
        if invert_yaxis:
            ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')

        ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
        ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)

        ax.set_xticklabels(xticks, minor=False, fontsize=tick_size, rotation=tick_rotation)
        ax.set_yticklabels(yticks, minor=False, fontsize=tick_size)
        plt.colorbar(hm)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_dicts(self, data_map, title="", xlab="", ylab="", marker='-', legend_loc='best', legend=True, hist_style=False, cmap=plt.cm.rainbow, alpha=0.6):
        """Plot elements of multiple dictionaries on a single plot.

        :param self: object.
        :param data_map: A dictionary with labels as keys and dictionaries as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param marker: Marker passed to the plot function.
        :param legend_loc: Location of legend.
        :param legend: Hide legend if False.
        :param hist_style: Plot histogram-style bar plots.
        :param cmap: Colormap for histogram plots.
        :param alpha: Transparency value for histograms.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        if not hist_style:
            for label, d in data_map.iteritems():
                x, y = d.keys(), d.values()
                plt.plot(x, y, marker, label=label)
        else:
            color = iter(cmap(np.linspace(0, 1, len(data_map))))
            for label, d in data_map.iteritems():
                x, y = d.keys(), d.values()
                plt.bar(x, y, label=label,
                        align='center', color=next(color), alpha=alpha)

        if legend:
            plt.legend(loc=legend_loc)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_histograms(self, data_map, title="", xlab="", ylab="", bins=50, alpha=0.7, legend_loc='best', legend=True):
        """Plot histograms of multiple data arrays.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data arrays as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param bins: Number of bins.
        :param alpha: Transparency value for histograms.
        :param legend_loc: Location of legend.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, data in data_map.iteritems():
            if len(data) > 0:
                plt.hist(data, bins=bins, label=label, alpha=alpha)
        if legend:
            plt.legend(loc=legend_loc)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def close(self):
        """Close PDF backend. Do not forget to call this at the end of your script or your output will be damaged!

        :param self: object
        :returns: None
        :rtype: object
        """
        self.pages.close()

    def plot_line(self, data, x, y, title="",  xlab="", ylab=""):
        '''
        Generate a line plot from pandas dataframe

        :param data: pandas dataframe
        :param x: X axis data
        :param y: Y axis data
        :param title: Figure title
        :param xlab: X axis label
        :param ylab: Y axis label
        :return: None
        :rtype: object
        '''
        fig = plt.figure()

        plt.plot(data[x], data[y], 'k-', linewidth=1.5)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_scatter(self, data, x, y, title="",  xlab="", ylab="", alpha=0.5, ylim=None, xlim=None):
        '''
        Generates a scatter plot from a pandas dataframe

        :param data: Pandas dataframe
        :param x: X axis data
        :param y: Y axis data
        :param title: Figure title
        :param xlab: X axis label
        :param ylab: Y axis label
        :param alpha: opacity of data pionts
        :param ylim: Y axis limit
        :param xlim: X axis limit
        :return: None
        :rtype: object
        '''

        fig = plt.figure()
        plt.scatter(data[x], data[y], alpha=alpha)

        plt.ylim(ylim)
        plt.xlim(xlim)
        self._set_properties_and_close(fig, title, xlab, ylab)
