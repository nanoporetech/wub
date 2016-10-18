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

    def plot_arrays(self, data_map, title="", xlab="", ylab="", marker='.', legend_loc='upper right', legend=True):
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

    def plot_bars_simple(self, data_map, title="", xlab="", ylab="", alpha=0.6):
        """Plot simple bar chart from input dictionary.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param alpha: Alpha value.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        labels = data_map.keys()
        positions = np.arange(labels)
        plt.bar(positions, data_map.values(), align='center', alpha=alpha)
        plt.xticks(positions, labels)

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

    def plot_dicts(self, data_map, title="", xlab="", ylab="", marker='-', legend_loc='upper right'):
        """Plot elements of multiple dictionaries on a single plot.

        :param self: object.
        :param data_map: A dictionary with labels as keys and dictionaries as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param marker: Marker passed to the plot function.
        :param legend_loc: Location of legend.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, d in data_map.iteritems():
            plt.plot(d.keys(), d.values(), marker, label=label)

        plt.legend(loc=legend_loc)
        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_histograms(self, data_map, title="", xlab="", ylab="", bins=50, alpha=0.7, legend_loc='upper right', legend=True):
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
