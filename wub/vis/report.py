#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Report:

    """ Message here """

    def __init__(self, pdf):
        """ Sample message here """
        self.pdf = pdf
        self.pages = PdfPages(pdf)

    def _set_properties_and_close(self, fig, title, xlab, ylab):
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.close(fig)

    def plot_arrays(self, a, b, title="", xlab="", ylab=""):
        """
        """
        fig = plt.figure()

        plt.plot(a, b, 'b.')

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_heatmap(self, z, title="", xlab="", ylab=""):
        """
        """
        fig = plt.figure()

        p = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_hashes(self, data_map, title="", xlab="", ylab="", plot_type='lines'):
        """
        """
        fig = plt.figure()

        for label, d in data_map.iteritems():
            plt.plot(d.keys(), d.values(), '-', label=label)

        plt.legend(loc='upper right')
        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_heatmap(self, z, title="", xlab="", ylab=""):
        """
        """
        fig = plt.figure()

        p = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_hists(self, data_map, title="", xlab="", ylab="", bins=50):
        """
        """
        fig = plt.figure()

        for label, data in data_map.iteritems():
            plt.hist(data, bins=bins, label=label, alpha=0.7)
        plt.legend(loc='upper right')

        self._set_properties_and_close(fig, title, xlab, ylab)

    def close(self):
        """
        """
        self.pages.close()
