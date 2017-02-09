#! /usr/bin/env python
# -*- coding: utf-8 -*-

from core import _V

import h5py

pc_list = ['BT_PCCOMP', 'PCSCORES', 'TOTAL_PCCOMP']
misc_list = ['CLOUD_COEF_FILENAME', 'COEF_FILENAME', 'ID_COMMON_NAME',
             'INSTRUMENT', 'NCHANNELS', 'PC_COEF_FILENAME', 'SATELLITE',
             'SCAT_COEF_FILENAME', 'WAVENUMBERS']


class PCMISC(dict, _V):

    def __init__(self):
        dict.__init__(self)
        for key in misc_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}


class PCCOMP(dict, _V):
    """
    The PCclass
    PC object is Read-Only and contains minimum information
    """

    def __init__(self):
        dict.__init__(self)
        for key in pc_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}


class pc():
    """
      combine PC and MISC
    """

    def __init__(self, fileName):
        f = h5py.File(fileName, 'r')
        h5 = f['/PCCOMP/']
        h5misc = f['/MISC/']
        self.pc = PCCOMP()
        self.misc = PCMISC()
        self.pc.loadh5(h5)
        self.misc.loadh5(h5misc)

    def display(self):
        print "PC-----------------------"
        self.pc.display()
        print "MISC---------------------"
        self.misc.display()


if __name__ == '__main__':

    import h5py
    # Open file ReadOnly
    fileName = '../rttov_tests/pc.h5'
    f = h5py.File(fileName, 'r')

    # get the Dataset
    h5 = f['/PCCOMP/']
    h5misc = f['/MISC/']

    pc = PCCOMP()
    pcmisc = PCMISC()
    # Load
    pc.loadh5(h5)
    pcmisc.loadh5(h5misc)

    # Close HDF file
    f.close()

    # Display option
    print "PCCOMP-------------------"
    pc.display()
    print "MISC---------------------"
    pcmisc.display()

    # Plots radiance
    import matplotlib.pyplot as plt

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = 0
    for v in ['PCSCORES']:
        color = colors[i]
        i += 1
        plt.plot(pc[v], color + '-', label=v)
        plt.xlabel('scores (minus one)')
        plt.ylabel('PC' + '   (' + pc[v + '_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(fileName)
        plt.show()

    for v in ['BT_PCCOMP']:
        color = colors[i]
        i += 1
        plt.plot(pc[v], color + '-', label=v)
        plt.xlabel('channel (minus one)')
        plt.ylabel('Reconstitued BT' +
                   '   (' + pc[v + '_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(fileName)
        plt.show()

    for v in ['TOTAL_PCCOMP']:
        color = colors[i]
        i += 1
        plt.plot(pc[v], color + '-', label=v)
        plt.xlabel('channel (minus one)')
        plt.ylabel('Total_PCCOMP' +
                   '   (' + pc[v + '_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(fileName)
        plt.show()
