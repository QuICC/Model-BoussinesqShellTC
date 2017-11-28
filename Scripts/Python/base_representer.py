# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:20 2017

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
import os
import pandas as pd
from matplotlib import pyplot as pp
import xml.etree.ElementTree as ET


class BaseRepresenter:

    def __init__(self):
        self.data = None

        # use an incrementing indices so that draw can be used for different plots
        self.idx_draw = 2

        pass

    # define the name columns static variable
    _name_columns = None
    @property
    def name_columns(self):
        if self._name_columns == None:
            raise NotImplementedError('Name Columns not implemented in the base class')
        return self._name_columns

    def open(self):

        try:
            argv = sys.argv
            # print argv
            self.filename = argv[1]
        except RuntimeError as e:
            print(e)
            print('Supposed usage: python represent_energies.py filename')
            sys.exit()

        try:

            folder_name = os.path.relpath(".", "..")
            data = pd.read_csv(self.filename, sep='\t', skiprows=3, names=self.name_columns)

        except IOError as e:
            folder_name = os.path.relpath(".", "..")
            data = []
            for folder in os.listdir('.'):
                if (os.path.isdir(folder)):
                    # print(folder)
                    try:
                        datatemp = pd.read_csv(folder + '/' + self.filename, sep='\t', skiprows=3,
                                               names=self.name_columns)
                        print(folder, data)
                        data.append(datatemp)
                    except IOError as e:
                        #print(e)
                        pass
            print(data)
            data = pd.concat(data, ignore_index=True)

            data.reindex()
            # print(data.head())
            data.sort_values([r'$t$'], inplace=True)

        self.data = data

    def draw(self):

        try:
            fout = sys.argv[self.idx_draw]
            pp.savefig(fout)
            self.idx_draw+=1
        except IndexError as e:
            pp.show()

    def search_in_parameter(self, param):

        dirname = os.path.dirname(self.filename)
        print(dirname)
        cfg_file = open(dirname + '/parameters.cfg', 'r')
        header = cfg_file.readline()
        root = ET.fromstring(header + '<root>' + cfg_file.read() + '</root>')
        leaf = root.find('*/*/'+param)
        return leaf.text
