# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:20 2017

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
import os
import pandas as pd
from matplotlib import pyplot as pp


class BaseRepresenter:

    def __init__(self):
        self.data = None
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
            filename = argv[1]
        except RuntimeError as e:
            print(e)
            print('Supposed usage: python represent_energies.py filename')
            sys.exit()

        try:

            folder_name = os.path.relpath(".", "..")
            data = pd.read_csv(filename, sep='\t', skiprows=3, names=self.name_columns)

        except IOError as e:
            folder_name = os.path.relpath(".", "..")
            data = []
            for folder in os.listdir('.'):
                if (os.path.isdir(folder)):
                    # print(folder)
                    try:
                        datatemp = pd.read_csv(folder + '/' + filename, sep='\t', skiprows=3,
                                               names=self.name_columns)
                        data.append(datatemp)
                    except IOError as e:
                        #print(e)
                        pass

            data = pd.concat(data, ignore_index=True)
            data.reindex()
            # print(data.head())
            data.sort_values(['time'], inplace=True)

        self.data = data

    def draw(self):

        try:
            fout = sys.argv[2]
            pp.savefig(fout)
        except IndexError as e:
            pp.show()
