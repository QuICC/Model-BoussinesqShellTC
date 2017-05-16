# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
import os
import pandas as pd
from matplotlib import pyplot as pp
from matplotlib.ticker import OldScalarFormatter


if __name__=="__main__":



    # open files

    try:
        argv = sys.argv
        #print argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python represent_energies.py filename')
        sys.exit()



    try:

        #folder_name = os.path.relpath("..","../..")
        folder_name = os.path.relpath(".", "..")
        data = pd.read_csv(filename, sep='\t', skiprows=3, names=['time', 'total', 'toroidal', 'poloidal'])

    except IOError as e:
        folder_name = os.path.relpath(".","..")
        data = []
        for folder in os.listdir('.'):
            if (os.path.isdir(folder)):
                #print(folder)
                try:
                    datatemp = pd.read_csv(folder + '/' + filename, sep='\t', skiprows=3, names=['time', 'total', 'toroidal', 'poloidal'])
                    data.append(datatemp)
                except IOError as e:
                    #print(e)
                    pass



        data = pd.concat(data, ignore_index=True)
        data.reindex()
        #print(data.head())
        data.sort_values(['time'], inplace=True)
        pass

    data['time'] -= min(data['time'])
    data = data[data['total']<1e2]
    
    idx = max(data.index)
    string_ratio = str("%.2f" % (data['toroidal'][idx]/data['total'][idx]*100))
    #print('Final Toroidal to Total energy ratio of: '+string_ratio)

    ax = data.plot(x='time',y=['total', 'toroidal','poloidal'])
    # set parameters for plotting
    ax.yaxis.set_major_formatter(OldScalarFormatter())
    pp.rcParams['font.family'] = 'ubuntu'
    #pp.rcParams['font.size'] = 12

    #ax.set_title(folder_name)#+',  toroidal/total energy ratio: '+ string_ratio+'%')
    ax.set_xlabel('t')
    ax.set_ylabel('E')


    try:
        fout = argv[2]
        pp.savefig(fout)
    except IndexError as e:
        pp.show()
    
    
    
    
