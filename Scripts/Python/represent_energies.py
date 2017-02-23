# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
import os
import pandas as pd
from matplotlib import pyplot as pp

if __name__=="__main__":

    try:
        argv = sys.argv
        print argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python represent_energies.py filename')
        sys.exit()



    try:

        data = pd.read_csv(filename, sep='\t', skiprows=3, names=['time', 'total', 'toroidal', 'poloidal'])

    except IOError as e:
        data = []
        for folder in os.listdir('.'):
            if (os.path.isdir(folder)):
                print(folder)
                try:
                    datatemp = pd.read_csv(folder + '/' + filename, sep='\t', skiprows=3, names=['time', 'total', 'toroidal', 'poloidal'])
                except IOError as e:
                    print(e)


                data.append(datatemp)
        data = pd.concat(data, ignore_index=True)
        data.reindex()
        print(data.head())
        data.sort(['time'], inplace=True)
        pass

    data = data[data['total']<1e2]

    data.plot(x='time',y=['total', 'toroidal','poloidal'])
    
    idx = max(data.index)
    print('Final Toroidal to Total energy ratio of: '+str(data['toroidal'][idx]/data['total'][idx]))
    
    pp.show()
    
    
    
    
