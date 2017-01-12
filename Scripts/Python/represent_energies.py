# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
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
    
    table = pd.read_csv(filename, skiprows=3, names=['time', 'total', 'toroidal', 'poloidal'], sep = '\t')
    table.plot(x='time',y=['total', 'toroidal','poloidal'])
    
    idx = max(table.index)
    print('Final Toroidal to Total energy ratio of: '+str(table['toroidal'][idx]/table['total'][idx]))
    
    pp.show()
    
    
    
    
