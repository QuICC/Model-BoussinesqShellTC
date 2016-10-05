# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""

import pandas as pd
from matplotlib import pyplot as pp

if __name__="__main__":
    
    table = pd.read_csv(filename,skiprows=3,header=['time', 'total' 'toroidal', 'poloida'], sep = '\t')
    table.plot(x='time',y=['total', 'toroidal','poloidal'])
    
    print('Final Toroidal to Total energy ratio of: '+str(table['toroidal'][-1]/table['total'][-1]))
    
    
    