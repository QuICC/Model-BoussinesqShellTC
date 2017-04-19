import numpy as np
from matplotlib import pyplot as pp
import sys
import os
import pandas as pd
from quicc.geometry.spherical import shell_radius
from pandas.tools.plotting import autocorrelation_plot

if __name__=="__main__":

    try:
        argv = sys.argv
        #print argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python represent_torque.py filename')
        sys.exit()

    try:

        folder_name = os.path.relpath("..","../..")
        data = pd.read_csv(filename, sep='\t', skiprows=3, names=['time', 'torque'])

    except IOError as e:
        folder_name = os.path.relpath(".","..")
        data = []
        for folder in os.listdir('.'):
            if (os.path.isdir(folder)):
                #print(folder)
                try:
                    datatemp = pd.read_csv(folder + '/' + filename, sep='\t', skiprows=3, names=['time', 'Torque'])
                    data.append(datatemp)
                except IOError as e:
                    #print e
                    pass


        data = pd.concat(data, ignore_index=True)
        data.reindex()
        #print(data.head())
        data.sort(['time'],inplace=True)
        pass
    data['time']-= min(data['time'])
    data = data[abs(data['Torque'])<1.0]

    ax = data.plot(x='time', y='Torque', title=folder_name)
    ax.set_title(folder_name)
    ax.set_xlabel('t')
    """
    pp.figure()

    data.set_index('time')
    autocorrelation_plot(data['Torque'])
    """

    try:
        fout = argv[2]
        pp.savefig(fout)
    except IndexError as e:
        pp.show()
