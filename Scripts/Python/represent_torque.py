import numpy as np
from matplotlib import pyplot as pp
import sys
import os
import pandas as pd
import seaborn as sns
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels  import RBF, WhiteKernel
#from statsmodels.api.nonparametric import lowess
from sklearn.kernel_ridge import KernelRidge
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
        data = pd.read_csv(filename, sep='\t', skiprows=3, names=['time', 'Torque'])

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
        data.sort_values(['time'],inplace=True)
        pass
    data['time']-= min(data['time'])
    data = data[abs(data['Torque'])<1.0]


    #e# ax = data.plot(x='time', y='Torque', title=folder_name)
    ax = data.plot(x='time', y = 'Torque', alpha = 0.25)

    """
    # train a Gaussian Process Regressor (expensive)
    data = data.values
    kernel = 1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-10, 1e+1))
    gp = GaussianProcessRegressor()
    gp.fit(data[:,[0]],data[:,1])
    y_predict , y_std = gp.predict(data[:,[0]],return_std = True)

    ax.plot(data[:,0],y_predict)
    ax.fill_between(data[:,0], y_predict - y_std, y_predict + y_std,alpha=.25)
    """
    """
    # train a Kernel Ridge regression
    data = data
    kr = KernelRidge()
    kr.fit(data[:,[0]],data[:,1])
    y_predict , y_std = kr.predict(data[:,[0]],return_std = True)

    ax.plot(data[:,0],y_predict)
    ax.fill_between(data[:,0], y_predict - y_std, y_predict + y_std,alpha=.25)
    """
    alpha = 0.05
    data.set_index('time', inplace = True)
    # forward
    ewm = data['Torque'].ewm(alpha=alpha, adjust=True)
    m = ewm.agg(['mean','std'])
    # backward
    ewm_bwd= data['Torque'][::-1].ewm(alpha = alpha,adjust =True)
    m_bwd = ewm_bwd.agg(['mean','std'])
    print(m[::-1].head())
    m = (m+m_bwd)/2.
    ax = m['mean'].plot()
    ax.fill_between(m.index, m['mean'] - m['std'], m['mean'] + m['std'], alpha=.5)


    ax.set_title(folder_name)
    ax.set_xlabel('t')


    try:
        fout = argv[2]
        pp.savefig(fout)
    except IndexError as e:
        pp.show()
