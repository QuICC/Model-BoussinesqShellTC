import numpy as np
from matplotlib import pyplot as pp
import sys
import pandas as pd
from geomhdiscc.geometry.spherical import shell_radius

if __name__=="__main__":

    try:
        argv = sys.argv
        print argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python represent_spectra.py filename')
        sys.exit()

    datafull = np.loadtxt(filename, skiprows=3)[:,1:]

    #data = np.mean(data,axis=0)

    #pp.figure()
    #pp.ion()
    """ for i in range(datafull.shape[0]-1,datafull.shape[0]):
        pp.clf()
        """
    data = datafull[-1,:]

    print (data.shape)
    Lmax = data.shape[0]/4

    pp.loglog(data[0:Lmax], label='L spectrum, toroidal')
    pp.loglog(data[Lmax:2*Lmax], label='M spectrum, toroidal')
    pp.loglog(data[2*Lmax:3*Lmax], label='L spectrum, poloidal')
    pp.plot(data[3*Lmax:-1], label='M spectrum, poloidal')

    pp.legend()
    pp.show()




