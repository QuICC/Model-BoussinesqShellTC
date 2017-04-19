import numpy as np
from matplotlib import pyplot as pp
import sys, os
import pandas as pd
from quicc.geometry.spherical import shell_radius

if __name__=="__main__":

    try:
        argv = sys.argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python represent_spectra.py filename')
        sys.exit()


    try:
        #folder_name = os.path.relpath("..","../..")
        folder_name = os.path.relpath(".", "..")
        datafull = np.loadtxt(filename, skiprows=3)
    except IOError as e:

        folder_name = os.path.relpath(".", "..")
        datafull = []
        for folder in os.listdir('.'):
            if (os.path.isdir(folder)):
                #print(folder)
                try:
                    datatemp =  pd.DataFrame(np.loadtxt(folder + '/' + filename, skiprows=3))
                    datafull.append(datatemp)
                except IOError as e:
                    #print(e)
                    pass


        #datafull = np.vstack(data)
        #datafull = pd.DataFrame(datafull).sort([0])
        datafull = pd.concat(datafull, ignore_index=True)
        datafull.reindex()
        datafull.sort([0],inplace=True)

        datafull = datafull.as_matrix()
    pass
    pp.figure()
    pp.clf()

    data = datafull[-1,1:]

    Lmax = data.shape[0]/4

    pp.loglog(data[0:Lmax], label='L spectrum, toroidal')
    pp.loglog(data[Lmax:2*Lmax], label='M spectrum, toroidal')
    pp.loglog(data[2*Lmax:3*Lmax], label='L spectrum, poloidal')
    pp.plot(data[3*Lmax:-1], label='M spectrum, poloidal')

    pp.title(folder_name)
    pp.xlabel('l/m')
    pp.ylabel('E')
    pp.legend()
        #pp.draw()
    try:
        fout = argv[2]
        pp.savefig(fout)
    except IndexError as e:
        pp.show()




