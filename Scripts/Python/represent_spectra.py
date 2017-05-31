import numpy as np
from matplotlib import pyplot as pp
import sys, os
import pandas as pd
from quicc.geometry.spherical import shell_radius


def draw(datafull):
    #pp.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    pp.rcParams['font.family'] = 'ubuntu'
    # pp.rcParams['font.size'] = 14
    # set parameters for plotting
    #pp.ticklabel_format(style='sci', axis='y')

    data = datafull[-1, 1:]

    Lmax = data.shape[0] / 4

    pp.loglog(data[0:Lmax], label='L spectrum, toroidal')
    pp.loglog(data[Lmax:2 * Lmax], label='M spectrum, toroidal')
    pp.loglog(data[2 * Lmax:3 * Lmax], label='L spectrum, poloidal')
    pp.loglog(data[3 * Lmax:-1], label='M spectrum, poloidal')

    # pp.title(folder_name)
    pp.xlabel('l/m')
    pp.ylabel('E')
    pp.legend()
    # pp.draw()
    try:
        fout = argv[2]
        pp.savefig(fout)
    except IndexError as e:
        pp.show()

        return

if __name__=="__main__":


    # open files

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
        datafull.sort_values([0],inplace=True)

        datafull = datafull.as_matrix()
    pass

    try:
        if argv[2] == "-v":
            by = 10
            for i in range(0,datafull.shape[0], by):
                data = datafull[i:i+1,:]
                pp.clf()
                pp.ion()
                draw(data)
                pp.pause(0.05)

        else:
            #pp.imshow(datafull)
            Lmax= datafull.shape[1]/4
            pp.plot(datafull[:,0],np.log(datafull[:,3+Lmax])/np.log(10))
            pp.plot(datafull[:,0],np.log(datafull[:,2+Lmax])/np.log(10))
            pp.show()


    except IndexError as e:
        draw(datafull)
        pass








