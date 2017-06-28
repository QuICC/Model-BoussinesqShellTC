import numpy as np
from matplotlib import pyplot as pp
import sys, os
import pandas as pd
from quicc.geometry.spherical import shell_radius
from base_representer import BaseRepresenter

class SpectraRepresenter(BaseRepresenter):

    _name_columns = ['L spectrum, toroidal', 'M spectrum, toroidal', 'L spectrum, poloidal', 'M spectrum, poloidal']
    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def open(self):

        try:
            argv = sys.argv
            filename = argv[1]
        except RuntimeError as e:
            print(e)
            print('Supposed usage: python represent_spectra.py filename')
            sys.exit()

        try:
            # folder_name = os.path.relpath("..","../..")
            folder_name = os.path.relpath(".", "..")
            datafull = np.loadtxt(filename, skiprows=3)
        except IOError as e:

            folder_name = os.path.relpath(".", "..")
            datafull = []
            for folder in os.listdir('.'):
                if (os.path.isdir(folder)):
                    # print(folder)
                    try:
                        datatemp = pd.DataFrame(np.loadtxt(folder + '/' + filename, skiprows=3))
                        datafull.append(datatemp)
                    except IOError as e:
                        # print(e)
                        pass

            datafull = pd.concat(datafull, ignore_index=True)
            datafull.reindex()
            datafull.sort_values([0], inplace=True)

            datafull = datafull.as_matrix()
            self.datafull = datafull

    def draw(self):

        pp.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 14
        # set parameters for plotting
        pp.ticklabel_format(style='sci', axis='y')

        print(type(self.datafull))


        self.draw_snapshot(self.datafull)

    def draw_snapshot(self, datafull):
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

        BaseRepresenter.draw(self)

if __name__=="__main__":
    reader = SpectraRepresenter()
    print(reader.name_columns)
    reader.open()
    reader.draw()










