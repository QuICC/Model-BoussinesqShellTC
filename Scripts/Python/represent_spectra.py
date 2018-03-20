import numpy as np
from matplotlib import pyplot as pp
import sys, os
import pandas as pd
#from quicc.geometry.spherical import shell_radius
from base_representer import BaseRepresenter
import xml.etree.ElementTree as ET

class SpectraRepresenter(BaseRepresenter):

    _name_columns = ['L spectrum, toroidal', 'M spectrum, toroidal', 'L spectrum, poloidal', 'M spectrum, poloidal']
    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def open(self, filename=None):

        if filename==None:
            try:
                argv = sys.argv
                self.filename = argv[1]
            except RuntimeError as e:
                print(e)
                print('Supposed usage: python represent_spectra.py filename')
                sys.exit()
        else:
            self.filename=filename


        try:
            #folder_name = os.path.relpath("..","../..")
            #folder_name = os.path.relpath(".", "..")
            folderpath = os.path.dirname(self.filename)
            datafull = np.loadtxt(self.filename, skiprows=3)
        except IOError as e:

            folder_name = os.path.relpath(".", "..")
            datafull = []
            for folder in os.listdir(folderpath):
                if (os.path.isdir(folderpath+'/'+folder)):
                    print(folder)
                    try:
                        print(folderpath+'/'+folder + '/' + os.path.basename(self.filename))
                        datatemp = pd.DataFrame(np.loadtxt(folderpath+'/'+folder + '/' + os.path.basename(self.filename), skiprows=3))
                        datafull.append(datatemp)
                        print(datatemp)
                    except IOError as e:
                        # print(e)
                        pass

            datafull = pd.concat(datafull, ignore_index=True)
            datafull.reindex()
            datafull.sort_values([0], inplace=True)

            datafull = datafull.as_matrix()
        self.datafull = datafull

    def draw(self,type='both'):

        pp.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 14
        # set parameters for plotting
        pp.ticklabel_format(style='sci', axis='y')

        #print(type(self.datafull))
        print(type)
        self.draw_snapshot(self.datafull,type=type)

    def draw_snapshot(self, datafull, **kwargs):
        #pp.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 14
        # set parameters for plotting
        #pp.ticklabel_format(style='sci', axis='y')

        data = datafull[-1, 1:]
        #data = datafull[ 1:]

        try:

            dirname = os.path.dirname(self.filename)
            print(dirname)
            cfg_file = open(dirname + '/parameters.cfg', 'r')
            header = cfg_file.readline()
            root = ET.fromstring(header + '<root>' + cfg_file.read() + '</root>')
            Lmax = int(root[1][0][1].text)+1
            Mmax = int(root[1][0][2].text)+1
            print('Parameter file found', Lmax, Mmax)
        except BaseException as e:
            Lmax = Mmax = data.shape[0] / 4
            pass


        if kwargs['type']!='l':
            pp.loglog(np.cumsum(np.ones_like(data[:Lmax])), data[:Lmax], label='L spectrum, toroidal')
            pp.loglog(np.cumsum(np.ones_like(data[Mmax + Lmax:Mmax + 2 * Lmax])), data[Mmax + Lmax:Mmax + 2 * Lmax],
                      label='L spectrum, poloidal')

        if kwargs['type']!='m':
            pp.loglog(np.cumsum(np.ones_like(data[Lmax:Lmax + Mmax])), data[Lmax:Lmax + Mmax], label='M spectrum, toroidal')
            pp.loglog(np.cumsum(np.ones_like(data[Mmax+ 2 * Lmax:])), data[Mmax+ 2 * Lmax:], label='M spectrum, poloidal')
            
        pp.ylim(ymin=1e-20)

        """
        pp.figure()
        pp.semilogy(data[0:Lmax], label='L spectrum, toroidal')
        pp.semilogy(data[Lmax:Lmax + Mmax], label='M spectrum, toroidal')
        pp.semilogy(data[Mmax + Lmax:Mmax+ 2 * Lmax], label='L spectrum, poloidal')
        pp.semilogy(data[Mmax+ 2 * Lmax:], label='M spectrum, poloidal')
        # pp.title(folder_name)
        """
        pp.xlabel(r'$l+1 \quad m+1$')
        pp.ylabel(r'$E_{kin}$')
        pp.legend(prop={'size': 14})

        BaseRepresenter.draw(self)

if __name__=="__main__":
    reader = SpectraRepresenter()
    print(reader.name_columns)
    reader.open()
    reader.draw()










