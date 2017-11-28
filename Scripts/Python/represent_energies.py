# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
from base_representer import BaseRepresenter
from matplotlib.ticker import OldScalarFormatter
from matplotlib import pyplot as pp

class EnergyRepresenter(BaseRepresenter):

    _name_columns = [r'$t$', 'total', 'toroidal', 'poloidal']

    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def draw(self):
        data = self.data
        data[r'$t$'] -= min(data[r'$t$'])
        data = data[data['total'] < 1e2]

        idx = max(data.index)
        string_ratio = str("%.2f" % (data['toroidal'][idx] / data['total'][idx] * 100))
        # print('Final Toroidal to Total energy ratio of: '+string_ratio)

        ax = data.plot(x=r'$t$', y=['total', 'toroidal', 'poloidal'])
        # set parameters for plotting
        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 12

        # ax.set_title(folder_name)#+',  toroidal/total energy ratio: '+ string_ratio+'%')
        ax.set_xlabel('t')
        ax.set_ylabel('E')

        BaseRepresenter.draw(self)

    pass


if __name__=="__main__":

    reader = EnergyRepresenter()
    reader.open()
    reader.draw()







    
    
    
    
