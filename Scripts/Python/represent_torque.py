import numpy as np
from matplotlib import pyplot as pp

from matplotlib.ticker import OldScalarFormatter

from base_representer import BaseRepresenter

class TorqueRepresenter(BaseRepresenter):

    _name_columns = [r'$t$', 'value']
    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def draw(self):
        data = self.data
        data[r'$t$'] -= min(data[r'$t$'])
        data = data[abs(data['value']) < 1.0]

        # post-processing on the value, to be removed once the value calculations are correct
        #ri = 0.35 / (1. - 0.35)
        #T = np.sqrt(4 * np.pi / 3) * ri
        #data['value'] = (data['value'] + ri * T * 8 * np.pi / 3. * 1e-5) * ri

        # e# ax = data.plot(x='time', y='value', title=folder_name)
        ax = data.plot(x=r'$t$', y='value', alpha=0.25)

        # set parameters for plotting

        # pp.rcParams['font.size'] = 14

        alpha = 0.05
        data.set_index(r'$t$', inplace=True)
        # forward
        ewm = data['value'].ewm(alpha=alpha, adjust=True)
        m = ewm.agg(['mean', 'std'])
        # backward
        ewm_bwd = data['value'][::-1].ewm(alpha=alpha, adjust=True)
        m_bwd = ewm_bwd.agg(['mean', 'std'])
        print(m[::-1].head())
        m = (m + m_bwd) / 2.
        ax = m['mean'].plot(legend='mean')
        ax.fill_between(m.index, m['mean'] - m['std'], m['mean'] + m['std'], alpha=.1, label='std')

        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        ax.set_ylabel(r'$G$')
        # ax.set_title(folder_name)
        #ax.set_xlabel('t')

        BaseRepresenter.draw(self)



if __name__=="__main__":

    reader = TorqueRepresenter()
    reader.open()
    reader.draw()