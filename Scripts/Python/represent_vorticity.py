from base_representer import BaseRepresenter
from matplotlib.ticker import OldScalarFormatter
from matplotlib import pyplot as pp

class VorticityRepresenter(BaseRepresenter):

    _name_columns = ['time', 'omegax', 'omegay', 'omegaz']

    def __init__(self):
        BaseRepresenter.__init__()


    def __init__(self):
        BaseRepresenter.__init__(self)
        pass


    def draw(self):
        data = self.data
        data['time'] -= min(data['time'])
        data = data[data['omegaz'] < 1e2]
        print(data.columns)

        ax = data.plot(x='time', y=['omegax', 'omegay', 'omegaz'])
        # set parameters for plotting
        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 12

        # ax.set_title(folder_name)#+',  toroidal/total energy ratio: '+ string_ratio+'%')
        ax.set_xlabel('t')
        ax.set_ylabel('omega')

        BaseRepresenter.draw(self)


    pass

if __name__ =="__main__":
    reader = VorticityRepresenter()
    reader.open()
    reader.draw()