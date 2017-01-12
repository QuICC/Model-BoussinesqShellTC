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

    data = np.loadtxt(filename, skiprows=3)[1:-1,:]

    Lmax = data.shape[1]/4

    pp.loglog(data[-1,0:Lmax], label='L spectrum, toroidal')
    pp.loglog(data[-1,Lmax:2*Lmax], label='M spectrum, toroidal')
    pp.loglog(data[-1,2*Lmax:3*Lmax], label='L spectrum, poloidal')
    pp.loglog(data[-1,3*Lmax:-1], label='M spectrum, poloidal')

    pp.legend()
    pp.show()




    """ old implementation
    # call linear_r2x
    ro = file['physical/ro'].value
    rratio = file['physical/rratio'].value
    a, b = shell_radius.linear_r2x(ro, rratio)

    # get the matrices
    nr = file['truncation/spectral/dim1D'].value+2
    bc = dict()
    bc[0]=0
    IntgOp = shell_radius.r2(nr, a, b, bc)
    SphIntgOp = shell_radius.integral(nr, a, b)


    pp.figure()
    to_select = ['velocity/velocity_tor', 'velocity/velocity_pol']
    df = pd.DataFrame()
    for field in to_select:
        a = file[field][:]*file[field][:]

        for aa in a:
            print(aa)
            print(np.dot(SphIntgOp,aa))
            raise RuntimeError
        vec = (a**2).sum(1).sum(1)
        index = 0
        lmax = int((vec.shape[0]*2+.25)**.5-.5)
        result = np.zeros(lmax)


        # for different m
        for l in range(lmax):
            indexnext = index+l
            temp = np.array(vec[index:indexnext])
            temp[0:0]/=2.
            result[0:l] += temp *l*(l+1)
            index = indexnext+1


        # for different l
        for l in range(lmax):
                    indexnext = index+l+1
                    ll = np.array(range(l+1))
                    ll = ll*(ll+1)
                    result[l] = (vec[index:indexnext]*ll).sum()

    """