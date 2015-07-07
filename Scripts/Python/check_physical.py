#!/usr/bin/env python

from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function

import sys, getopt

import h5py
import numpy as np
import re

def getOrdering(scheme):
    """Get data ordering for scheme"""

    if scheme in [b'TTT']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'TFT']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'TFF']:
        idSlow = 'z'
        idMid = 'x'
        idFast = 'y'
    elif scheme in [b'FFF']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'AFT', b'CFT', b'WFT']:
        idSlow = 'r'
        idMid = 'theta'
        idFast = 'z'
    elif scheme in [b'SLFm', b'SLFl', b'BLFl', b'BLFm', b'WLF']:
        idSlow = 'r'
        idMid = 'theta'
        idFast = 'phi'

    return (idFast, idMid, idSlow)

def getGrid(scheme, h5_file):
    """Get the grid for scheme"""

    idFast, idMid, idSlow = getOrdering(scheme)

    gFast = h5_file['mesh']['grid_'+idFast]
    gMid = h5_file['mesh']['grid_'+idMid]
    gSlow = h5_file['mesh']['grid_'+idSlow]

    return (gFast, gMid, gSlow)

def getGridSize(scheme, h5_file):
    """Get the grid for scheme"""

    idFast, idMid, idSlow = getOrdering(scheme)

    nFast = h5_file['mesh']['grid_'+idFast].size
    nMid = h5_file['mesh']['grid_'+idMid].size
    nSlow = h5_file['mesh']['grid_'+idSlow].size

    return (nFast, nMid, nSlow)

def checkShellBenchmarkC0(h5_file):
    """Check initial state for sphere benchmark C0"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']

    phi = gFast[()]
    ro = h5_file['/physical/ro'][()]
    ri = h5_file['/physical/rratio'][()]*ro
    itS = np.nditer(gSlow, flags=['c_index'])
    while not itS.finished:
        r = itS[0]
        x = 2.0*r - ri  - ro
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = np.max(np.abs(dataT[itS.index, itM.index,:] - (21.0/np.sqrt(17920*np.pi)*(1.0 - 3*x**2 + 3.0*x**4 - x**6)*np.sin(theta)**4*np.cos(4.0*phi))))

            # Check velocity field
            errVr = np.max(np.abs(dataVr[itS.index, itM.index,:]))
            errVth = np.max(np.abs(dataVth[itS.index, itM.index,:]))
            errVph = np.max(np.abs(dataVph[itS.index, itM.index,:]))
            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in temperature field:     {:g}".format(errT))

def checkShellBenchmarkC1(h5_file):
    """Check initial state for sphere benchmark C2"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']
    dataBr = h5_file['/magnetic/magnetic_r']
    dataBth = h5_file['/magnetic/magnetic_theta']
    dataBph = h5_file['/magnetic/magnetic_phi']

    phi = gFast[()]
    ro = h5_file['/physical/ro'][()]
    ri = h5_file['/physical/rratio'][()]*ro
    itS = np.nditer(gSlow, flags=['c_index'])
    while not itS.finished:
        r = itS[0]
        x = 2.0*r - ri  - ro
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = np.max(np.abs(dataT[itS.index, itM.index,:] - (21.0/np.sqrt(17920*np.pi)*(1.0 - 3*x**2 + 3.0*x**4 - x**6)*np.sin(theta)**4*np.cos(4.0*phi))))

            # Check velocity field
            errVr = np.max(np.abs(dataVr[itS.index, itM.index,:]))
            errVth = np.max(np.abs(dataVth[itS.index, itM.index,:]))
            errVph = np.max(np.abs(dataVph[itS.index, itM.index,:]))

            # Check magnetic field
            scale = 0.5
            errBr = np.max(np.abs(dataBr[itS.index, itM.index,:] - (scale*5.0/8.0*(8.0*ro - 6.0*r - 2.0*ri**4/r**3)*np.cos(theta))))
            errBth = np.max(np.abs(dataBth[itS.index, itM.index,:] - (scale*5.0/8.0*(9.0*r - 8.0*ro - ri**4/r**3)*np.sin(theta))))
            errBph = np.max(np.abs(dataBph[itS.index, itM.index,:] - (scale*5.0*np.sin(np.pi*(r - ri))*np.sin(2.0*theta))))

            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in radial magnetic field: {:g}".format(errBr))
    print("Error in theta magnetic field:  {:g}".format(errBth))
    print("Error in phi magnetic field:    {:g}".format(errBph))
    print("Error in temperature field:     {:g}".format(errT))

def checkSphereBenchmarkC1(h5_file):
    """Check initial state for sphere benchmark C1"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataDTr = h5_file['/temperature_grad/temperature_grad_r']
    dataDTth = h5_file['/temperature_grad/temperature_grad_theta']
    dataDTph = h5_file['/temperature_grad/temperature_grad_phi']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']

    phi = gFast[()]
    itS = np.nditer(gSlow, flags=['c_index'])
    while not itS.finished:
        r = itS[0]
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = np.max(np.abs(dataT[itS.index, itM.index,:] - (0.5*(1 - r**2) + 1e-5/8.*np.sqrt(35/np.pi)*r**3*(1 - r**2)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))))

            # Check temperature gradient
            errDTr = np.max(np.abs(dataDTr[itS.index, itM.index,:] - (-r + 1e-5/8.*np.sqrt(35/np.pi)*(3.0*r**2 - 5.0*r**4)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))))
            errDTth = np.max(np.abs(dataDTth[itS.index, itM.index,:] - (-3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.cos(theta)*np.sin(theta)**2*(np.cos(3.0*phi) + np.sin(3.0*phi)))))
            errDTph = np.max(np.abs(dataDTph[itS.index, itM.index,:] - (3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.sin(theta)**2*(np.cos(3.0*phi) - np.sin(3.0*phi)))))

            # Check velocity field
            errVr = np.max(np.abs(dataVr[itS.index, itM.index,:]))
            errVth = np.max(np.abs(dataVth[itS.index, itM.index,:]))
            errVph = np.max(np.abs(dataVph[itS.index, itM.index,:]))
            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in temperature field:     {:g}".format(errT))
    print("Error in radial temperature gradient:     {:g}".format(errDTr))
    print("Error in theta temperature gradient:     {:g}".format(errDTth))
    print("Error in phi temperature gradient:     {:g}".format(errDTph))

def checkSphereBenchmarkC2(h5_file):
    """Check initial state for sphere benchmark C2"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    # Temperature field and gradient
    dataT = h5_file['/temperature/temperature']
    dataDTr = h5_file['/temperature_grad/temperature_grad_r']
    dataDTth = h5_file['/temperature_grad/temperature_grad_theta']
    dataDTph = h5_file['/temperature_grad/temperature_grad_phi']
    # Velocity field and curl
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']
    dataCVr = h5_file['/velocity_curl/velocity_curl_r']
    dataCVth = h5_file['/velocity_curl/velocity_curl_theta']
    dataCVph = h5_file['/velocity_curl/velocity_curl_phi']
    # Magnetic field and curl
    dataBr = h5_file['/magnetic/magnetic_r']
    dataBth = h5_file['/magnetic/magnetic_theta']
    dataBph = h5_file['/magnetic/magnetic_phi']
    dataCBr = h5_file['/magnetic_curl/magnetic_curl_r']
    dataCBth = h5_file['/magnetic_curl/magnetic_curl_theta']
    dataCBph = h5_file['/magnetic_curl/magnetic_curl_phi']

    phi = gFast[()]
    itS = np.nditer(gSlow, flags=['c_index'])
    while not itS.finished:
        r = itS[0]
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = np.max(np.abs(dataT[itS.index, itM.index,:] - (0.5*(1 - r**2) + 1e-5/8.*np.sqrt(35/np.pi)*r**3*(1 - r**2)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))))

            # Check temperature gradient
            errDTr = np.max(np.abs(dataDTr[itS.index, itM.index,:] - (-r + 1e-5/8.*np.sqrt(35/np.pi)*(3.0*r**2 - 5.0*r**4)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))))
            errDTth = np.max(np.abs(dataDTth[itS.index, itM.index,:] - (-3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.cos(theta)*np.sin(theta)**2*(np.cos(3.0*phi) + np.sin(3.0*phi)))))
            errDTph = np.max(np.abs(dataDTph[itS.index, itM.index,:] - (3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.sin(theta)**2*(np.cos(3.0*phi) - np.sin(3.0*phi)))))

            # Check velocity field
            errVr = np.max(np.abs(dataVr[itS.index, itM.index,:]))
            errVth = np.max(np.abs(dataVth[itS.index, itM.index,:] - (-10*r**2/(7.0*np.sqrt(3.0))*np.cos(theta)*(3.0*(-147.0 + 343.0*r**2 - 217.0*r**4 + 29*r**6)*np.cos(phi) + 14.0*(-9.0 - 125.0*r**2 + 39*r**4 + 27*r**6)*np.sin(phi)))))
            errVph = np.max(np.abs(dataVph[itS.index, itM.index,:] - (-5.0*r/5544.*(7.0*((43700.0 - 58113*r**2 - 15345*r**4 + 1881*r**6 + 20790*r**8)*np.sin(theta) + 1485*r**2*(-9.0 + 115.0*r**2 - 167.0*r**4 + 70*r**6)*np.sin(3.0*theta)) + 528.0*np.sqrt(3)*r*np.cos(2.0*theta)*(14.0*(-9.0 - 125.0*r**2 + 39.0*r**4 + 27.0*r**6)*np.cos(phi) + 3.0*(147.0 - 343.0*r**2 + 217*r**4 - 29.0*r**6)*np.sin(phi))))))
            
            # Check curl of velocity field
            errCVr = np.max(np.abs(dataCVr[itS.index, itM.index,:] - ((1.0/1386.0)*(-35.0*(21850.0 + 99.0*r**2*(-361.0 + 785.0*r**2 - 1243.0*r**4 + 630.0*r**6))*np.cos(theta) + 495.0*r*(-105.0*r*(-9.0 + 115.0*r**2 - 167.0*r**4 + 70.0*r**6)*np.cos(3.0*theta) + 4.0*np.sqrt(3)*np.sin(2.0*theta)*(14.0*(-9.0 - 125.0*r**2 + 39.0*r**4 + 27.0*r**6)*np.cos(phi) + 3.0*(147.0 - 343.0*r**2 + 217.0*r**4 - 29.0*r**6)*np.sin(phi)))))))
            errCVth = np.max(np.abs(dataCVth[itS.index, itM.index,:] - ((54625.0*np.sin(theta))/99.0 + 5.0/84.0*r*(21.0*r*((-1174.0 - 465.0*r**2 + 76.0*r**4 + 1050.0*r**6)*np.sin(theta) + 15.0*(-18.0 + 345.0*r**2 - 668.0*r**4 + 350.0*r**6)*np.sin(3.0*theta)) + 8.0*np.sqrt(3)*np.cos(2.0*theta)*(14.0*(-27.0 - 625.0*r**2 + 273.0*r**4 + 243.0*r**6)*np.cos(phi) + 3.0*(441.0 - 1715.0*r**2 + 1519.0*r**4 - 261.0*r**6)*np.sin(phi))))))
            errCVph = np.max(np.abs(dataCVph[itS.index, itM.index,:] - (10.0*r/(7.0*np.sqrt(3))*np.cos(theta)*((1323.0 - 5145.0*r**2 + 4557.0*r**4 - 783.0*r**6)*np.cos(phi) + 14.0*(27.0 + 625.0*r**2 - 273.0*r**4 - 243.0*r**6)*np.sin(phi)))))
            

            # Check magnetic field
            errBr = np.max(np.abs(dataBr[itS.index, itM.index,:]))
            errBth = np.max(np.abs(dataBth[itS.index, itM.index,:] - (-3.0/2.0*r*(-1.0 + 4.0*r**2 - 6.0*r**4 + 3.0*r**6)*(np.cos(phi) + np.sin(phi)))))
            errBph = np.max(np.abs(dataBph[itS.index, itM.index,:] - (-3.0/4.0*r*(-1.0 + r**2)*np.cos(theta)*(3.0*r*(2.0 - 5.0*r**2 + 4.0*r**4)*np.sin(theta) + 2.0*(1.0 - 3.0*r**2 + 3.0*r**4)*(np.cos(phi) - np.sin(phi))))))

            # Check curl of magnetic field
            errCBr = np.max(np.abs(dataCBr[itS.index, itM.index,:] - (-9.0/8.0*r*(-2.0 + 7.0*r**2 - 9.0*r**4 + 4.0*r**6)*(1.0 + 3.0*np.cos(2.0*theta)) + 3.0*(-1.0 + 4*r**2 - 6.0*r**4 + 3.0*r**6)*np.sin(theta)*np.cos(phi) + 3.0*(1.0 - 4.0*r**2 + 6.0*r**4 - 3.0*r**6)*np.sin(theta)*np.sin(phi))))
            errCBth = np.max(np.abs(dataCBth[itS.index, itM.index,:] - (3.0/4.0*np.cos(theta)*(3.0*r*(-6.0 + 35.0*r**2 - 63.0*r**4 + 36.0*r**6)*np.sin(theta) + 4.0*(-1.0 + 8.0*r**2 - 18.0*r**4 + 12.0*r**6)*(np.cos(phi) - np.sin(phi))))))
            errCBph = np.max(np.abs(dataCBph[itS.index, itM.index,:] - (-3.0*(-1.0 + 8.0*r**2 - 18.0*r**4 + 12.0*r**6)*(np.cos(phi) + np.sin(phi)))))

            itM.iternext()
        itS.iternext()

    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in radial curl of velocity field: {:g}".format(errCVr))
    print("Error in theta curl of velocity field:  {:g}".format(errCVth))
    print("Error in phi curl of velocity field:    {:g}".format(errCVph))

    print("Error in radial magnetic field: {:g}".format(errBr))
    print("Error in theta magnetic field:  {:g}".format(errBth))
    print("Error in phi magnetic field:    {:g}".format(errBph))
    print("Error in radial curl of magnetic field: {:g}".format(errCBr))
    print("Error in theta curl of magnetic field:  {:g}".format(errCBth))
    print("Error in phi curl of magnetic field:    {:g}".format(errCBph))

    print("Error in temperature field:     {:g}".format(errT))
    print("Error in radial temperature gradient:     {:g}".format(errDTr))
    print("Error in theta temperature gradient:     {:g}".format(errDTth))
    print("Error in phi temperature gradient:     {:g}".format(errDTph))

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        print('check_physical.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_physical.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg

    # Read HDF5 file
    h5_file = h5py.File(inputfile, 'r')
    # Check state
    #checkShellBenchmarkC0(h5_file)
    #checkShellBenchmarkC1(h5_file)
    #checkSphereBenchmarkC1(h5_file)
    checkSphereBenchmarkC2(h5_file)

    # Close file
    h5_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
