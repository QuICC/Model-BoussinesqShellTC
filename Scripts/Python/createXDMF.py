#!/usr/bin/env python

from __future__ import print_function

import sys, getopt

from math import *
import h5py
import numpy as np
import re

def main(argv):
    inputfile = ''
    outputfile = ''
    snapshots = 1 
    with_components = False
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:", ['with-components'])
    except getopt.GetoptError:
        print('Single file: createXDMF.py (--with-components) -i <inputfile> -o <outputfile>')
        print('Timeseries: createXDMF.py (--with-components) -i <inputfile> -o <outputfile> -n <number of snapshots>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Single file: createXDMF.py (--with-components) -i <inputfile> -o <outputfile>')
            print('Timeseries: createXDMF.py (--with-components) -i <inputfile> -o <outputfile> -n <number of snapshots>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("--with-components"):
            with_components = True
        elif opt in ("-n"):
            snapshots = int(arg)
    # Extract file information
    h5_file = h5py.File(inputfile, 'r')
    scheme = h5_file['/'].attrs['type']
    if scheme in [b'TTT']:
        gSlow = 'x'
        gMid = 'y'
        gFast = 'z'
    elif scheme in [b'TFT']:
        gSlow = 'x'
        gMid = 'y'
        gFast = 'z'
    elif scheme in [b'TFF']:
        gSlow = 'z'
        gMid = 'x'
        gFast = 'y'
    elif scheme in [b'FFF']:
        gSlow = 'x'
        gMid = 'y'
        gFast = 'z'
    elif scheme in [b'AFT', b'CFT', b'WFT']:
        gSlow = 'r'
        gMid = 'theta'
        gFast = 'z'
    elif scheme in [b'SLFm', b'SLFl', b'BLF', b'WLF']:
        gSlow = 'r'
        gMid = 'theta'
        gFast = 'phi'
    nSlow = h5_file['mesh']['grid_'+gSlow].size
    nMid = h5_file['mesh']['grid_'+gMid].size
    nFast = h5_file['mesh']['grid_'+gFast].size
    time = h5_file['run']['time'].value
    sId = int(re.findall('\d+', inputfile.split('.')[0])[0])
    basename = inputfile.split(re.findall('\d+', inputfile.split('.')[0])[0])[0]
    # Set default output file
    if outputfile == '' and snapshots == 1:
        outputfile = basename+str(sId).zfill(4)+'.xdmf' 
    elif outputfile == '':
        outputfile = basename+'Series_'+str(sId).zfill(4)+'_'+str(sId+snapshots-1).zfill(4)+'.xdmf' 
    print("Input file: ", inputfile)
    print("Input ID: ", sId)
    print("Input base: ", basename)
    print('Output file: ', outputfile)
    print("Scheme: ", scheme)
    print(gSlow + "grid size: ", nSlow)
    print(gMid + " grid size: ", nMid)
    print(gFast + " grid size: ", nFast)
    print("Time: ", time)
    h5_file.close()

    # XDMF blocks and templates
    xdmfHead = '<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n<Xdmf Version="2.0">\n\t<Domain>'
    xdmfRevGrid = '\t\t\t\t<DataItem ItemType="Function" Function="-1.0*$0" Dimensions="{nD}">\n'
    xdmfRevGridEnd = '\t\t\t\t</DataItem>\n'
    xdmfVxVyVzGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DRectMesh" NumberOfElements="{nSlow} {nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="VxVyVz">\n{rFast}\t\t\t\t<DataItem Dimensions="{nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gFast}\n\t\t\t\t</DataItem>\n{rFastEnd}{rMid}\t\t\t\t<DataItem Dimensions="{nMid}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gMid}\n\t\t\t\t</DataItem>\n{rMidEnd}{rSlow}\t\t\t\t<DataItem Dimensions="{nSlow}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gSlow}\n\t\t\t\t</DataItem>\n{rSlowEnd}\t\t\t</Geometry>'
    xdmfXYZGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DSMesh" NumberOfElements="{nSlow} {nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="XYZ">\n\t\t\t\t<DataItem Dimensions="{nN} 3" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{gridfile}.hdf5:/mesh/grid_{gFast[0]}{gMid[0]}{gSlow[0]}\n\t\t\t\t</DataItem>\n\t\t\t</Geometry>'
    xdmfScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{sname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
    xdmfVScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{vname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
    xdmfTime = '\t\t\t<Time Value="{time}" />'
    xdmfEnd = '\t\t</Grid>\n\t</Domain>\n</Xdmf>'
    xdmfSeries = '\t\t<Grid Name="Timeseries" GridType="Collection" CollectionType="Temporal">'
    xdmfEndSeries = '\t\t</Grid>'

    # Open file
    out_file = open(outputfile, 'w')
    print(xdmfHead, file=out_file)
    if snapshots > 1:
        print(xdmfSeries, file=out_file)
    for fId in range(sId, sId+snapshots):
        current = basename+str(fId).zfill(4)+'.hdf5'
        h5_file = h5py.File(current, 'r')
        if scheme in [b'TTT', b'TFT',b'TFF',b'FFF']:
            rFast = False
            rMid = False
            rSlow = False
            if scheme in [b'TFF', b'TFT', b'TTT']:
                rSlow = True
            if scheme in [b'TFT', b'TTT']:
                rFast = True
            if scheme in [b'TTT']:
                rMid = True
            #if fId == sId:
            #    gridfunc = boxXYZ
            #    gridfname = 'box'
            #    makeGridFile(h5_file, gFast, gMid, gSlow, gridfname, gridfunc)
            #print(xdmfXYZGrid.format(nFast = nFast, nMid = nMid, nSlow = nSlow, nN = nFast*nMid*nSlow, gFast = gFast, gMid = gMid, gSlow = gSlow, gridfile = gridfname+'_grid'), file=out_file)

            print(xdmfVxVyVzGrid.format(nFast = nFast, nMid = nMid, nSlow = nSlow, fid = fId, basename = basename, gFast = gFast, gMid = gMid, gSlow = gSlow, rFast = xdmfRevGrid.format(nD = nFast) if rFast else '', rFastEnd = xdmfRevGridEnd if rFast else '', rMid = xdmfRevGrid.format(nD = nMid) if rMid else '', rMidEnd = xdmfRevGridEnd if rMid else '', rSlow = xdmfRevGrid.format(nD = nSlow) if rSlow else '', rSlowEnd = xdmfRevGridEnd if rSlow else ''), file=out_file)

        elif scheme in [b'CFT', b'AFT', b'BLF', b'SLFm', b'SLFl']:
            if fId == sId:
                if scheme in [b'CFT']:
                    gridfunc = cylinderXYZ
                    gridfname = 'cylinder'
                elif scheme in [b'AFT']:
                    gridfunc = annulusXYZ
                    gridfname = 'annulus'
                elif scheme in [b'BLF']:
                    gridfunc = sphereXYZ
                    gridfname = 'sphere'
                elif scheme in [b'SLFm', b'SLFl']:
                    gridfunc = shellXYZ
                    gridfname = 'shell'

                makeGridFile(h5_file, gFast, gMid, gSlow, gridfname, gridfunc)

            print(xdmfXYZGrid.format(nFast = nFast, nMid = nMid, nSlow = nSlow, nN = nFast*nMid*nSlow, gFast = gFast, gMid = gMid, gSlow = gSlow, gridfile = gridfname+'_grid'), file=out_file)

        # Create scalars  
        for s in list(h5_file):
            if s in list(h5_file[s]):
                print(xdmfScalar.format(nFast = nFast, nMid = nMid, nSlow = nSlow, sname = s, fid = fId, basename = basename), file=out_file)
        # Create vectors as scalars if requested
        if with_components:
            for v in list(h5_file):
                for ext in [gFast, gMid, gSlow]:
                    if v +  '_' + ext in list(h5_file[v]):
                        print(xdmfVScalar.format(nFast = nFast, nMid = nMid, nSlow = nSlow, vname = v, sname = v +  '_' + ext, fid = fId, basename = basename), file=out_file)
        time = h5_file['run']['time'].value
        print(xdmfTime.format(time = time), file=out_file)
        if snapshots > 1:
            print(xdmfEndSeries, file=out_file)
        h5_file.close()
    print(xdmfEnd, file=out_file)
    out_file.close()

def makeGridFile(h5_file, gFast, gMid, gSlow, fname, func):
    g_fast = h5_file['mesh']['grid_' + gFast]
    g_mid = h5_file['mesh']['grid_' + gMid]
    g_slow = h5_file['mesh']['grid_' + gSlow]
    size = g_fast.size*g_mid.size*g_slow.size
    grid_file = h5py.File(fname+'_grid.hdf5', 'w')
    if 'mesh' in grid_file and grid_file['mesh'].attrs['n_'+gFast[0]] == g_fast.size and grid_file['mesh'].attrs['n_'+gMid[0]] == g_mid.size and grid_file['mesh'].attrs['n_'+gSlow[0]] == g_slow.size:
        grid_file.close()
    else:
        if 'mesh' in grid_file:
            del grid_file['mesh']

        box = np.zeros([size,3])
        i = 0
        for ps in np.nditer(g_slow):
            for pm in np.nditer(g_mid):
                for pf in np.nditer(g_fast):
                    box[i,:] = func(pf, pm, ps)
                    i = i + 1
        mesh = grid_file.create_group('mesh')
        mesh.attrs['n_'+gFast[0]] = g_fast.size
        mesh.attrs['n_'+gMid[0]] = g_mid.size
        mesh.attrs['n_'+gSlow[0]] = g_slow.size
        dset = mesh.create_dataset('grid_'+gFast[0]+gMid[0]+gSlow[0], (size, 3), '=f8')
        dset[:,:] = box
        grid_file.close()

def boxXYZ(pFast, pMid, pSlow):
    return [pFast, pMid, pSlow]

def cylinderXYZ(pz, pth, pr):
    return [pz, pr*cos(pth), pr*sin(pth)]

def annulusXYZ(pz, pth, pr):
    return [pz, pr*cos(pth), pr*sin(pth)]

def sphereXYZ(pph, pth, pr):
    return [pr*cos(pth), pr*sin(pth)*cos(pph), pr*sin(pth)*sin(pph)]

def shellXYZ(pFast, pMid, pSlow):
    return [pr*cos(pth), pr*sin(pth)*cos(pph), pr*sin(pth)*sin(pph)]

if __name__ == "__main__":
    main(sys.argv[1:])
