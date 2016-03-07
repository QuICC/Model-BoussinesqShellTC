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
    with_energy = False
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:", ['with-components', 'with-energy'])
    except getopt.GetoptError:
        print('Single file: createXDMF.py (--with-components) (--with-energy) -i <inputfile> -o <outputfile>')
        print('Timeseries: createXDMF.py (--with-components) (--with-energy) -i <inputfile> -o <outputfile> -n <number of snapshots>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Single file: createXDMF.py (--with-components) (--with-energy) -i <inputfile> -o <outputfile>')
            print('Timeseries: createXDMF.py (--with-components) (--with-energy) -i <inputfile> -o <outputfile> -n <number of snapshots>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("--with-components"):
            with_components = True
        elif opt in ("--with-energy"):
            with_energy = True
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
    elif scheme in [b'SLFm', b'SLFl', b'BLFl', b'BLFm', b'WLF']:
        gSlow = 'r'
        gMid = 'theta'
        gFast = 'phi'
    elif scheme in [b'TT']:
        gSlow = None
        gMid = 'x'
        gFast = 'z'
    elif scheme in [b'TF']:
        gSlow = None
        gMid = 'z'
        gFast = 'x'
    if gSlow is not None:
        nSlow = h5_file['mesh']['grid_'+gSlow].size
    else:
        nSlow = 1
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
    if gSlow is not None:
        print(gSlow + " grid size: ", nSlow)
    print(gMid + " grid size: ", nMid)
    print(gFast + " grid size: ", nFast)
    print("Time: ", time)
    h5_file.close()

    # XDMF blocks and templates
    xdmfHead = '<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n<Xdmf Version="2.0">\n\t<Domain>'
    xdmfRevGrid = '\t\t\t\t<DataItem ItemType="Function" Function="-1.0*$0" Dimensions="{nD}">\n'
    xdmfRevGridEnd = '\t\t\t\t</DataItem>\n'
    if gSlow is not None:
        xdmfVxVyVzGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DRectMesh" NumberOfElements="{nSlow} {nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="VxVyVz">\n{rFast}\t\t\t\t<DataItem Dimensions="{nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gFast}\n\t\t\t\t</DataItem>\n{rFastEnd}{rMid}\t\t\t\t<DataItem Dimensions="{nMid}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gMid}\n\t\t\t\t</DataItem>\n{rMidEnd}{rSlow}\t\t\t\t<DataItem Dimensions="{nSlow}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gSlow}\n\t\t\t\t</DataItem>\n{rSlowEnd}\t\t\t</Geometry>'
        xdmfXYZGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DSMesh" NumberOfElements="{nSlow} {nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="XYZ">\n\t\t\t\t<DataItem Dimensions="{nN} 3" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{gridfile}.hdf5:/mesh/grid_{gFast[0]}{gMid[0]}{gSlow[0]}\n\t\t\t\t</DataItem>\n\t\t\t</Geometry>'
        xdmfScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{sname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
        xdmfVScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{vname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
        xdmfEFunction ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem ItemType="Function" Function="$0*$0 + $1*$1 + $2*$2" Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8">'
        xdmfEScalar ='\t\t\t\t\t<DataItem Dimensions="{nSlow} {nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t\t{basename}{fid:04d}.hdf5:/{vname}/{sname}\n\t\t\t\t\t</DataItem>'
        xdmfEEnd ='\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
    else:
        xdmfVxVyVzGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="2DRectMesh" NumberOfElements="{nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="VxVy">\n{rMid}\t\t\t\t<DataItem Dimensions="{nMid}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gMid}\n\t\t\t\t</DataItem>\n{rMidEnd}{rFast}\t\t\t\t<DataItem Dimensions="{nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/mesh/grid_{gFast}\n\t\t\t\t</DataItem>\n{rFastEnd}\t\t\t</Geometry>'
        xdmfXYZGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="2DSMesh" NumberOfElements="{nMid} {nFast}"/>\n\t\t\t<Geometry GeometryType="XY">\n\t\t\t\t<DataItem Dimensions="{nN} 2" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{gridfile}.hdf5:/mesh/grid_{gFast[0]}{gMid[0]}\n\t\t\t\t</DataItem>\n\t\t\t</Geometry>'
        xdmfScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{sname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
        xdmfVScalar ='\t\t\t<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="{nMid} {nFast}" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t{basename}{fid:04d}.hdf5:/{vname}/{sname}\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
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
        if scheme in [b'TTT', b'TFT',b'TFF',b'FFF',b'TT',b'TF']:
            rFast = False
            rMid = False
            rSlow = False
            if scheme in [b'TFF', b'TFT', b'TTT']:
                rSlow = True
            if scheme in [b'TFT', b'TTT', b'TT']:
                rFast = True
            if scheme in [b'TTT', b'TT', b'TF']:
                rMid = True
            #if fId == sId:
            #    gridfunc = boxXYZ
            #    gridfname = 'box'
            #    makeGridFile(h5_file, gFast, gMid, gSlow, gridfname, gridfunc)
            #print(xdmfXYZGrid.format(nFast = nFast, nMid = nMid, nSlow = nSlow, nN = nFast*nMid*nSlow, gFast = gFast, gMid = gMid, gSlow = gSlow, gridfile = gridfname+'_grid'), file=out_file)

            print(xdmfVxVyVzGrid.format(nFast = nFast, nMid = nMid, nSlow = nSlow, fid = fId, basename = basename, gFast = gFast, gMid = gMid, gSlow = gSlow, rFast = xdmfRevGrid.format(nD = nFast) if rFast else '', rFastEnd = xdmfRevGridEnd if rFast else '', rMid = xdmfRevGrid.format(nD = nMid) if rMid else '', rMidEnd = xdmfRevGridEnd if rMid else '', rSlow = xdmfRevGrid.format(nD = nSlow) if rSlow else '', rSlowEnd = xdmfRevGridEnd if rSlow else ''), file=out_file)

        elif scheme in [b'CFT', b'AFT', b'BLFl', b'BLFm', b'SLFm', b'SLFl']:
            if fId == sId:
                if scheme in [b'CFT']:
                    gridfunc = cylinderXYZ
                    gridfname = 'cylinder'
                elif scheme in [b'AFT']:
                    gridfunc = annulusXYZ
                    gridfname = 'annulus'
                elif scheme in [b'BLFl',b'BLFm']:
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
                # Extract vectors
                for ext in [gFast, gMid, gSlow]:
                    if (ext is not None) and (v +  '_' + ext in list(h5_file[v])):
                        print(xdmfVScalar.format(nFast = nFast, nMid = nMid, nSlow = nSlow, vname = v, sname = v +  '_' + ext, fid = fId, basename = basename), file=out_file)
                # Extract 2D tensors
                for ext1 in [gFast, gMid, gSlow]:
                    for ext2 in [gFast, gMid, gSlow]:
                        if (ext1 is not None and ext2 is not None) and (v +  '_' + ext1 + ext2 in list(h5_file[v])):
                            print(xdmfVScalar.format(nFast = nFast, nMid = nMid, nSlow = nSlow, vname = v, sname = v +  '_' + ext1 + ext2, fid = fId, basename = basename), file=out_file)
        # Create energy density function if requested
        if with_energy:
            for v in list(h5_file):
                if v +  '_' + gFast in h5_file[v] and len(h5_file[v]) == 3:
                    print(xdmfEFunction.format(nFast = nFast, nMid = nMid, nSlow = nSlow, sname = v +  '_energy'), file=out_file)
                    for ext in [gFast, gMid, gSlow]:
                        if (ext is not None) and (v +  '_' + ext in list(h5_file[v])):
                            print(xdmfEScalar.format(nFast = nFast, nMid = nMid, nSlow = nSlow, vname = v, sname = v +  '_' + ext, fid = fId, basename = basename), file=out_file)
                    print(xdmfEEnd, file=out_file)
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
    if gSlow is not None:
        g_slow = h5_file['mesh']['grid_' + gSlow]
        size = g_fast.size*g_mid.size*g_slow.size
    else:
        size = g_fast.size*g_mid.size
    grid_file = h5py.File(fname+'_grid.hdf5', 'w')
    if 'mesh' in grid_file and grid_file['mesh'].attrs['n_'+gFast[0]] == g_fast.size and grid_file['mesh'].attrs['n_'+gMid[0]] == g_mid.size and grid_file['mesh'].attrs['n_'+gSlow[0]] == g_slow.size:
        grid_file.close()
    else:
        if 'mesh' in grid_file:
            del grid_file['mesh']
        
        # Create mesh dataset
        mesh = grid_file.create_group('mesh')
        mesh.attrs['n_'+gFast[0]] = g_fast.size
        mesh.attrs['n_'+gMid[0]] = g_mid.size
        if gSlow is not None:
            mesh.attrs['n_'+gSlow[0]] = g_slow.size
            dset = mesh.create_dataset('grid_'+gFast[0]+gMid[0]+gSlow[0], (size, 3), '=f8')
        else:
            dset = mesh.create_dataset('grid_'+gFast[0]+gMid[0], (size, 2), '=f8')

        j = 0
        np_fast = np.array(g_fast[:])
        n_fast = g_fast.shape[0]
        ds_size = n_fast*g_mid.shape[0]
        box = np.zeros([ds_size,3])
        for ps in np.nditer(g_slow):
            i = 0
            for pm in np.nditer(g_mid):
                box[i:i+n_fast,:] = func(np_fast, pm, ps)
                i = i + n_fast
            dset[j:j+ds_size,:] = box
            j = j + ds_size
        grid_file.close()

def boxXYZ(pFast, pMid, pSlow):
    return np.array([pFast, pMid*np.ones(pFast.shape), pSlow*np.ones(pFast.shape)])

def cylinderXYZ(pz, pth, pr):
    return np.array([pz, pr*np.cos(pth)*np.ones(pz.shape), pr*np.sin(pth)*np.ones(pz.shape)])

def annulusXYZ(pz, pth, pr):
    return np.array([pz, pr*np.cos(pth)*np.ones(pz.shape), pr*np.sin(pth)*np.ones(pz.shape)])

def sphereXYZ(pph, pth, pr):
    return np.array([pr*cos(pth)*np.ones(pph.shape), pr*np.sin(pth)*np.cos(pph), pr*np.sin(pth)*np.sin(pph)]).T

def shellXYZ(pph, pth, pr):
    return np.array([pr*cos(pth)*np.ones(pph.shape), pr*np.sin(pth)*np.cos(pph), pr*np.sin(pth)*np.sin(pph)]).T

if __name__ == "__main__":
    main(sys.argv[1:])
