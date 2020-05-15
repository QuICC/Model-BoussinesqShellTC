#!/usr/bin/env python

from __future__ import print_function

import sys, getopt

from math import *
import h5py
import numpy as np
import re

tab = ' '*2
endl = '\n'

def main(argv):
    inputfile = ''
    outputfile = ''
    snapshots = 1 
    with_components = False
    with_energy = False
    with_cylradius = False
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:", ['with-components', 'with-energy', 'with-cylradius'])
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
        elif opt in ("--with-cylradius"):
            with_cylradius = True
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
    elif scheme in [b'SLFm', b'SLFl', b'WLFl', b'WLFm']:
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
    time = h5_file['run']['time'][()]
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

    # Open file
    out_file = open(outputfile, 'w')
    if snapshots > 1:
        print(xdmfSeries(), file=out_file)
    for fId in range(sId, sId+snapshots):
        current = basename+str(fId).zfill(4)+'.hdf5'
        h5_file = h5py.File(current, 'r')
        gridfname = None
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

            print(xdmfHead((basename,fId), gridfname, nFast = nFast, nMid = nMid, nSlow = nSlow), file=out_file)
            print(xdmfVGrid([(nFast, gFast, rFast),(nMid, gMid, rMid),(nSlow, gSlow, rSlow)]), file=out_file)

        elif scheme in [b'CFT', b'WFT', b'AFT', b'SLFm', b'SLFl', b'WLFl', b'WLFm']:
            if fId == sId:
                efuncs = None
                if scheme in [b'CFT',b'WFT']:
                    gridfunc = cylinderXYZ
                    gridfname = 'cylinder'
                elif scheme in [b'AFT']:
                    gridfunc = annulusXYZ
                    gridfname = 'annulus'
                elif scheme in [b'WLFl',b'WLFm']:
                    gridfunc = {'n':genericGName,'h5':genericH5Grid,'g':sphereXYZ}
                    efuncs = [{'n':sphCThetaGName,'h5':sphH5CSThetaGrid,'f':np.cos},
                              {'n':sphSThetaGName,'h5':sphH5CSThetaGrid,'f':np.sin},
                              {'n':sphCPhiGName,'h5':sphH5CSPhiGrid,'f':np.cos},
                              {'n':sphSPhiGName,'h5':sphH5CSPhiGrid,'f':np.sin}]
                    gridfname = 'sphere'
                    vcompfunc = sphVComps
                elif scheme in [b'SLFm', b'SLFl']:
                    gridfunc = shellXYZ
                    gridfunc = {'n':genericGName,'h5':genericH5Grid,'g':shellXYZ}
                    efuncs = [{'n':sphCThetaGName,'h5':sphH5CSThetaGrid,'f':np.cos},
                              {'n':sphSThetaGName,'h5':sphH5CSThetaGrid,'f':np.sin},
                              {'n':sphCPhiGName,'h5':sphH5CSPhiGrid,'f':np.cos},
                              {'n':sphSPhiGName,'h5':sphH5CSPhiGrid,'f':np.sin}]
                    gridfname = 'shell'
                    vcompfunc = sphVComps

                makeGridFile(h5_file, gFast, gMid, gSlow, gridfname, gridfunc, efuncs)

            print(xdmfHead((basename,fId), gridfname, nFast, nMid, nSlow), file=out_file)
            print(xdmfGrid(comps = [gFast[0], gMid[0], gSlow[0]]), file=out_file)

        # Create scalars  
        for s in list(h5_file):
            if s in list(h5_file[s]):
                print(xdmfScalar(s), file=out_file)
        # Create vectors
        for v in list(h5_file):
            # Extract vectors
            if v +  '_' + gFast in h5_file[v] and len(h5_file[v]) == 3:
                vcomps = vcompfunc(v, gFast, gMid, gSlow)
                print(xdmfVector(vname = v, comps = vcomps), file=out_file)
        # Create vectors as scalars if requested
        if with_components:
            for v in list(h5_file):
                # Extract vectors
                for ext in [gFast, gMid, gSlow]:
                    if (ext is not None) and (v +  '_' + ext in list(h5_file[v])):
                        print(xdmfVScalar(vname = v, sname = v +  '_' + ext), file=out_file)
                # Extract 2D tensors
                for ext1 in [gFast, gMid, gSlow]:
                    for ext2 in [gFast, gMid, gSlow]:
                        if (ext1 is not None and ext2 is not None) and (v +  '_' + ext1 + ext2 in list(h5_file[v])):
                            print(xdmfVScalar(vname = v, sname = v +  '_' + ext1 + ext2), file=out_file)
        # Create energy density function if requested
        if with_energy:
            for v in list(h5_file):
                if v +  '_' + gFast in h5_file[v] and len(h5_file[v]) == 3:
                    print(xdmfEScalar(ename = v +  '_energy', vname = v, comps = [gSlow,gMid,gFast]), file=out_file)
        # Create cylindrical radial component if requested
        if with_cylradius:
            for v in list(h5_file):
                if v +  '_' + gFast in h5_file[v] and len(h5_file[v]) == 3:
                    comps = sphCylRadiusComps(v, gFast, gMid, gSlow)
                    print(xdmfFScalar(fname = v +  '_s', func = xdmfFuncCylRad, comps = comps), file=out_file)
        time = h5_file['run']['time'][()]
        print(xdmfTime(time), file=out_file)
        if snapshots > 1:
            print(xdmfEndSeries(), file=out_file)
        h5_file.close()
    print(xdmfEnd(), file=out_file)
    out_file.close()

def xdmfFuncEnergy(size):
    func = ""
    for i in range(0, size):
        func += f'${i}*${i}'
        if i < size - 1:
            func += ' + '
    return func

def xdmfFuncSph2X(size = 7):
    if size != 7:
        raise RuntimeError("Conversion from spherical to X component requires 7 fields!")
    func = '$0*$4*$5 + $1*$3*$5 - $2*$6'
    return func

def xdmfFuncSph2Y(size = 7):
    if size != 7:
        raise RuntimeError("Conversion from spherical to Y component requires 7 fields!")
    func = '$0*$4*$6 + $1*$3*$6 + $2*$5'
    return func

def xdmfFuncSph2Z(size = 4):
    if size != 4:
        raise RuntimeError("Conversion from spherical to Z component requires 4 fields!")
    func = '$0*$2 - $1*$3'
    return func

def xdmfFuncCylRad(size):
    if size != 4:
        raise RuntimeError("Cylindrical radius function expect 4 components!")
    func = '$0*$3 + $1*$2'
    return func

# XDMF blocks and templates
def xdmfHead(data, grid, nFast, nMid, nSlow):
    s  = '<?xml version="1.0" ?>' + endl
    s += '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [' + endl
    if nSlow is not None:
        s += tab + f'<!ENTITY geoType "XYZ">' + endl
        s += tab + f'<!ENTITY geoVType "VxVyVz">' + endl
        s += tab + f'<!ENTITY topoType "3DSMesh">' + endl
        s += tab + f'<!ENTITY topoVType "3DRectMesh">' + endl
        nN = nSlow*nMid*nFast
        s += tab + f'<!ENTITY sDimsGrid "{nSlow} {nMid} {nFast}">' + endl
        s += tab + f'<!ENTITY vDimsGrid "{nSlow} {nMid} {nFast} 3">' + endl
        s += tab + f'<!ENTITY gDimsGrid "{nN} 3">' + endl
    else:
        s += tab + f'<!ENTITY geoType "XY">' + endl
        s += tab + f'<!ENTITY geoVType "VxVy">' + endl
        s += tab + f'<!ENTITY topoType "2DSMesh">' + endl
        s += tab + f'<!ENTITY topoVType "2DRectMesh">' + endl
        nN = nMid*nFast
        s += tab + f'<!ENTITY sDimsGrid "{nMid} {nFast}">' + endl
        s += tab + f'<!ENTITY vDimsGrid "{nMid} {nFast} 2">' + endl
        s += tab + f'<!ENTITY gDimsGrid "{nN} 2">' + endl
    if data is not None:
        dataFile = f'{data[0]}{data[1]:04d}.hdf5'
        s += tab + f'<!ENTITY dataFile "{dataFile}">' + endl
    if grid is not None:
        gridFile = f'{grid}_grid.hdf5'
        s += tab + f'<!ENTITY gridFile "{gridFile}">' + endl
    s += ']>' + endl
    s += '<Xdmf Version="2.0">' + endl
    s += tab + '<Domain>'
    return s

def xdmfRevGrid(nD, bt = 4):
    s = tab*bt + f'<DataItem ItemType="Function" Function="-1.0*$0" Dimensions="{nD}">'
    return s

def xdmfRevGridEnd(bt = 4):
    s = tab*bt + '</DataItem>'
    return s

def xdmfVGrid(comps, bt = 2):
    s  = tab*bt + '<Grid Name="grid" GridType="Uniform">' + endl
    s += tab*(bt+1) + '<Topology TopologyType="&topoVType;" NumberOfElements="&sDimsGrid;"/>' + endl
    s += tab*(bt+1) + '<Geometry GeometryType="&geoVType;">' + endl
    for c in comps:
        if c[2]:
            s += xdmfRevGrid(c[0], bt = bt+2) + endl
        s += tab*(bt+2) + f'<DataItem Dimensions="{c[0]}" NumberType="Float" Precision="8" Format="HDF">' + endl
        s += tab*(bt+3) + '&dataFile;:/mesh/grid_{c[1]}' + endl
        s += tab*(bt+2) + '</DataItem>' + endl
        if c[2]:
            s += xdmfRevGridEnd(bt = bt+2) + endl
    s += tab*(bt+1) + '</Geometry>'
    return s

def xdmfGrid(comps, bt = 2):
    s  = tab*bt + '<Grid Name="grid" GridType="Uniform">' + endl
    s += tab*(bt+1) + '<Topology TopologyType="&topoType;" NumberOfElements="&sDimsGrid;"/>' + endl
    s += tab*(bt+1) + '<Geometry GeometryType="&geoType;">' + endl
    s += tab*(bt+2) + '<DataItem Dimensions="&gDimsGrid;" NumberType="Float" Precision="8" Format="HDF">' + endl
    s += tab*(bt+3) + '&gridFile;:/mesh/grid_' + "".join(comps) + endl
    s += tab*(bt+2) + '</DataItem>' + endl
    s += tab*(bt+1) + '</Geometry>'
    return s

def xdmfScalar(sname, bt = 3):
    s  = tab*bt + f'<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">' + endl
    s += tab*(bt+1) + '<DataItem Dimensions="&sDimsGrid;" NumberType="Float" Precision="8" Format="HDF">' + endl
    s += tab*(bt+2) + f'&dataFile;:/{sname}/{sname}' + endl
    s += tab*(bt+1) + '</DataItem>' + endl
    s += tab*bt + '</Attribute>'
    return s

def xdmfVScalar(vname, sname, bt = 3):
    s  = tab*bt + f'<Attribute Name="{sname}" AttributeType="Scalar" Center="Node">' + endl
    s += tab*(bt+1) + '<DataItem Dimensions="&sDimsGrid;" NumberType="Float" Precision="8" Format="HDF">' + endl
    s += tab*(bt+2) + f'&dataFile;:/{vname}/{sname}' + endl
    s += tab*(bt+1) + '</DataItem>' + endl
    s += tab*bt + '</Attribute>'
    return s

def xdmfFunction(func, comps, bt = 4):
    f = func(len(comps))
    s = tab*bt + f'<DataItem ItemType="Function" Function="{f}" Dimensions="&sDimsGrid;" NumberType="Float" Precision="8">' + endl
    for c in comps:
        s += tab*(bt+1) + '<DataItem Dimensions="&sDimsGrid;" NumberType="Float" Precision="8" Format="HDF">' + endl
        s += tab*(bt+2) + f'{c["f"]}:/{c["g"]}/{c["s"]}' + endl
        s += tab*(bt+1) + '</DataItem>' + endl
    s += tab*(bt) + '</DataItem>'
    return s

def xdmfFScalar(fname, func, comps, bt = 3):
    s  = tab*bt + f'<Attribute Name="{fname}" AttributeType="Scalar" Center="Node">' + endl
    s += xdmfFunction(func, comps, bt = bt+1) + endl
    s += tab*bt + '</Attribute>'
    return s

def xdmfVector(vname, comps, bt = 3):
    s =  tab*bt + f'<Attribute Name="{vname}" AttributeType="Vector" Center="Node">' + endl
    s += tab*(bt+1) + '<DataItem Dimensions="&vDimsGrid;" Function="JOIN($2, $0, $1)" ItemType="Function">' + endl
    for c in comps:
        s += xdmfFunction(c['func'], c['c'], bt = bt+2) + endl
    s += tab*(bt+1) + '</DataItem>' + endl
    s += tab*bt + '</Attribute>'
    return s

def xdmfEScalar(ename, vname, comps, bt = 3):
    s  = tab*bt + f'<Attribute Name="{ename}" AttributeType="Scalar" Center="Node">' + endl
    fcomps = list()
    for c in comps:
        fcomps.append({'f':'&dataFile;',
            'g':vname,
            's':f'{vname}_{c}'})
    s += xdmfFunction(xdmfFuncEnergy, fcomps, bt = bt+1) + endl
    s += tab*bt + '</Attribute>'
    return s

def xdmfTime(time, bt = 3):
    s = tab*bt + f'<Time Value="{time}" />'
    return s

def xdmfEnd(bt = 0):
    s  = tab*(bt+2) + '</Grid>' + endl
    s += tab*(bt+1) + '</Domain>' + endl
    s += tab*bt + '</Xdmf>'
    return s

def xdmfSeries(bt = 2):
    s = tab*bt + '<Grid Name="Timeseries" GridType="Collection" CollectionType="Temporal">'
    return s

def xdmfEndSeries(bt = 2):
    s = tab*bt + '</Grid>'
    return s

def makeGridFile(h5_file, gFast, gMid, gSlow, fname, gfunc, efuncs = None):
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
            g_dset = mesh.create_dataset(gfunc['n'](gFast, gMid, gSlow), (size, 3), '=f8')
            gfunc['h5'](g_dset, gfunc['g'], g_fast, g_mid, g_slow)

            if efuncs is not None:
                for efunc in efuncs:
                    e_dset = mesh.create_dataset(efunc['n'](gFast, gMid, gSlow), (g_slow.size, g_mid.size, g_fast.size), '=f8')
                    efunc['h5'](e_dset, efunc['f'], g_fast, g_mid, g_slow)

        else:
            g_dset = mesh.create_dataset(gfunc['n'](gFast, gMid), (size, 2), '=f8')
            gfunc['h5'](g_dset, gfunc['g'], g_fast, g_mid)

            if efuncs is not None:
                for efunc in efuncs:
                    e_dset = mesh.create_dataset(efunc['n'](gFast, gMid), (g_mid.size, g_fast.size), '=f8')
                    efunc['h5'](e_dset, efunc['f'], g_fast, g_mid)


        grid_file.close()
        
def boxXYZ(pFast, pMid, pSlow):
    return np.array([pFast, pMid*np.ones(pFast.shape), pSlow*np.ones(pFast.shape)])

def cylinderXYZ(pz, pth, pr):
    return np.array([pz, pr*np.cos(pth)*np.ones(pz.shape), pr*np.sin(pth)*np.ones(pz.shape)]).T

def annulusXYZ(pz, pth, pr):
    return np.array([pz, pr*np.cos(pth)*np.ones(pz.shape), pr*np.sin(pth)*np.ones(pz.shape)]).T

def sphereXYZ(pph, pth, pr):
    return np.array([pr*cos(pth)*np.ones(pph.shape), pr*np.sin(pth)*np.cos(pph), pr*np.sin(pth)*np.sin(pph)]).T

def shellXYZ(pph, pth, pr):
    return np.array([pr*cos(pth)*np.ones(pph.shape), pr*np.sin(pth)*np.cos(pph), pr*np.sin(pth)*np.sin(pph)]).T

def sphVComps(vname, gFast, gMid, gSlow):
    x = {'c': [{'f':'&dataFile;','g':vname,'s':f'{vname}_{gSlow}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gMid}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gFast}'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_p'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_p'},],
         'func':xdmfFuncSph2X}
    y = {'c': [{'f':'&dataFile;','g':vname,'s':f'{vname}_{gSlow}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gMid}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gFast}'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_p'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_p'},],
         'func':xdmfFuncSph2Y}
    z = {'c': [{'f':'&dataFile;','g':vname,'s':f'{vname}_{gSlow}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gMid}'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_t'},],
         'func':xdmfFuncSph2Z}
    return [x,y,z]

def sphCylRadiusComps(vname, gFast, gMid, gSlow):
    c = [{'f':'&dataFile;','g':vname,'s':f'{vname}_{gSlow}'},
               {'f':'&dataFile;','g':vname,'s':f'{vname}_{gMid}'},
               {'f':'&gridFile;','g':'mesh','s':'grid_cos_t'},
               {'f':'&gridFile;','g':'mesh','s':'grid_sin_t'},]
    return c

def sphH5CSRGrid(e_dset, e_func, g_fast, g_mid, g_slow):
    for k in range(0,g_slow.size):
        e_dset[k,:,:] = e_func(g_slow[k])

def sphH5CSThetaGrid(e_dset, e_func, g_fast, g_mid, g_slow):
    for j in range(0,g_mid.size):
        e_dset[:,j,:] = e_func(g_mid[j])

def sphH5CSPhiGrid(e_dset, e_func, g_fast, g_mid, g_slow):
    for i in range(0,g_fast.size):
       e_dset[:,:,i] = e_func(g_fast[i])

def genericGName(gFast, gMid, gSlow = None):
    if gSlow is not None:
        n = 'grid_'+gFast[0]+gMid[0]+gSlow[0]
    else:
        n = 'grid_'+gFast[0]+gMid[0]
    return n

def sphCThetaGName(gFast, gMid, gSlow = None):
    n = 'grid_cos_'+gMid[0]
    return n

def sphSThetaGName(gFast, gMid, gSlow = None):
    n = 'grid_sin_'+gMid[0]
    return n

def sphCPhiGName(gFast, gMid, gSlow = None):
    n = 'grid_cos_'+gFast[0]
    return n

def sphSPhiGName(gFast, gMid, gSlow = None):
    n = 'grid_sin_'+gFast[0]
    return n

def genericH5Grid(g_dset, g_func, g_fast, g_mid, g_slow):
    j = 0
    np_fast = np.array(g_fast[:])
    n_fast = g_fast.shape[0]
    ds_size = n_fast*g_mid.shape[0]
    box = np.zeros([ds_size,3])
    for ps in np.nditer(g_slow):
        i = 0
        for pm in np.nditer(g_mid):
            box[i:i+n_fast,:] = g_func(np_fast, pm, ps)
            i = i + n_fast
        g_dset[j:j+ds_size,:] = box
        j = j + ds_size

if __name__ == "__main__":
    main(sys.argv[1:])
