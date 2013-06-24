#!/usr/bin/env python

import sys, getopt

from math import *
import h5py
import numpy as np
import re

def main(argv):
   inputfile = ''
   outputfile = ''
   snapshots = 1 
   try:
      opts, args = getopt.getopt(argv,"hi:o:n:")
   except getopt.GetoptError:
      print 'Single file: createXDMF.py -i <inputfile> -o <outputfile>'
      print 'Timeseries: createXDMF.py -i <inputfile> -o <outputfile> -n <number of snapshots>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'Single file: createXDMF.py -i <inputfile> -o <outputfile>'
         print 'Timeseries: createXDMF.py -i <inputfile> -o <outputfile> -n <number of snapshots>'
         sys.exit()
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-o"):
         outputfile = arg
      elif opt in ("-n"):
         snapshots = int(arg)
   # Extract file information
   h5_file = h5py.File(inputfile, 'r')
   scheme = h5_file['/'].attrs['type']
   g1D = 'x'
   g2D = 'y'
   g3D = 'z'
   n1D = h5_file['mesh']['grid_'+g1D].size
   n2D = h5_file['mesh']['grid_'+g2D].size
   n3D = h5_file['mesh']['grid_'+g3D].size
   time = h5_file['run']['time'].value
   sId = int(re.findall('\d+', inputfile.split('.')[0])[0])
   basename = inputfile.split(re.findall('\d+', inputfile.split('.')[0])[0])[0]
   # Set default output file
   if outputfile == '' and snapshots == 1:
      outputfile = basename+str(sId).zfill(4)+'.xdmf' 
   elif outputfile == '':
      outputfile = basename+'Series_'+str(sId).zfill(4)+'_'+str(sId+snapshots-1).zfill(4)+'.xdmf' 
   print 'Input file: ', inputfile
   print 'Input ID: ', sId
   print 'Input base: ', basename
   print 'Output file: ', outputfile
   print 'Scheme: ', scheme
   print '1D grid size: ', n1D
   print '2D grid size: ', n2D
   print '3D grid size: ', n3D
   print 'Time: ', time
   h5_file.close()

   # Open file
   xdmfHead = '<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n<Xdmf Version="2.0">\n\t<Domain>'
   xdmfVxVyVzGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DRectMesh" NumberOfElements="%(n1D)u %(n2D)u %(n3D)u"/>\n\t\t\t<Geometry GeometryType="VxVyVz">\n\t\t\t\t<DataItem Dimensions="%(n3D)u" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t%(basename)s%(fid)04u.hdf5:/mesh/grid_%(g3D)s\n\t\t\t\t</DataItem>\n\t\t\t\t<DataItem Dimensions="%(n2D)u" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t%(basename)s%(fid)04u.hdf5:/mesh/grid_%(g2D)s\n\t\t\t\t</DataItem>\n\t\t\t\t<DataItem Dimensions="%(n1D)u" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t%(basename)s%(fid)04u.hdf5:/mesh/grid_%(g1D)s\n\t\t\t\t</DataItem>\n\t\t\t</Geometry>'
   xdmfXYZGrid = '\t\t<Grid Name="grid" GridType="Uniform">\n\t\t\t<Topology TopologyType="3DSMesh" NumberOfElements="%(n1D)u %(n2D)u %(n3D)u"/>\n\t\t\t<Geometry GeometryType="XYZ">\n\t\t\t\t<DataItem Dimensions="%(nN)u 3" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t%(gridfile)s.hdf5:/mesh/grid_%(g1D)s%(g2D)s%(g3D)s\n\t\t\t\t</DataItem>\n\t\t\t</Geometry>'
   xdmfScalar ='\t\t\t<Attribute Name="%(sname)s" AttributeType="Scalar" Center="Node">\n\t\t\t\t<DataItem Dimensions="%(n1D)u %(n2D)u %(n3D)u" NumberType="Float" Precision="8" Format="HDF">\n\t\t\t\t\t%(basename)s%(fid)04u.hdf5:/%(sname)s/%(sname)s\n\t\t\t\t</DataItem>\n\t\t\t</Attribute>'
   xdmfTime = '\t\t\t<Time Value="%(time)e" />'
   xdmfEnd = '\t\t</Grid>\n\t</Domain>\n</Xdmf>'
   xdmfSeries = '\t\t<Grid Name="Timeseries" GridType="Collection" CollectionType="Temporal">'
   xdmfEndSeries = '\t\t</Grid>'

   out_file = open(outputfile, 'w')
   print >>out_file, xdmfHead
   if snapshots > 1:
      print >>out_file, xdmfSeries
   for fId in range(sId, sId+snapshots):
      current = basename+str(fId).zfill(4)+'.hdf5'
      h5_file = h5py.File(current, 'r')
      if scheme in ('TFT'):
         #if fId == 1:
         #   boxXYZ(h5_file)
         #print >>out_file, xdmfXYZGrid % {'n1D': n1D, 'n2D': n2D, 'n3D': n3D, 'nN': n1D*n2D*n3D, 'g1D': g1D, 'g2D': g2D, 'g3D': g3D, 'gridfile': 'box_grid'}
         print >>out_file, xdmfVxVyVzGrid % {'n1D': n1D, 'n2D': n2D, 'n3D': n3D, 'fid': fId, 'basename': basename, 'g1D': g1D, 'g2D': g2D, 'g3D': g3D}
      elif scheme in ('CFT'):
         if fId == 1:
            annulusXYZ(h5_file)
            annulusWedgeXYZ(h5_file)
         print >>out_file, xdmfXYZGrid % {'n1D': n1D, 'n2D': n2D, 'n3D': n3D, 'nN': n1D*n2D*n3D, 'g1D': g1D, 'g2D': g2D, 'g3D': g3D, 'gridfile': 'annulus_grid'}
      elif scheme in ('SLF'):
         if fId == 1:
            shellXYZ(h5_file)
         print >>out_file, xdmfXYZGrid % {'n1D': n1D, 'n2D': n2D, 'n3D': n3D, 'nN': n1D*n2D*n3D, 'g1D': g1D, 'g2D': g2D, 'g3D': g3D, 'gridfile': 'shell_grid'}
      for s in list(h5_file):
         if s in list(h5_file[s]):
            print >>out_file, xdmfScalar % {'n1D': n1D, 'n2D': n2D, 'n3D': n3D, 'sname': s, 'fid': fId, 'basename': basename}
      time = h5_file['run']['time'].value
      print >>out_file, xdmfTime % {'time': time}
      if snapshots > 1:
         print >>out_file, xdmfEndSeries
      h5_file.close()
   print >>out_file, xdmfEnd
   out_file.close()

def boxXYZ(h5_file):
   g_x = h5_file['mesh']['grid_x']
   g_y = h5_file['mesh']['grid_y']
   g_z = h5_file['mesh']['grid_z']
   size = g_x.size*g_y.size*g_z.size
   box = np.zeros([size,3])
   i = 0
   for px in np.nditer(g_x):
      for py in np.nditer(g_y):
         for pz in np.nditer(g_z):
            box[i,:] = [pz, py, px]
            i = i + 1
   grid_file = h5py.File('box_grid.hdf5', 'w')
   mesh = grid_file.create_group('mesh')
   dset = mesh.create_dataset('grid_xyz', (size, 3), '=f8')
   dset[:,:] = box
   grid_file.close()

def annulusXYZ(h5_file):
   g_r = h5_file['mesh']['grid_x']
   g_th = h5_file['mesh']['grid_y']
   g_z = h5_file['mesh']['grid_z']
   size = g_r.size*g_th.size*g_z.size
   annulus = np.zeros([size,3])
   i = 0
   for pr in np.nditer(g_r):
      for pth in np.nditer(g_th):
         for pz in np.nditer(g_z):
            annulus[i,:] = [pz, pr*cos(pth), pr*sin(pth)]
            i = i + 1
   grid_file = h5py.File('annulus_grid.hdf5', 'w')
   mesh = grid_file.create_group('mesh')
   dset = mesh.create_dataset('grid_xyz', (size, 3), '=f8')
   dset[:,:] = annulus
   grid_file.close()

def annulusWedgeXYZ(h5_file):
   g_r = h5_file['mesh']['grid_x']
   g_th = h5_file['mesh']['grid_y']
   g_z = h5_file['mesh']['grid_z']
   size = g_r.size*2*g_z.size
   annulus = np.zeros([size,3])
   i = 0
   for pr in np.nditer(g_r):
      pth = g_th[-1]
      for pz in np.nditer(g_z):
         annulus[i,:] = [pz, pr*cos(pth), pr*sin(pth)]
         i = i + 1
      pth = g_th[0]
      for pz in np.nditer(g_z):
         annulus[i,:] = [pz, pr*cos(pth), pr*sin(pth)]
         i = i + 1
   grid_file = h5py.File('annulus_grid.hdf5','r+')
   dset = grid_file['mesh'].create_dataset('grid_wedge', (size, 3), '=f8')
   dset[:,:] = annulus
   grid_file.close()

def shellXYZ(h5_file):
   g_r = h5_file['mesh']['grid_x']
   g_th = h5_file['mesh']['grid_y']
   g_ph = h5_file['mesh']['grid_z']
   size = g_r.size*g_th.size*g_ph.size
   shell = np.zeros([size,3])
   i = 0
   for pr in np.nditer(g_r):
      for pth in np.nditer(g_th):
         for pph in np.nditer(g_ph):
            shell[i,:] = [pr*cos(pth), pr*sin(pth)*cos(pph), pr*sin(pth)*sin(pph)]
            i = i + 1
   grid_file = h5py.File('shell_grid.hdf5', 'w')
   mesh = grid_file.create_group('mesh')
   dset = mesh.create_dataset('grid_xyz', (size, 3), '=f8')
   dset[:,:] = shell
   grid_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
