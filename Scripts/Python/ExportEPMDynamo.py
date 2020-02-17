import sys, getopt

import numpy as np
import h5py

input_file = ''
output_file = ''
argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print('ExportEPMDynamo.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('ExportEPMDynamo.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i"):
        input_file = arg
    elif opt in ("-o"):
        output_file = arg

# load QuICC state
quicc_file = h5py.File(input_file, 'r')
# Get truncation
spec1D = quicc_file['truncation']['spectral']['dim1D'][()]
spec2D = quicc_file['truncation']['spectral']['dim2D'][()]
spec3D = quicc_file['truncation']['spectral']['dim3D'][()]
# Get physical parameters
physE = quicc_file['physical']['ekman'][()]
physRa = quicc_file['physical']['rayleigh'][()]
physPr = quicc_file['physical']['prandtl'][()]
physPm = quicc_file['physical']['magnetic_prandtl'][()]
physQ = physPm/physPr
physRo = physE/physPm
print('QuICC physical parameters:')
print('\t'+'E = ' + str(physE))
print('\t'+'Q = ' + str(physQ))
print('\t'+'Ra = ' + str(physRa))
print('\t'+'Ro = ' + str(physRo))
# Get run parameters
runT = quicc_file['run']['time'][()]
runDt = quicc_file['run']['timestep'][()]

# create QuICC state file
epm_file = h5py.File(output_file,'w')
epm_file.attrs.create('StateFile', data=0, dtype='i4')
epm_file.attrs.create('Version', data=1, dtype='f4')

# Create run group
group = epm_file.create_group('RunParameters')
group.create_dataset('Time', (), 'f8', data = runT)
group.create_dataset('Step', (), 'f8', data = runDt)

# Create physical group
group = epm_file.create_group('PhysicalParameters')
group.create_dataset('E', (), 'f8', data = physE)
group.create_dataset('q', (), 'f8', data = physQ)
group.create_dataset('Ra', (), 'f8', data = physRa)
group.create_dataset('Ro', (), 'f8', data = physRo)

# Create truncation group
group = epm_file.create_group('Truncation')
group.create_dataset('N', (), 'i4', data = spec1D)
group.create_dataset('L', (), 'i4', data = spec2D)
group.create_dataset('M', (), 'i4', data = spec3D)
group.create_dataset('Mp', (), 'i4', data = 1)

def rescale(d, l, m):
    out = np.zeros(d.shape)
    lm = -1
    for k in range(0, l+1):
        norm0 = 1.0/ np.sqrt(4.0*np.pi/(2.0*k + 1.0))
        normm = 1.0/np.sqrt(8.0*np.pi/(2.0*k + 1.0))
        for j in range(0, k+1):
            lm += 1
            if j == 0:
                norm = norm0
            else:
                norm = normm
            out[lm,:,:] = norm*d[lm,:,:]
    return out

def add_background(d):
    print("### Adding thermal background state ###")
    d[0,0,0] += 1.110720734539592
    d[0,1,0] += -0.7853981633974483
    return d


# Create temperature field
group = epm_file.create_group('Codensity')
quicc_data = quicc_file['temperature']['temperature'][:]
ds = group.create_dataset('Codensity', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
ds[:] = rescale(add_background(quicc_data), spec2D, spec3D)

## Create velocity field
group = epm_file.create_group('Velocity')
quicc_data = quicc_file['velocity']['velocity_tor'][:]
ds = group.create_dataset('VelocityTor', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
ds[:] = rescale(quicc_data, spec2D, spec3D)
quicc_data = quicc_file['velocity']['velocity_pol'][:]
ds = group.create_dataset('VelocityPol', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
ds[:] = rescale(quicc_data, spec2D, spec3D)

## Create magnetic field
group = epm_file.create_group('Magnetic')
quicc_data = quicc_file['magnetic']['magnetic_tor'][:]
ds = group.create_dataset('MagneticTor', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
ds[:] = rescale(quicc_data, spec2D, spec3D)
quicc_data = quicc_file['magnetic']['magnetic_pol'][:]
ds = group.create_dataset('MagneticPol', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
ds[:] = rescale(quicc_data, spec2D, spec3D)

# Finished
epm_file.close()
quicc_file.close()
