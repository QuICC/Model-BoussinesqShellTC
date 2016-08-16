import numpy as np
import numpy.fft as fft
import matplotlib.pylab as pl
import tables

nnx = 48
nny = 48
nnz = 64
Ek = 1e-7
Ro = Ek**(1./3.)
Ra = 20
Pr = 1
scale1d = 2.0

boxX = 10
boxY = 10

kc = 1.3048
kx = kc/np.sqrt(2.0)
ky = kc/np.sqrt(2.0)
tmp = np.genfromtxt("single_mode.dat", skip_header = 2)
nz  = tmp.shape[0]
grid = np.array([np.cos(i*np.pi/(nz - 1.0)) for i in range(0,nz)])
fields = ["mean_temperature", "dz_meantemperature", "w", "psi", "temperature"]
data = dict()
for i,field in enumerate(fields):
    data[field] = tmp[:,i]

for f in fields:
    phys = np.zeros(2*(nz - 1))
    phys[0:nz] = data[f][0:nz][::-1]
    if f == "mean_temperature":
        phys[0:nz] -= 1.0 - 0.5*(grid + 1.0)
        phys[0:nz] /= Ro
    phys[nz:] = phys[1:nz-1][::-1]
    data[f] = fft.rfft(phys).real/phys.shape[0]

hdf5_fields = ["mean_temperature", "velocity_tor", "velocity_pol", "temperature"]
hdf5_data = dict()
hdf5_data["mean_temperature"] = data["mean_temperature"]
hdf5_data["temperature"] = data["temperature"]
hdf5_data["velocity_tor"] = data["psi"]
hdf5_data["velocity_pol"] = (1.0/(kx**2 + ky**2))*data["w"]

hdf5_file = tables.open_file("single_state.hdf5", mode = 'w')
hdf5_file.set_node_attr('/', 'header', 'StateFile'.encode('ascii')) 
hdf5_file.set_node_attr('/', 'type', 'TFF'.encode('ascii')) 
hdf5_file.set_node_attr('/', 'version', '1.0'.encode('ascii')) 

#group = hdf5_file.create_group("/","mean_temperature")
#tmp = np.zeros((2*nny-1,nnx,nnz), dtype=np.complex128)
#tmp[0,0,0:nnz].real = hdf5_data['mean_temperature'][0:nnz]
#hdf5_file.create_array(group, "mean_temperature", tmp)

group = hdf5_file.create_group("/","physical")
hdf5_file.create_array(group, "bc_velocity", 1)
hdf5_file.create_array(group, "bc_temperature", 0)
hdf5_file.create_array(group, "ekman", Ek)
hdf5_file.create_array(group, "prandtl", Pr)
hdf5_file.create_array(group, "rayleigh", Ra)
hdf5_file.create_array(group, "scale1d", scale1d)

group = hdf5_file.create_group("/", "run")
hdf5_file.create_array(group, "time", 0.0)
hdf5_file.create_array(group, "timestep", 1e-4)

group = hdf5_file.create_group("/", "temperature")
tmp = np.zeros((2*nny-1,nnx,nnz), dtype=np.complex128)
tmp[0,0,0:nnz].real = hdf5_data['mean_temperature'][0:nnz]
tmp[boxY,0,0:nnz].real = 0.25*hdf5_data['temperature'][0:nnz]
tmp[0,boxX,0:nnz].real = 0.25*hdf5_data['temperature'][0:nnz]
hdf5_file.create_array(group, "temperature", tmp)

group = hdf5_file.create_group("/","truncation")
group2 = hdf5_file.create_group(group, "physical")
hdf5_file.create_array(group2, "dim1D", 1)
hdf5_file.create_array(group2, "dim2D", 1)
hdf5_file.create_array(group2, "dim3D", 1)
group2 = hdf5_file.create_group(group, "spectral")
hdf5_file.create_array(group2, "dim1D", nnz)
hdf5_file.create_array(group2, "dim2D", nnx)
hdf5_file.create_array(group2, "dim3D", nny)
group2 = hdf5_file.create_group(group, "transform")
hdf5_file.create_array(group2, "dim1D", 1)
hdf5_file.create_array(group2, "dim2D", 1)
hdf5_file.create_array(group2, "dim3D", 1)

group = hdf5_file.create_group("/","velocity")
tmp = np.zeros((2*nny-1,nnx,nnz), dtype=np.complex128)
tmp[boxY,0,0:nnz].real = 0.25*hdf5_data['velocity_tor'][0:nnz]
tmp[0,boxX,0:nnz].real = 0.25*hdf5_data['velocity_tor'][0:nnz]
hdf5_file.create_array(group, "velocity_tor", tmp)
tmp = np.zeros((2*nny-1,nnx,nnz), dtype=np.complex128)
tmp[boxY,0,0:nnz].real = 0.25*hdf5_data['velocity_pol'][0:nnz]
tmp[0,boxX,0:nnz].real = 0.25*hdf5_data['velocity_pol'][0:nnz]
hdf5_file.create_array(group, "velocity_pol", tmp)

hdf5_file.close()
