#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset

# physical constants
D   = 4000   # depth of water

# grid dimensions
nx = 400
ny = 250
nz = 72

# grid resolution (10km)
rx = 1.e4
dx = np.ones(nx)*rx
dy = np.ones(ny)*rx

Lx = nx*rx
Ly = ny*rx

x = np.cumsum(dx) - rx/2.0 -Lx/2.0
y = np.cumsum(dy) - rx/2.0 -Ly/2.0

X,  Y  = np.meshgrid(x, y)

# create vertical grid
dzmin = 5
dzmax = 98
decay = 512.3207
z = np.empty(nz)
zw = np.empty(nz+1)
dz = np.empty(nz)
z[0] = -dzmin / 2
zw[0] = 0
for i in range(0, nz - 1):
    dz[i] = dzmax - (dzmax - dzmin) / np.cosh(z[i] / decay)**2
    z[i+1] = z[i] - dz[i]
    zw[i+1] = -(z[i] + z[i+1]) / 2
dz[nz - 1] = D - dz.sum()
zw[nz] = D


# wind stress 
tau_x = 0.2*(np.sin(np.pi*((Y+Ly/2.0)/Ly)))

# SST relaxation - surface restoring
SST1 = np.where(y>y[-1]-5.e5)

ind_sst_nor = SST1[0][0]
ind_sst_sou = 0

SST = np.full((ny,nx),np.nan)
SST[ind_sst_nor:,:] = 12.

for j in range(ind_sst_sou,ind_sst_nor):
    SST[j,:] = 12*(y[j]-y[ind_sst_sou])/(y[ind_sst_nor]-y[ind_sst_sou])

# sponge layer parameter 1
# Idamp: the inverse damping rate 
Idamp = np.zeros((ny,nx))
# width: 100km 
# unit of the inverse damping rate:s-1
index = np.where(y>y[-1]-1.e5)
ind = index[0][0]
num = ny - ind
# time scales !!!! not the rate
mdamp = 7*24*3600   # north
for m in np.arange(ind,ny):
    Idamp[m,:] = np.linspace(0.,1./mdamp,num)[m-ind]
### Post processing of surface restoring above sponge
    SST[m,:] = SST[ind,0]

# sponge layer parameter 2
# ETA: the interface height
ETA = np.nan*np.ones((nz+1,ny,nx))
for k in range(0,nz+1):
    ETA[k,:,:] = -zw[k]

# sponge layer parameter 3
# PTEMP: the potential temperature profile
PTEMP = np.nan*np.ones((nz,ny,nx))

# PTEMP e-folding depth = ~500m
index2 = np.where(-z<1000)
ind2 = np.max(index2)

# PTEMP should only vary with depth for a specific y value, be non-zero in the sponge layer

for j in range(ny):
    # PTEMP at the surface should match the surface restoring temperature SST
    PTEMP[0,j,:] = SST[ny-1,0]
    for k in range(1,nz):
        #PTEMP[k,j,:] = (PTEMP[0,j,:]-2.0)*np.exp(-((-z[k])/1000))+2.0
        PTEMP[k,j,:] = (PTEMP[0,j,:]-2.0)*(np.exp(z[k]/500)-np.exp(-D/500))/(1-np.exp(-D/500))+2.0

# sponge layer parameter 4
# SALT: the salinity 
SALT = 35.0*np.ones((nz,ny,nx))

# coord.nc - initial salinity
isali = 35.0*np.ones((nz,ny,nx))

# coord.nc - initial temperature
itemp = np.ones(isali.shape)

for j in range(nz):
    itemp[j,:,:] = PTEMP[j,ind+1,0]

# Make base topography
h = -D*np.ones(X.shape)
# First, add a ridge-like feature.
h = h+(2500 + 300*np.sin(10*np.pi*Y/Ly) + 400*np.sin(8*np.pi*Y/Ly)+ 300*np.sin(25*np.pi*Y/Ly) )/np.cosh((X-0.2*Y+3e5)/1.2e5);
# Now add another ridge!!
h = h+(2000 + 600*np.sin(11*np.pi*Y/Ly) + 300*np.sin(7*np.pi*Y/Ly)+ 500*np.sin(21*np.pi*Y/Ly) )/np.cosh((X+0.1*Y+1.5e6)/1.2e5);
# Next add a continental slope next to Antarctica
h = np.maximum(h,-0.01*(Y+Ly/2))

# Unify meridional depth within the sponge layer
for p in np.arange(ind,ny):
    for q in range(nx):
        h[p,q] = h[ind,q]
       
# write out data
with Dataset('/short/v45/lxy581/mom6/input/exp1_w2/indata.nc', 'w') as d:
    d.createDimension('x', nx)
    d.createDimension('y', ny)
    d.createDimension('z', nz)
    d.createDimension('zt', nz)
    d.createDimension('zw', nz+1)

    d.createVariable('x', 'f8', ('x',))[:] = x
    d.createVariable('y', 'f8', ('y',))[:] = y
    d.createVariable('z', 'f8', ('z',))[:] = z
    d.createVariable('zt', 'f8', ('zt',))[:] = -z
    d.createVariable('zw', 'f8', ('zw',))[:] = zw
    d.createVariable('dz', 'f8', ('z',))[:] = dz
    d.createVariable('topog', 'f8', ('y', 'x'))[:] = -h
    d.createVariable('taux',  'f8', ('y', 'x'))[:] = tau_x
    d.createVariable('sst',   'f8', ('y', 'x'))[:] = SST
    d.createVariable('Idamp', 'f8', ('y', 'x'))[:] = Idamp
    d.createVariable('ETA', 'f8', ('zw','y', 'x'))[:] = ETA
    d.createVariable('PTEMP', 'f8', ('zt','y', 'x'))[:] = PTEMP
    d.createVariable('SALT', 'f8', ('zt','y', 'x'))[:] = SALT
    d.createVariable('isali', 'f8', ('zt','y', 'x'),fill_value=-1e+20)[:] = isali
    d.createVariable('itemp', 'f8', ('zt','y', 'x'),fill_value=-1e+20)[:] = itemp

