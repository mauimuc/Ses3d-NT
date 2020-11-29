#! /usr/bin/python
# -*- coding: utf-8 -*-

import ses3dpy
import numpy as np
from sys import stdout


def gauss_3d( X, Y, Z, mx, my, mz, sx, sy, sz ):
    from numpy import exp
    r = ((X-mx)/(2.0*sx))**2 + ((Y-my)/(2.0*sy))**2 + ((Z-mz)/(2.0*sz))**2
    return exp( -r )

# Read configuration from file
config = ses3dpy.read_config_file( './l_aquila.mod' )
# Read model-parameter file
rhoinv = ses3dpy.read_field( config, './data/rhoinv')

# generate coordinates
t, p, r = ses3dpy.generate_coordinates( config, no_dup=False )
t = t * np.pi / 180.0
p = p * np.pi / 180.0
# 3d coordinate mesh
T,P,R = ses3dpy.mesh3d(t,p,r)
# transform to Cartesian coordinates
X = R * np.sin(T) * np.cos(P)
Y = R * np.sin(T) * np.sin(P)
Z = R * np.cos(T)


# Rank 3 field's shape
shape = rhoinv.shape

# Number of random seeds
seeds = 150

# Random coordinate indices in lateral direction
ix_rand = np.random.random_integers( 0, shape[0]-1, seeds )
# Random coordinate indices in longitudinal direction
iy_rand = np.random.random_integers( 0, shape[1]-1, seeds )
# Random coordinate indices in radial direction
iz_rand = np.random.random_integers( 0, shape[2]-1, seeds )

# Average grid spacing in lateral direction
dx_mean = (X[1:,:,:]-X[:-1,:,:]).mean()
# Average grid spacing in longitudinal direction
dy_mean = (Y[:,1:,:]-Y[:,:-1,:]).mean()
# Average grid spacing in radial direction
dz_mean = (Z[:,:,1:]-Z[:,:,:-1]).mean()
# With Gaussian blob
s = 10

# Allocate field to store random perturbations
pert = np.zeros( shape ); a=0
# Iterate through random coordinate indices
for ix, iy, iz in zip( ix_rand, iy_rand, iz_rand ):
    # Put some randomness in lateral standard deviation
    sx = s * dx_mean * (np.random.random() + 1.0)/2.0
    # Put some randomness in longitudinal standard deviation
    sy = s * dy_mean * (np.random.random() + 1.0)/2.0
    # Put some randomness in radial standard deviation
    sz = s * dz_mean * (np.random.random() + 1.0)/2.0
    # Add 'random' Gaussian blob to field of random perturbations
    # and scale it by a random number
    pert += gauss_3d( X=X, mx=X[ix,iy,iz], sx=sx, \
                      Y=Y, my=Y[ix,iy,iz], sy=sy, \
                      Z=Z, mz=Z[ix,iy,iz], sz=sz )*np.random.random()
    a = a+1
    print '\r%f%%' % (a*100.0/seeds),
    stdout.flush()

# Scale the whole field to the range of 0.0 to 1.0
pert = (pert - pert.min())/pert.ptp()

# Finally scale rhoinv by 10% of those random perturbations
rhoinv = rhoinv*(1.0 + 0.1*pert )

# Insert a discontinuity between first and second radial layer of elements 
# by taking the element's averages 
rhoinv[:,:,:5] = rhoinv[:,:,:5].mean()  

# Write field to disk
ses3dpy.write_field( config, rhoinv, './data/rhoinv_pert' )


