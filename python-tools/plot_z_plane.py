#!/usr/bin/python
# -*- coding: utf-8 -*-

# plot_z_plane.py -- v0.4
# So 13. Feb 11:57:26 CET 2011
# by Stefan Mauerberger

#####################################################################
# this script reads in and plots a SES3D data files with r constant #
# all plotting parameters are to be set below                       #
#####################################################################

import SES3dPy
import sys
import numpy                as np
import matplotlib.pyplot    as plt
from   mpl_toolkits.basemap import Basemap

# read parameter file
par = SES3dPy.read_config_file( sys.argv[1] ) # read parameters of the simulation
plt.rcParams['text.usetex'] = False     # format text with TeX

### script parameters ####
path_f_x = 'DATA/mess/vz_'+str(sys.argv[2])+'.vol'      # file prefix for x-velocity field
src_r    = 6371000.0                   # depth to radius
gl_lon   = 5                           # lon grid lines
gl_lat   = 5                           # lat grid lines


### read data
t, p, r  = SES3dPy.generate_coordinates( par, pcolor=True )
lat, lon = np.meshgrid( 90. - t, p )            # generate 2d coordinate grid
src_r    = np.argmin( np.abs( r - src_r ) )-1   # determine r-index for slice
f_x      = SES3dPy.read_field( par, path_f_x )  # read binary files
f_x      = f_x[ :, :, src_r ]                   # cut out slice at theta, phi

# set global parameters of matplotlib
plt.rcParams['axes.formatter.limits'] = (-4,4)  # power-limits
## figure
fig = plt.figure( figsize=(4.5 , 2.8)  )
fig.subplots_adjust( top=0.99, bottom=0.0, left=0.08, right=0.98 )
ax = fig.add_subplot( (111) )

## plot base-map & projection
# projection parameters for lamberts conformal conic
height = r[src_r] * np.pi * ( par['maxX'][0] - par['minX'][0] ) / 180.
width  = r[src_r] * np.pi * ( par['maxX'][1] - par['minX'][1] ) / 180. \
         * np.sin( np.pi * par['maxX'][0] / 180. )
lat_0  = 90.0 - ( par['maxX'][0] + par['minX'][0] ) / 2
lat_1  = 90.1 - ( par['maxX'][0] )
lat_2  = 90.0 - ( par['minX'][0] )
lon_0  = par['minX'][1] + ( par['maxX'][1] - par['minX'][1] ) / 2
m = Basemap( resolution='l', projection='lcc', ax=ax, \
             width=width, height=height, rsphere=r[src_r], \
             lat_0=lat_0, lat_1=lat_1, lat_2=lat_2, lon_0=lon_0 )
gl_lat = 90.0 - np.linspace( par['maxX'][0], par['minX'][0], gl_lat + 1 )[::-1]
gl_lon =        np.linspace( par['maxX'][1], par['minX'][1], gl_lon + 1 )
m.drawcoastlines()
ax.yaxis.set_label_position('right')
ax.yaxis.set_label_text( r'$\theta \, [deg]$' )
m.drawparallels( gl_lat, labels=[1,0,0,0], dashes=[2,2],\
                 labelstyle='+/-', linewidth=0.5, fmt='%.0f')
ax.xaxis.set_label_position('top')
ax.xaxis.set_label_text( r'$\phi \, [deg]$' )
m.drawmeridians( gl_lon, labels=[0,0,0,1], dashes=[2,2], \
                 labelstyle='+/-', linewidth=0.5, fmt='%.0f')
## plot wave field with colorbar
vmax = max( f_x.max(), -f_x.min() ) # max value for symmetric visualisation
x, y = m( lon, lat )                # transform lat, lon to Cartesian coord.
a = m.pcolormesh( x, y, f_x.T, vmin=-vmax, vmax=vmax, cmap='seismic')
a.set_rasterized( True )                 # stores wave field as raster grafix
cbar = fig.colorbar( a, shrink=0.6, aspect=30,pad=0.09 ) # generates a color bar
## labelling the color bar
cbar.set_label( r'$ v_\theta \, [m/s]$' )

# export and/or display
#fig.savefig( '/home/mauerberger/thesis/images/ses3dpy/slice_z%s.eps' % nt , dpi=150 )
plt.show()

