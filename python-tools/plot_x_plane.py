#!/usr/bin/python
# -*- coding: utf-8 -*-

# ex_plot_x_plane.py -- v0.4
# Mo 7. Feb 21:16:00 CET 2011
# Stefan Mauerberger

###########################################################################
# this python script reads in and plots a SES3D data files in theta plain #
# parameters are set below                                                #
###########################################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import ses3dpy
from   matplotlib.transforms                 import Affine2D
from   matplotlib.projections                import PolarAxes
import mpl_toolkits.axisartist.floating_axes as     floating_axes
from   mpl_toolkits.axisartist.grid_finder   import FixedLocator

# read parameter file
par = ses3dpy.read_config_file( sys.argv[1] ) # read parameters of the simulation

### script parameters ####
path_f_x = 'DATA/mess/vz_'+sys.argv[2]+'.vol'       # file prefix for x-velocity field
src_t    = 90.0                         # depth to radius
gl_lon   = 5                            # lon grid lines
gl_r     = 5                            # lat grid lines

### read data
t, p, r = ses3dpy.generate_coordinates( par, pcolor=True )
p       = p * np.pi / 180                        # transform from grad to rad
R, P    = np.meshgrid( r, p )                     # generate 2d coordinate mesh
src_t   = np.argmin( np.abs( t - src_t ) )       # determine r-index for slice
f_x     = ses3dpy.read_field( par, path_f_x )    # read binary files
f_x     = f_x[ src_t, :, : ]                     # cut out slice at theta, phi

# set global parameters of matplotlib
plt.rcParams['axes.formatter.limits'] = (-4,4)  # powerlimits
plt.rcParams['text.usetex'] = False             # format text with TeX

# matplotlib formatter - for floating_axes plot
class MyFormatter(object):
    def __init__(self, fmt='$%f$',factor=1):
            self.fmt    = fmt
            self.factor = factor
    def __call__(self, direction, factor, values):
            values = values * self.factor
            return [self.fmt % v for v in values]
## figure
fig = plt.figure( figsize=(4.5,1.9), dpi=150 )
fig.subplots_adjust( bottom=0.11, top=1.0, right=0.95, left=0.08 )
# transform from polar koordinates to carthesien
tr =      Affine2D().scale( -np.pi/180., 1.)
tr = tr + Affine2D().translate( np.pi/2, 0)
tr = tr + PolarAxes.PolarTransform()
# grid
grid_locator_phi = FixedLocator( np.linspace( par['minX'][1], par['maxX'][1], 1 + gl_lon ) )
grid_locator_r   = FixedLocator( np.linspace( par['minX'][2], par['maxX'][2], 1 + gl_r ) )
# generate curved, floating axes
grid_helper = floating_axes.GridHelperCurveLinear(tr,
                         extremes=( par['minX'][1], par['maxX'][1],
                                    par['minX'][2], par['maxX'][2] ),\
                         grid_locator1   = grid_locator_phi,\
                         tick_formatter1 = MyFormatter('$%3.1f \, ^\circ$'),\
                         grid_locator2   = grid_locator_r ,\
                         tick_formatter2 = MyFormatter('$%i $',1e-3) )
# add axes to figurs
ax1 = floating_axes.FloatingSubplot( fig, 111, grid_helper=grid_helper )
fig.add_subplot( ax1 )
# labling options
ax1.axis['left'  ,'top'  ].toggle( all=True, label=True )
ax1.axis['bottom','right'].toggle( all=False, label=True )
ax1.axis['bottom'        ].set_label(r'$\phi \, [deg]$')
ax1.axis['right'         ].set_label(r'$r \, [km]$')
## plot wave field with colorbar
vmax = max( f_x.max(), -f_x.min() ) # max value for symmetric visualisation
x, y = R * np.sin(P), R * np.cos(P) # transform R, P to Cartesian coord.
cbar = ax1.pcolormesh( x, y, f_x, cmap = 'seismic', vmax=vmax, vmin=-vmax )
cbar.set_rasterized(True)
cbar = fig.colorbar( cbar, pad=-0.08, orientation='horizontal', shrink=0.6, aspect=40  )
cbar.set_label( r'$ v_\theta \, \mathrm{[m/s]}\!\times\!10^{-6}$' )  # label colorbar
ticks = ticks = np.linspace(-vmax,vmax,5)
cbar.set_ticks( ticks )
ticklabels=[]
for i in ticks:
    ticklabels.append( '$%.2f$' % ( i*1e6 ) )
cbar.set_ticklabels( ticklabels )
ax1.grid( True )                           # show grid

# export and/or display
#fig.savefig( '/home/mauerberger/thesis/images/ses3dpy/slice_x%s.eps' % nt , dpi=150 )
plt.show()

