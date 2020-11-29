#!/usr/bin/python
# -*- coding: utf-8 -*-

# plot_mode_z_1d.py -- v0.4.2
# Sat Jan  5 10:26:31 CET 2013
# by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>

########################################################################
# this python script reads in and plots SES3Ds radial model parameters #
# all plotting parameters are to be set below                          #
########################################################################

import sys
import ses3dpy
import numpy as np
import matplotlib.pyplot as plt

# set global parameters of matplotlib
plt.rcParams['axes.formatter.limits'] = (-4,4)  # powerlimits
plt.rcParams['text.usetex'] = False             # format text with TeX

# parse parameter file
par = ses3dpy.read_config_file( sys.argv[1] )

## read data
t, p, r   = ses3dpy.generate_coordinates( par, no_dup=False )
# reads data (no_dup because of potentially discontinuities)
rho = 1. / ses3dpy.read_field( par, par['rhoinv']) # read binary files
rho = rho[ 1, 1, : ] # cut out slice at theta, phi
lam = ses3dpy.read_field( par, par['lambda']) # read binary files
lam = lam[ 1, 1, : ] # cut out slice at theta, phi
mu  = ses3dpy.read_field( par, par['mu']) # read binary files
mu  =  mu[ 1, 1, : ] # cut out slice at theta, phi

## plot model
fig  = plt.figure( figsize=(10,6))
fig.subplots_adjust(right=0.98,left=0.1,top=0.92,bottom=0.07)
ax   = fig.add_subplot( 111 )
ax.set_ylabel( r'$r \, [m]$' )
#ax.set_ylim( r[1] * 0.999, r[-1] * 1.001 )
ax.set_yticks( np.linspace(par['minX'][2],par['maxX'][2],5) )
ax.grid()

# plot rho
ax.plot( rho * 1e-3, r , label=r'$\rho(r) \,10^{3} [kg/m^3]$', color='b' )
ax.plot( rho * 1e-3, r , '.', color='b' ) # optional - adds dots on each point
# plot mu
ax.plot( mu * 1e-11 , r , label=r'$\mu(r) \, 10^{11} [Nm]$', color='r' )
ax.plot( mu * 1e-11, r , '.', color='r' )  # optional - adds dots on each point
# plot lambda
ax.plot( lam * 1e-11 , r , label=r'$\lambda(r) \, 10^{11} [Nm]$', color='g' )
ax.plot( lam * 1e-11, r , '.', color='g' ) # optional - adds dots on each point

## plot c_p
#v_p = np.sqrt( (lam+2*mu)/rho )
#ax.plot( v_p*1e-3, r, '--', label=r'$v_p(r) \, [km/s]$', color='k' )
##ax.plot( rho , r , '.', color='b' ) # optional - adds dots on each point
## plot c_s
#v_s =  np.sqrt( (mu)/rho )
#ax.plot( v_s*1e-3, r, '--',  label=r'$v_s(r) \, [km/s]$', color='k' )
##ax.plot( mu, r , '.', color='r' )  # optional - adds dots on each point

# add legend
leg = ax.legend( loc=8, borderaxespad=0 )
leg.draw_frame( False )

# export and/or display
#fig.savefig( './model_z.pdf' , dpi=150 )
plt.show()
