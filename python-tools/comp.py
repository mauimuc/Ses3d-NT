#!/usr/bin/python
# -*- coding: utf-8 -*-

import obspy.core
import numpy as np
import matplotlib.pyplot as plt

path = './DATA/mess/'

st_new = obspy.core.read( './DATA/mess/Ses3D.X*.sac')

st_old = st_new.copy()

for tr in st_old:
    tr.data = (0,)

st_old.select( station='XX02', location='a', channel='vx')[0].data = np.loadtxt( path+'XX02_x', skiprows=7 )
st_old.select( station='XX02', location='a', channel='vy')[0].data = np.loadtxt( path+'XX02_y', skiprows=7 )
st_old.select( station='XX02', location='a', channel='vz')[0].data = np.loadtxt( path+'XX02_z', skiprows=7 )

st_old.select( station='XX02', location='b', channel='vx')[0].data = np.loadtxt( path+'XX22_x', skiprows=7 )
st_old.select( station='XX02', location='b', channel='vy')[0].data = np.loadtxt( path+'XX22_y', skiprows=7 )
st_old.select( station='XX02', location='b', channel='vz')[0].data = np.loadtxt( path+'XX22_z', skiprows=7 )

st_old.select( station='XX01', location='a', channel='vx')[0].data = np.loadtxt( path+'XX01_x', skiprows=7 )
st_old.select( station='XX01', location='a', channel='vy')[0].data = np.loadtxt( path+'XX01_y', skiprows=7 )
st_old.select( station='XX01', location='a', channel='vz')[0].data = np.loadtxt( path+'XX01_z', skiprows=7 )

st_old.select( station='XX01', location='b', channel='vx')[0].data = np.loadtxt( path+'XX21_x', skiprows=7 )
st_old.select( station='XX01', location='b', channel='vy')[0].data = np.loadtxt( path+'XX21_y', skiprows=7 )
st_old.select( station='XX01', location='b', channel='vz')[0].data = np.loadtxt( path+'XX21_z', skiprows=7 )

st_max = max(st_new.max() + st_old.max() )
st_min = min(st_new.max() + st_old.max() )
max = max( st_max, -st_min )

for loc, col in zip( ('a', 'b'), (1,2) ):
    for chan, row in zip( ('vx', 'vy', 'vz' ), (1,2,3) ):
        pos = row + col + (row-1) -1
        plt.subplot('%i%i%i' %( 3,2,pos))
        if ( col == 1 ):
            plt.ylabel( chan )
        if ( row == 1 ):
            plt.title( loc )
        for tr_new, tr_old in zip( st_new.select( channel=chan, location=loc ), \
                                   st_old.select( channel=chan, location=loc ) ):
            plt.plot( tr_new.data )
            plt.plot( tr_old.data )
            plt.ylim(-max, max)

plt.show()

