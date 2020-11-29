#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import Beach
import ses3dpy
import namelist
import sys

c_file = sys.argv[1]

config = ses3dpy.read_config_file(c_file)

plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)

llcrnrlon, urcrnrlat = ( config['minX'][1], 90.0 - config['minX'][0] )
urcrnrlon, llcrnrlat = ( config['maxX'][1], 90.0 - config['maxX'][0] )
lon_0    , lat_0     = ( (llcrnrlon+urcrnrlon)/2.0, (llcrnrlat+urcrnrlat)/2.0 ) 

m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, \
            resolution='i', projection='merc', \
            lon_0=lon_0, lat_0=lat_0, lat_ts=lat_0)


m.drawcoastlines()
m.fillcontinents()

# Draw a computational domain
m.drawparallels( [llcrnrlat, urcrnrlat] )
m.drawmeridians( [llcrnrlon, urcrnrlon] )


# Print receivers
receivers = namelist.namelist2dict(c_file)['receiver']

lats = []
lons = []
names = []

for rec in receivers:
    lon = rec['lon']
    lat = rec['lat']
    name = rec['station']
    if ( lon in lons and lat in lats ):
        pass
    else:
        lats.append(lat)
        lons.append(lon)
        names.append(name)

x, y = m(lons, lats)
m.scatter(x, y, 150, color="b", marker=".", edgecolor="b", zorder=3)
for x_, y_, name_ in zip(x,y,names):
    plt.text(x_, y_, name_)#, va="top", family="monospace", weight="bold")


# Tapering boundaries 

# Width of an element 
dlon = ( urcrnrlon - llcrnrlon )/config['nX_global'][1]
dlat = ( urcrnrlat - llcrnrlat )/config['nX_global'][0]
# Taper width
dlon *= config['pml']
dlat *= config['pml']

# Southern taper 
lons =  np.linspace(llcrnrlon+dlon, urcrnrlon-dlon, 100).tolist()
lats =  np.linspace(llcrnrlat+dlat, llcrnrlat+dlat, 100).tolist()
# Eastern taper
lons += np.linspace(urcrnrlon-dlon, urcrnrlon-dlon, 100).tolist()
lats += np.linspace(llcrnrlat+dlat, urcrnrlat-dlat, 100).tolist()
# Nothern taper
lons += np.linspace(urcrnrlon-dlon, llcrnrlon+dlon, 100).tolist()
lats += np.linspace(urcrnrlat-dlat, urcrnrlat-dlat, 100).tolist()
# Western taper
lons += np.linspace(llcrnrlon+dlon, llcrnrlon+dlon, 100).tolist()
lats += np.linspace(urcrnrlat-dlat, llcrnrlat+dlat, 100).tolist()

x, y = m(lons,lats)
m.plot(x,y)

# Source
source = namelist.namelist2dict(c_file)['source_mt'][0]
x, y = m(source['lon'], source['lat']) 
# tt, tp, tr, pp, pr, rr
focmecs = source['moment_tensor']
# rr, tt, pp, rt, rp, tp
focmecs = [focmecs[5],focmecs[0],focmecs[3],focmecs[2],focmecs[4],focmecs[1]] 
ax = plt.gca()
b = Beach(focmecs, xy=(x, y), width=100000, linewidth=1, facecolor='r')
b.set_zorder(10)
ax.add_collection(b)


plt.savefig('basemap_plot.png', transparent=True)
plt.show()
