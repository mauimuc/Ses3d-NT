#!/usr/bin/python
# -*- coding: utf-8 -*-

import glob
import read_SES3D_trace

for file in glob.glob('S*.[x,y,z]'):
    network, station, location, channel = file.rsplit('.', 4)
    tr = read_seismogram(file, station=station, location=location )
    filename = 'SES3D.%s.%s.%s.sac' % (tr.stats.station, \
                                       tr.stats.location, \
                                       tr.stats.channel )
    tr.write( filename, 'SAC' )

