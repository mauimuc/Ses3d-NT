#!/usr/bin/python
# -*- coding: utf-8 -*-

def read_seismogram( filename, station='', location='' ):
    """Reads Andreas' legacy SES3D ASCII synthetic seismograms. 
An ObsPy-trace object is returned.  

Parameters:
-----------
filename : string 
station: string 
location: string """

    from re         import compile
    from obspy.core import Trace
    from numpy      import loadtxt

    with open( filename, 'r' ) as f:    # open the file
        channel = f.readline()  # channel 
        if channel.startswith(' theta'):
            channel = 'vx'
        elif channel.startswith( ' phi' ):
            channel = 'vy'
        elif channel.startswith( ' r' ):
            channel = 'vz'
        else: 
            raise Exception('Channel not known!')
        f.readline()            # read 2ed line
        dt = float( f.readline().split('=')[1] ) # sampling-rate
        f.readline()            # 4. ...
        pattern = compile( r'([a-z]+)=\s*(\-?\d+\.\d+)' ) # reg-ex
        coord = {}              # receiver coordinates
        for match in pattern.finditer( f.readline() ):
            if ( match.group(1) == 'x' ):
                coord['lat'] = 90.0 - float( match.group(2) )
            elif ( match.group(1) == 'y' ):
                coord['lon'] = float( match.group(2) )
            elif ( match.group(1) == 'z' ):
                coord['depth'] = float( match.group(2) )
            else:
                coord[match.group(1)] = float( match.group(2) )
        f.readline() 			# 6. ...
        source = {}             # source coordinates
        for match in pattern.finditer( f.readline() ):
            if ( match.group(1) == 'x' ):
               source['lat'] = 90.0 - float( match.group(2) )
            elif ( match.group(1) == 'y' ):
               source['lon'] = float( match.group(2) )
            elif ( match.group(1) == 'z' ):
               source['depth'] = float( match.group(2) )
            else:
               source[match.group(1)] = float( match.group(2) )
        # reads data from file 
        data = loadtxt( f )
        # reads the data and stores it in a ObsPy-Trace within parameters
        return Trace( data, { 'delta': dt, \
                              'network': 'SES3D', \
                              'station': station, \
                              'location': location, \
                              'channel': channel, \
                              'coord': coord, \
                              'source': source } ) 

