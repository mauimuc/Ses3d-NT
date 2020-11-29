#! /usr/bin/python
# -*- coding: utf-8 -*-

import ses3dpy

# Read configuration file
config = ses3dpy.read_config_file('config')

# Create empty array
field = ses3dpy.raw_empty_field( config )

# rho
field.fill( 1.0/3000.0 )
ses3dpy.write_raw_output( './rhoinv', field, dtype='float32' )

# lambda 
field.fill( 1.1e11 )
ses3dpy.write_raw_output( './lambda', field, dtype='float32' )

# mu
field.fill( 7.0e10 )
ses3dpy.write_raw_output( './mu', field, dtype='float32' )

