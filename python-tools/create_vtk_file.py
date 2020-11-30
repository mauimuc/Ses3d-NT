#!/usr/bin/python
# -*- coding: utf-8 -*-

# create_vtk_file.py

# by Stefan Mauerberger
#    Michael Haas

# last modified: $Date: 2012-10-26 18:31:45 +0200 (Fri, 26 Oct 2012) $

#######################################################################
# this python script converts SES3Ds binary files into the VTK-format #
#######################################################################

from numpy import pi, sin, cos
import ses3dpy
from argparse import ArgumentParser

#### argument passed to cmd-line declaring which time-step to use ####
parser = ArgumentParser("Create vtk files for parsed volume-snapshot files")
parser.add_argument("-i",'--input', help="Rank 6 volume-snapshot files",\
                    nargs='+', type=str, required=True )
parser.add_argument("-o", '--output', help="Output file-name",\
                    type=str, required=True )
parser.add_argument("-c",'--config', help="Ses3d-NT config. file name",\
                    type=str, required=True )
args=parser.parse_args()

output_file = args.output
config_file = args.config
file_names = args.input

# read parameter file
par = ses3dpy.read_config_file( config_file ) # read parameters of the simulation

# generate coordinates
t,p,r = ses3dpy.generate_coordinates( par, no_dup=False )
t     = t * pi / 180
p     = p * pi / 180
# 3d coordinate mesh
T,P,R = ses3dpy.mesh3d(t,p,r)
# transform to Cartesian coordinates
X = R * sin(T) * cos(P)
Y = R * sin(T) * sin(P)
Z = R * cos(T)
X = ses3dpy.pop_duplicates( par, X )
Y = ses3dpy.pop_duplicates( par, Y )
Z = ses3dpy.pop_duplicates( par, Z )

fields = {} # init dictionary
for file_name in file_names:
    field = ses3dpy.read_field( par, file_name )
    field = ses3dpy.pop_duplicates( par, field )
    fields[file_name] = ( field, )
# write vtk file
ses3dpy.write_vtk( (X,Y,Z), fields, output_file )
# output location
print 'Created %s' % output_file

