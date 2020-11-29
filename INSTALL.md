Build Instructions
==================

This file gives an overview on how to build and installing Ses3d-NT. In the 
second half of this document there are actual examples given how to build 
Ses3d-NT on some selected environments. For further information please visit 
the Ses3d-NT web site:

http://www.ses3d-nt.org


Prerequisites
=============

Before you try to build Ses3d-NT make sure your system matches all the 
prerequisites. This is a list of general requirements:

- A current Linux based system
- A commonly used shell (e.g. bash)
- GCC ≥ 4.7.x (preferably gfortran 4.8.2)
- Any MPI ≥ 2.x implementation (e.g. OpenMPI, MPICH, ... )
- GNU Make
- SVN ≥ 1.5

If you want to write output in the NetCDF file format you will additionally need

- zlib (≥ 1.2.7)
- hdf5 parallel (≥ 1.8.9)
- netCDF with MPI support (C and Fortran) (≥ 4.2.x)


Getting the current most source code
====================================

The development takes place under strict revision control. If you meet all
prerequisites you may want to check out the latest version of Ses3d-NT.
The current most source code can be obtained by downloading it from the 
SVN-repository:

    svn checkout https://svn.geophysik.uni-muenchen.de/svn/ses3d/trunk/ ./ses3d-nt

This command will copy all the content of the current development branch named 
'trunk' into the local directory ./ses3d-nt.


Building the Binary
===================

To build the binary, at first switch to the source code directory by typing:

    cd ses3d-nt/src

Building Ses3d-NT is typically a combination of setting some environment
variables and `make`. Environment variables as well as build parameters are 
described below.

The easiest approach is to compile Ses3d-NT **without NetCDF support** by 
executing the following command:

    make ses3d-nt

If everything goes well the executable `ses3d-nt` can be found in the directory 
`ses3d-nt/bin/`. To check whether everything went well you can simply run one 
of the toy simulations. Please see the getting-started guide.

Remarks:
- Without preparations the `make` command will fail on almost all systems. For
  some environments we provide some selected build-examples. See below.
- In case you already tried to run `make` without success or in case you want
  to change something then you have to delete all already compiled object files 
  by executing `make clean` in the first place.
- By default NetCDF support is disabled. To build Ses3d-NT **with parallel
  NetCDF enabled** an environmental variable must be set. Before building
  binaries you should check 
  https://svn.geophysik.uni-muenchen.de/trac/ses3d/wiki/build_instructions/netCDF

Environment Variables
---------------------

As mentioned above some environment variables are recognized. That are:
 - MPI_FC: ...
 - NetCDF: ...
 - MPI_INCLUDE: ...
 - CPPFLAGS: ...
 - FFLAGS: ...
 - LDFLAGS: ...
 - LIBS: ...

There are three ways to set environment variables. Either you put them in front 
of the make command e.g. 

    CPPFLAGS='-DNetCDF' make ses3d-nt

or to 'export' ... 

    export CPPFLAGS='-DNetCDF'

or to alter the values inside of the Makefile. 

Please mind that in case you build with non standard libraries to execute 
Ses3d-NT with specifying the LD_LIBRARY_PATH environment variable.

Build Parameters
----------------

Amongst other parameters the floating-point representation can be specified in
line 45 of `./ses3d-nt/src/parameters.F90`. 

    INTEGER, PARAMETER :: real_kind = REAL64

The default value is REAL32. Possible values are REAL64 and REAL128.

Please mind that the entire code must be recompiled after changing the floating
point representation. Execute `make clean` first!


-------------------------------------------------------------------------------


Examples how to build Ses3d-NT
==============================

Here are some actual examples how to build Ses3d-NT on HPCs (High Performance
Computing) and on Debian based Systems.

SuperMUC
--------

The SuperMUC firewall permits only incoming SSH-connections. You can use port
forwarding to establish a connection between the subversion server and
SuperMUC, i.e., you may use one of the the following procedures.

    ssh -l <LoginName> -R <arbitraryPortNumber>:svn.geophysik.uni-muenchen.de:443 supermuc.lrz.de

    svn checkout https://<remoteLoginName>@localhost:<ForwardedPortNumber>/svn/ses3d/trunk ses3d-nt

To build Ses3d-NT On SuperMUC several modules need to be (un)loaded:

    module unload mpi.ibm
    module load gcc/4.7
    module load mpi.intel/4.1_gcc
    module load netcdf/mpi # optional

and in addition several environment variables must be set:

    export MP_COMPILER=gfortran
    export CPPFLAGS='-D NetCDF' # optional
    export FFLAGS=$NETCDF_INC # optional
    export LDFLAGS=$NETCDF_F90_SHLIB # optional

Ses3d-NT now should compile as usual:

    make -C src ses3d-nt


Segfault & Buserror
-------------------

On segfault and buserror there are current versions GCC, OpenMPI, HDF5 and
netCDF available:

    source /export/data/my_gccvars.sh
    export CPPFLAGS='-D NetCDF'
    export LIBS='-lnetcdff -lnetcdf -lnetcdf'

Ses3d-NT now should compile as usual:

    make -C src ses3d-nt

Note:
In case of previous builds or attempts make sure to enter:

    make -C src clean

before compiling.

Institute (without netCDF)
--------------------------

Considering some restrictions you can build Ses3d-NT on Debian Squeeze.
One Flaw is not having netCDF support.

    source /home/SOFTWARE/GCC/gccvars.sh
    export OMPI_FC=gfortran-4.7.1
    export CPPFLAGS='-D MPI_INCLUDE'


Debian based Systems
====================

In the following we assume you run a apt-based Linux system (e.g. Debian/Mint/
Ubuntu). For Linux Distributions based on other package managers (e.g. RPM), 
usually substituting the term `apt-get` by a term like `YUM` (Fedora) or
`urpmi` (Mageia) should work, please refer to the documentation of your
package manager for details.

During the progress there might be a message that you already have the newest
version of the specific program/library, because you installed it before or it
was shipped with your Linux installation. In that case just go ahead.

First of all update our package list:

    sudo apt-get update

You need to install the following:
- Subversion: The SVN revision control system 
- GNU make: A build automation tool. 
- gfortran: The Fortran compiler of the GNU Compiler Collection.
- OpenMPI: A free MPI implementation. 

To install those packages execute:

    sudo apt-get install make subversion gfortran openmpi-bin libopenmpi-dev
 
Again, make sure having a reasonable new version of GCC. Version 4.7.1 should
be sufficient. It is recommended having 4.8.2 at least. Double check with:

    gfortran --version

Hint:
 On some systems there is a more recent snapshot of GCC available. Try to
 install the package gcc-snapshot. Usually you will find that version in
 `/usr/lib/-gcc-snapshot/bin/`. To invoke this compiler you need to tell it to
 the MPI wrapper script. For OpenMPI with GCC say:

     export OMPI_FC=/usr/lib/gcc-snapshot/bin/gfortran


Optional:
 If you want to use the NetCDF file format to output 3d Volume snapshots, you 
 need to manually install the netCDF C and Fortran libraries and all its 
 dependencies **with parallel support enabled**. Please follow netCDF 
 installation guide 
 https://svn.geophysik.uni-muenchen.de/trac/ses3d/wiki/build_instructions/netCDF
 for further details 

