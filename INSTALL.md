Build Instructions
==================

This file gives an overview on how to build and installing Ses3d-NT. In the
second half of this document there are actual examples given how to build
Ses3d-NT on some selected environments.


Prerequisites
=============

Before you try to build Ses3d-NT make sure your system matches all the
prerequisites. This is a list of general requirements:

- A current Linux based system
- A commonly used shell (e.g. bash)
- Fortran 2008 compiler (preferably gfortran ≥ 4.8.2)
- Any MPI ≥ 2.x implementation (e.g. OpenMPI, MPICH, ... )
- GNU Make

If you want to write output in the NetCDF file format you will additionally need

- zlib (≥ 1.2.7)
- hdf5 parallel (≥ 1.8.9)
- netCDF with MPI support (C and Fortran) (≥ 4.2.x)

Having netCDF enabled is handy since visualizing volume snapshots is just few
mouse-clicks away. You may want render your results right away using Paraview.


Building the Binary
===================

To build the executable, switch to the source code directory by typing:

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
  some environments we provide selected build-examples. See below.
- In case you already tried to run `make` without success or in case you want
  to change something then you have to delete all already compiled object files
  by executing `make clean`.
- By default NetCDF support is disabled. To build Ses3d-NT **with parallel
  NetCDF enabled** an environmental variable must be set. See below.

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

use

    export CPPFLAGS='-DNetCDF'

or alter the values inside of the Makefile.

Please mind that in case you build with non standard libraries to execute
Ses3d-NT with specifying the LD_LIBRARY_PATH environment variable.

Build Parameters
----------------

Amongst other parameters the floating-point representation can be specified in
line 45 of /src/parameters.F90

    INTEGER, PARAMETER :: real_kind = REAL64

The default value is REAL32. Possible values are REAL64 and REAL128.

Please mind that the entire code must be recompiled after changing the floating
point representation. Execute `make clean` first!


Debian based Systems
====================

In the following we assume you run an apt-based Linux system (e.g. Debian/Mint/
Ubuntu). For Linux Distributions based on other package managers (e.g. RPM),
usually substituting the term `apt-get` by a term like `YUM` (Fedora) or
`urpmi` (Mageia) should work.

You need to install the following:
- GNU make: A build automation tool.
- gfortran: The Fortran compiler of the GNU Compiler Collection.
- OpenMPI: A free MPI implementation.

To install those packages execute:

    sudo apt-get install make gfortran openmpi-bin libopenmpi-dev

> Hint:
> On some systems there is a more recent snapshot of GCC available. Try to
> install the package gcc-snapshot. Usually you will find that version in
> `/usr/lib/-gcc-snapshot/bin/`. To invoke this compiler you need to tell it to
> the MPI wrapper script. For OpenMPI with GCC say:
> `export OMPI_FC=/usr/lib/gcc-snapshot/bin/gfortran`


> Optional:
> If you want to use the NetCDF file format to output 3d Volume snapshots, you
> need to manually install the netCDF C and Fortran libraries and all its
> dependencies **with parallel support enabled**.
