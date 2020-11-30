Getting Started
===============

This guide is intended to get you started using Ses3d-NT on commonly used
personal computers. The toy simulation in behind is wrapped around the 2009
L'Aquila Earthquake (for more Details on this example please see the
documentation of the L'Aquila example). We expect you being familiar with a
shell of your choice (e.g. bash), running on a Linux-based computer.


The L'Aquila toy-example
========================

Running the simulation is fairly simple. A configuration file with all the
necessary parameters set is provided in `examples/MESS_2013/config`.
The only ingredients needed, are the number of processes, the simulation shall
be parallelized onto and a MPI execution environment. It depends on your actual
MPI implementation but in most cases the commands are either `mpirun` or
`mpiexec`. Before running toy simulation switch to the `MESS_2013` folder. If
you installed Ses3d-NT in your home-folder, you should get there by typing:

    cd  ~/ses3d-nt/examples/MESS_2013

And then if you build Ses3D-NT **without NetCDF support**, simply run the
simulation by entering:

    mpirun -np 2  ../../bin/ses3d-nt config

where `ses3d-nt` is the actual executable in the `ses3d-nt/bin/` folder and
`config` is the configuration file in the `ses3d-nt/examples/MESS_2013` folder.
The flag `-np` to the MPI execution environment mpirun specifies the number of
processors invoked. We chose `-np 2`, which means two processes. Almost all
recent computers are expected to provide two cores. Feel free to change the
number if your hardware permits.

> Note: The number of processes cannot be chosen arbitrarily. It has to divide
> the number of elements in at least one direction as a whole number.

> Note: In case you built Ses3d-NT linking against your own libraries (e.g. NetCDF)
> keep in mind setting the appropriate library path (e.g.
> `export LD_LIBRARY_PATH=~/local/lib`).

If everything done correct while the simulation runs, you will see the
iterations running through on your standard output. After the end of iterations
you can find SAC and NetCDF file formats data in:
`ses3d-nt/examples/MESS_2013/`.

The example how to visualize these results and visualization requirements for
your Debian based Systems you can find in `VISUALIZATION_SUGGESTIONS.md`

