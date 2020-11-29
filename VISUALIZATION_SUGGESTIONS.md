Visualization Suggestions
=========================

At this section we provide you an example how to visualize Ses3d-NT’s output and
visualization requirements for your Debian based Systems.

Visualization Requirements
==========================

To use the Python toolbox for visualization your system needs to meet the
following requirements:

- Python (≥ 2.5), Python re
- Numpy (≥ 1.1.0)
- Matplotlib (≥ 1.0.0)
- Basemap
- Axisartist
- Obspy (≥ 3.6)
- ParaView
- vtk (≥ 5.0.4)

For Debian based Systems
========================

In the following we guide you through the build process of the dependencies to
meet the requirements to use the Python toolbox. Please follow the instructions
top to bottom.

Python
------
To install Python execute in a shell:

    sudo apt-get install python python2.7-dev

Note: This installs Python 2.7, the module Python re should be contained.

Numpy
-----
To get numpy simply execute:

    sudo apt-get install python-numpy

Matplotlib and Basemap
----------------------

    sudo apt-get install python-matplotlib

    sudo apt-get install python-mpltoolkits.basemap

ObsPy
-----
Informatiom how to install the most recent stable release `ObsPy` under
`Debian` based system you can find on the following link:
https://github.com/obspy/obspy/wiki/Installation-on-Linux-via-Apt-Repository

ParaView
--------
To install ParaView simply type and execute in your shell:

    sudo apt-get install paraview

This also installs vtk.


Visualization of the results
============================

The actual data for visualization can be created by following the instructions
of the L'Aquila toy-example in `GETTING_STARTED.md`. As the result you obtain
SAC and NetCDF file formats data in folder:`ses3d-nt/examples/MESS_2013/`.

Now that you have created the data you also might to take a look on it. There
are a huge variety of programs which can be used for visualizing and processing
Ses3d-NT’s output. In the following only two of them are introduced, namely
ObsPy and ParaView

Plotting Synthetic Seismograms
------------------------------

Ses3d-NT writes synthetic seismograms in the SAC file format, trace by trace.
SAC files can easily be read and processed with ObsPy. To do so, you should
stay in `/ses3d-nt/examples/MESS_2013/` folder (where the L’Aquila simulation
was carried out) and open a Python environment, by executing the following code:

    ipython

        import obspy
        st = obspy.read('./*.sac')
        stations = list (set ( tr.stats.station for tr in st ) )
        for station in stations:
            st.select( station=station ).plot( )

all SAC files in the directory `ses3d-nt/examples/MESS_2013/` will be read
into a stream-object. As a second step all station names are extracted and
stored in a python-list. Finally, the code iterates through all stations,
selecting all those traces corresponding to the current station name and
plotting these. Example of the seismic trace:
@image html Image1.png


Volume-snapshots
----------------

Probably the most user-friendly way to visualize Ses3d-NT’s 3d-datasets is the
3d visualization software ParaView, developed by Kitware. Volume snapshots
which are written in the NetCDF file format can be opened directly, without
conversion. The L’Aquila toy-example produces a whole bunch of output. One file
is `rho_mu_lambda_1.nc` with the simulation’s model-parameters ρ, λ and μ inside.
The other 74 files `vx_vy_vz_[#].nc` are volume-snapshots of displacement-
velocities 'vθ' , 'vφ' and 'vr' . To visualize it:

1. Load ParaView (version 3.98.1; old versions might not be able to load NetCDF),
   click 'open', move to /ses3d-nt/examples/MESS_2013/ and select the group
   `vx_vy_vz_..nc`.
2. In the occurring window, select the NetCDF generic file format.
3. In the left hand panel 'Properties' clicking 'Apply' will show an outline
   of the computational domain.
4. For better data representation, you should switch from 'Outline-view' to the
   'Surface-view' in the top menu.
5. To actually show data included in the NetCDF files, a dataset needs to be
   selected explicitly.

Select from the drop-down menu 'Solid Color' the dataset you are interested in.
ParaView comes with a huge variety of post-processing-, visualization and
filtering-tools. Some features are self-explanatory others are not.

Example of the snapshot:
@image html Image2.png
