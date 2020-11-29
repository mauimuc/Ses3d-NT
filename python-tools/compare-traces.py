#!/usr/bin/python
# -*- coding: utf-8 -*-

import obspy.core
import matplotlib.pyplot as plt

st = obspy.core.read( './*.sac')


stations = list( set( [ tr.stats.station for tr in st ] ) )
for station in stations:
    locations = list( set( [ tr.stats.location for tr in st.select(station=station) ] ) )
    for location in locations:
        fig = plt.figure()
        
        st_ = st.select( station=station, location=location )
        channels = list( set( [ tr.stats.channel for tr in st_ ] ) )
        n = len(channels)
        axs = [ fig.add_subplot(n,1,i+1) for i in range(n) ]

        for ax, cn in zip(axs,channels):
            ax.set_ylabel(cn)
            for tr in st_.select(channel=cn):
                ax.plot(tr.data, label=tr.stats.network)

        axs[0].set_title("%s.%s" % (station,location))
        axs[0].legend(frameon=False)

        plt.show()

