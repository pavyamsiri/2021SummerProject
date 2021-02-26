# Pavadol Yamsiri's 2021 USyd Summer Project

## Overview

* This repository stores data and scripts I used during my project and hence may not be entirely friendly to navigate/understand.
* The project aimed to compare asteroseismic measurements of angular diameters with interferometrically derived angular diameters to better understand the errors associated with interferometric measurements. The large majority of the project was the analysis of a large data set of stars using TESS data.

### Scripts

* ***plotter.py*** is a script used to plot TESS targets stored in the format outputted by the ***DHTTableProcessor.ipynb***
notebook. It will also output logs for debugging purposes and to track current progress in case of a large target list.

* ***utils.py*** is a set of utility functions used by ***plotter.py***.

* ***star.py*** contains the class **StarData** which is a container for target data, mainly used by ***plotter.py***.

### Notebooks

* ***DHTTableProcessor*** is a notebook that processes the raw target list by Dan Huber (**Angular Diameters December 2019.csv**) to produce dht_plot_data.csv which can be used to plot using ***plotter.py***.

* ***TWTTableProcessor*** is a notebook that processes the raw target list by Tim White (**Interferometric Targets - TimWhite.csv**) to produce twt_plot_data.csv which can be used to plot using ***plotter.py***. They are currently missing angular diameter measurements however.

* ***ResultsAnalyser.ipynb*** is a notebook that processes raw SYD-PYpline output (globalpars.csv) to create comparison plots, comparing seismic and interferometric angular diameters.

* ***DeltanuMeasurer.ipynb*** is a notebook used to edit and measure deltanus for *_comparison_data.csv files.

### Data

* ***RawTargetData*** is the folder containing Dan Huber's target list and Tim White's target list.

* ***Plots*** is the default output parent directory for ***plotter.py*** plots.
    * ***Dan_Huber_Targets*** is the folder of plots of Dan Huber's Targets.
    * ***Tim_White_Targets*** is the folder of plots of Tim White's Targets.
    * ***Comparison*** is the folder of comparison plots.

* ***PlotData*** is the folder used to store plot data that ***plotter.py*** uses to plot targets.
    * ***dht_plot_data.csv*** is the plot data of Dan Huber's targets.
    * ***twt_plot_data.csv*** is the plot data of Dan Huber's targets.
    * The other files are legacy versions of the current plot data format.

* ***Comparison*** is the folder that stores data that is used to make asteroseismology and interferometry measurement comparison plots.

* ***SYD-PYplineResults*** is the folder that the SYD-PYpline's fit background routine's results i.e. globalpars.csv be stored.
