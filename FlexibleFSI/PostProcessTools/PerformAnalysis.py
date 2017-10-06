#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use this script to produce graphics and extract information for a specific
# sheet and/or control point.

import sheetprocessing as myLib

sheetN = 0; # sheet (numbering of sheets from 0)
pointN = 19; # control point from sheet (numbering of cpoints from 0)

myLib.read_sheet_info( sheetN ) # just print available information on the sheet

myLib.energy_plot( sheetN ) # create+save energy evolution figure

myLib.length_plot( sheetN ) # create+save sheet's length evolution figure

myLib.point_plots( sheetN, pointN ) # create+save plots for position, velocity,
                                    # force and phase space in x,y directions 
                                    # for one control point.

#myLib.frequency_plot( sheetN, pointN ) # create+save plots of the FFT analysis
                                       # of the time series from both x and y 
                                       # position of one control point.

myLib.normalMode_frequency_plot( sheetN, pointN, mode=1, stretchRatio=0.011, align="y" )
# create+save the plot of the FFT analysis of the vibration (pinned-pinned)
# tracked by one control point. The mode (1,2,3,..) is needed as well as the 
# stretching ratio of the entire sheet and its alignment (x,y). 
# This information should be available as an output from READ_SHEET_INFO function.

#myLib.impulse_frequency_analysis( sheetN, pointN, newL=21, g=10, align="y" )
# create+save the plot of the FFT analysis of the vibration (pinned-free)
# tracked by one control point. The first three dominant frequencies are 
# identified and displayed. The stretched length newL is needed, the magnitude 
# of gravity g and the alignment (x,y) of the system.
# This information should be available as an output from READ_SHEET_INFO function.

#########################################################################
# Add automated selection process in normalMode and impulse functions.
# Points should be selected accordingly.
