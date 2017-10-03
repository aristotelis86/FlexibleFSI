#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use this script to produce graphics and extract information for a specific
# sheet or control point.

import readSpecificInfo as myLib

sheetN = 2; # sheet0
pointN = 12; # control point 3 from sheet0 (numbering of cpoints from 0)

myLib.read_sheet_info( sheetN )

#myLib.energy_plot( sheetN )
#myLib.length_plot( sheetN )
#myLib.point_plots( sheetN, pointN )
#myLib.frequency_plot(sheetN, pointN )
myLib.normalMode_frequency_plot( sheetN, pointN, 1, 0.1 )