#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use this script to produce graphics and extract information for a specific
# sheet or control point.

import readSpecificInfo as myLib

sheetN = 0; # sheet0
pointN = 3; # control point 3 from sheet0 (numbering of cpoints from 0)

simTime, energy, currLength = myLib.sheet_plots( sheetN )
xpos, ypos, xvel, yvel, xforce, yforce = myLib.point_plots( sheetN, pointN )
freq, amp  = myLib.frequency_analysis( sheetN, pointN )

