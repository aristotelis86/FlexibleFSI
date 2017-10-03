#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use this script to produce all the basic figures available for every sheet 
# and control point of the simulation. (Non-specific analysis)
# Sequential numbering for the sheets and control points starting from 0 is 
# assumed.

import glob
import sheetprocessing as myLib
import time

start = time.time()

sheetList = glob.glob("../info/sheet*/");
NSheets = len(sheetList);

for i in range(0,NSheets):
    myLib.energy_plot( i )
    myLib.length_plot( i )
    
    pointList = glob.glob(sheetList[i]+"cpoints*.txt");
    NPoints = len(pointList);
    
    for j in range(0,NPoints):
        myLib.point_plots( i, j )
        myLib.frequency_plot( i, j )
        

print("It took", time.time()-start,"s to run.");

    
