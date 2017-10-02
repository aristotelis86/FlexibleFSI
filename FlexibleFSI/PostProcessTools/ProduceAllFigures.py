#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use this script to produce all the figures available for every sheet and 
# control point of the simulation.

import glob
import readSpecificInfo as myLib


sheetList = glob.glob("../info/sheet*/");
NSheets = len(sheetList);

for i in range(0,NSheets):
    myLib.sheet_plots( i );
    
    pointList = glob.glob(sheetList[i]+"cpoints*.txt");
    NPoints = len(pointList);
    
    for j in range(0,NPoints):
        myLib.point_plots( i, j );
        myLib.frequency_analysis( i, j );
        
        


    
