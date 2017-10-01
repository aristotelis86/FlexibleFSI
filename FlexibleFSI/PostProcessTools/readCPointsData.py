#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 01:39:37 2017

@author: aristotelis
"""

import glob
import re
import matplotlib.pyplot as plt

sheetList = glob.glob("../info/sheet*/");

N = len(sheetList);

print(sheetList[0])

pointList = glob.glob(sheetList[0]+"cpoints*.txt");

print(pointList[0])

#for i in range(0,1):
#    s = filelist[i];
#    num = re.split('(\d+)',s);
#    
#    fileID = open(s,"r");
#    lines = fileID.readlines();
#    xpos =[];
#    ypos = [];
#    xvel = [];
#    yvel = [];
#    xforce = [];
#    yforce = [];
#    for x in lines:
#        ll = x.split(",");
#        xpos.append(ll[0]);
#        ypos.append(ll[1]);
#        xvel.append(ll[2]);
#        yvel.append(ll[3]);
#        xforce.append(ll[4]);
#        yforce.append(ll[5]);
#    
#    plt.figure(0)
#    plt.plot(xpos)
#    plt.figure(1)
#    plt.plot(ypos)
#    plt.figure(2)
#    plt.plot(xvel)
#    plt.figure(3)
#    plt.plot(yvel)
#    plt.figure(4)
#    plt.plot(xforce)
#    plt.figure(5)
#    plt.plot(yforce)

    
#f=open("./energy0.txt", "r")


#f1=f.readlines()

#for x in f1:
#    print(x)