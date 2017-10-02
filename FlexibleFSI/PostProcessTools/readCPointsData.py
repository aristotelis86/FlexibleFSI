#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import re
import numpy as np
import matplotlib.pyplot as plt

fSize = 14;

sheetList = glob.glob("../info/sheet*/");

NSheets = len(sheetList);

for i in range(0,NSheets):
    infoID = open(sheetList[i]+"generalInfo.txt","r");
    lines = infoID.readlines();
    info = [];

    for x in lines:
        ll = x.split();
        info.append(ll[1]);
    
    Length = float(info[0]);
    Mass = float(info[1]);
    Points = float(info[2]);
    Stiffness = float(info[3]);
    Damping = float(info[4]);
    dt = float(info[5]);
    
#    print("Length="+Length)
#    print("Mass="+Mass)
#    print("Points="+Points)
#    print("Stiffness="+Stiffness)
#    print("Damping="+Damping)
#    print("dt="+dt)
#    print(" ")

    
    enerID = open(sheetList[i]+"energy.txt","r");
    lines = enerID.readlines();
    simTime = [];
    energy = [];
    currLength = [];
    for x in lines:
        ll = x.split(",");
        simTime.append(float(ll[0]));
        energy.append(float(ll[1]));
        currLength.append(float(ll[2]));
    
    simTime = np.array(simTime);
    energy = np.array(energy);
    currLength = np.array(currLength);
    
    tit1 = "Length=%.2f, Mass=%.2f, Points=%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (Length, Mass, Points, Stiffness, Damping, dt)
        
    h = plt.figure(num=None, dpi=100)
    plt.plot(simTime,energy)
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("energy",fontsize = fSize)
    plt.grid()
    plt.title(tit1,fontsize = fSize)
    h.savefig("../info/FIGURES/energy_sheet"+str(i)+".png")
    
    h = plt.figure(num=None, dpi=100)
    plt.plot(simTime,currLength)
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("length of sheet",fontsize = fSize)
    plt.grid()
    plt.title(tit1,fontsize = fSize)
    h.savefig("../info/FIGURES/length_sheet"+str(i)+".png")
    
    pointList = glob.glob(sheetList[i]+"cpoints*.txt");
    NPoints = len(pointList);
    
    for j in range(0,NPoints):
        s = pointList[j];
        num = re.split('(\d+)',s);
        
        pointID = open(s,"r");
        lines = pointID.readlines();
        xpos =[];
        ypos = [];
        xvel = [];
        yvel = [];
        xforce = [];
        yforce = [];
        for x in lines:
            ll = x.split(",");
            xpos.append(float(ll[1]));
            ypos.append(float(ll[2]));
            xvel.append(float(ll[3]));
            yvel.append(float(ll[4]));
            xforce.append(float(ll[5]));
            yforce.append(float(ll[6]));
        
        xpos = np.array(xpos);
        ypos = np.array(ypos);
        xvel = np.array(xvel);
        yvel = np.array(yvel);
        xforce = np.array(xforce);
        yforce = np.array(yforce);
        
        tit2 = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (Length, Mass, j+1, Points, Stiffness, Damping, dt)
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xpos-xpos[0],label="streamwise")
        plt.plot(simTime,ypos-ypos[0],label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("displacement",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/displacement_sheet"+str(i)+"_cp"+str(j)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xvel,label="streamwise")
        plt.plot(simTime,yvel,label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("velocity",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/velocity_sheet"+str(i)+"_cp"+str(j)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xforce,label="streamwise")
        plt.plot(simTime,yforce,label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("force",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/force_sheet"+str(i)+"_cp"+str(j)+".png")
        
        
        phtit = "Point %d from %d" % (j+1, Points)
        h = plt.figure(num=None, dpi=100)
        sb1 = h.add_subplot(121)
        plt.plot(xpos-xpos[0],xvel,label="streamwise")
        plt.xlabel("displacement",fontsize = fSize)
        plt.ylabel("velocity",fontsize = fSize)
        plt.grid()
        plt.title("streamwise",fontsize = fSize)
        sb2 = h.add_subplot(122)
        plt.plot(ypos-ypos[0],yvel,label="cross-stream")
        plt.xlabel("displacement",fontsize = fSize)
        plt.ylabel("velocity",fontsize = fSize)
        sb2.yaxis.set_label_position("right")
        sb2.yaxis.tick_right()
        plt.grid()
        plt.title("cross-stream",fontsize = fSize)
        plt.suptitle(phtit, fontsize = fSize)
        h.savefig("../info/FIGURES/phase_space_sheet"+str(i)+"_cp"+str(j)+".png")
        
    
    
