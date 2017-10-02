#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt


fSize = 14;

def sheet_plots( sheetN ):
    # reads time, energy and instantaneous length of the sheet specified as input.
    # The input is an integer refering to the ID of the sheet under study.
    # Additional information of the sheet are read as well.
    
    filename = "../info/sheet"+str(sheetN)+"/generalInfo.txt";
    if os.path.isfile(filename):
        infoID = open(filename,"r");
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
    else:
        print("ERROR: \n Folder/File: ../info/sheet"+str(sheetN)+"/generalInfo.txt \n does not exist!")
        return 0,0,0
    
    filename = "../info/sheet"+str(sheetN)+"/energy.txt";
    if os.path.isfile(filename):
        enerID = open(filename,"r");
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
        h.savefig("../info/FIGURES/energy_sheet"+str(sheetN)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,currLength)
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("length of sheet",fontsize = fSize)
        plt.grid()
        plt.title(tit1,fontsize = fSize)
        h.savefig("../info/FIGURES/length_sheet"+str(sheetN)+".png")
        
        return simTime, energy, currLength
    else:
        print("ERROR: Folder/File does not exist!")
        return 0,0,0
        

def point_plots( sheetN, pointN ):
    # The time series of position, velocity and force for a single control
    # point of one sheet are extracted and plotted. A phase space plot is 
    # also generated. integer sheetN is the ID of the sheet and integer 
    # pointN is the ID of the control point.
    
    filename = "../info/sheet"+str(sheetN)+"/generalInfo.txt";
    if os.path.isfile(filename):
        infoID = open(filename,"r");
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
    else:
        print("ERROR: \n Folder/File: ../info/sheet"+str(sheetN)+"/generalInfo.txt \n does not exist!")
        return 0,0,0,0,0,0
    
    filename = "../info/sheet"+str(sheetN)+"/cpoints"+str(pointN)+".txt";
    if os.path.isfile(filename):
        pointID = open(filename,"r");
        lines = pointID.readlines();
        simTime = [];
        xpos =[];
        ypos = [];
        xvel = [];
        yvel = [];
        xforce = [];
        yforce = [];
        for x in lines:
            ll = x.split(",");
            simTime.append(float(ll[0]));
            xpos.append(float(ll[1]));
            ypos.append(float(ll[2]));
            xvel.append(float(ll[3]));
            yvel.append(float(ll[4]));
            xforce.append(float(ll[5]));
            yforce.append(float(ll[6]));
        
        simTime = np.array(simTime);
        xpos = np.array(xpos);
        ypos = np.array(ypos);
        xvel = np.array(xvel);
        yvel = np.array(yvel);
        xforce = np.array(xforce);
        yforce = np.array(yforce);
        
        tit2 = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (Length, Mass, pointN+1, Points, Stiffness, Damping, dt)
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xpos-xpos[0],label="streamwise")
        plt.plot(simTime,ypos-ypos[0],label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("displacement",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/displacement_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xvel,label="streamwise")
        plt.plot(simTime,yvel,label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("velocity",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/velocity_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(simTime,xforce,label="streamwise")
        plt.plot(simTime,yforce,label="cross-stream")
        plt.xlabel("simulation time",fontsize = fSize)
        plt.ylabel("force",fontsize = fSize)
        plt.grid();
        plt.legend();
        plt.title(tit2,fontsize = fSize)
        h.savefig("../info/FIGURES/force_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        
        
        phtit = "Point %d from %d" % (pointN+1, Points)
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
        h.savefig("../info/FIGURES/phase_space_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        return xpos, ypos, xvel, yvel, xforce, yforce
        
    else:
        print("ERROR: Folder/File of the requested control point does not exist!")
        return 0,0,0,0,0,0
    
    
def frequency_analysis( sheetN, pointN ):
    # Perform frequency analysis of the motion of control point (integer) pointN
    # from sheet (integer) sheetN.
    # The relevant figures are plotted.
    
    filename = "../info/sheet"+str(sheetN)+"/generalInfo.txt";
    if os.path.isfile(filename):
        infoID = open(filename,"r");
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
    else:
        print("ERROR: \n Folder/File: ../info/sheet"+str(sheetN)+"/generalInfo.txt \n does not exist!")
        return 0
    
    filename = "../info/sheet"+str(sheetN)+"/cpoints"+str(pointN)+".txt";
    if os.path.isfile(filename):
        pointID = open(filename,"r");
        lines = pointID.readlines();
        simTime = [];
        xpos =[];
        ypos = [];

        for x in lines:
            ll = x.split(",");
            simTime.append(float(ll[0]));
            xpos.append(float(ll[1]));
            ypos.append(float(ll[2]));
            
        simTime = np.array(simTime);
        xpos = np.array(xpos);
        ypos = np.array(ypos);
        
        figtit = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (Length, Mass, pointN+1, Points, Stiffness, Damping, dt)
        
        Nsize = simTime.size;
        
        xfft = fft(xpos);
        yfft = fft(ypos);
        freq = np.linspace(0.0, 1.0/(2.0*dt), Nsize//2)
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(freq, 2.0/Nsize * np.abs(xfft[0:Nsize//2]))
        plt.title(figtit, fontsize=fSize)
        plt.xlabel("frequency",fontsize=fSize)
        plt.ylabel("power",fontsize=fSize)
        plt.xlim(0,50)
        h.savefig("../info/FIGURES/fft_stream_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        
        h = plt.figure(num=None, dpi=100)
        plt.plot(freq, 2.0/Nsize * np.abs(yfft[0:Nsize//2]))
        plt.title(figtit, fontsize=fSize)
        plt.xlabel("frequency",fontsize=fSize)
        plt.ylabel("power",fontsize=fSize)
        plt.xlim(0,50)
        h.savefig("../info/FIGURES/fft_cross_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
        
        return freq, xfft, yfft;
    else:
        print("ERROR: Folder/File of the requested control point does not exist!")
        return 0
        
    