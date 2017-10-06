#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt


plt.ioff();
fSize = 14;
myDpi = 100;
figSize=(6,5);

def read_sheet_info( sheetN, show=True ):
    # Reads information file of (integer) sheetN.
    # Relevant parameters are printed to cmd line.
    
    filename = "../info/sheet"+str(sheetN)+"/generalInfo.txt";
    if os.path.isfile(filename):
        infoID = open(filename,"r");
        lines = infoID.readlines();
        info = [];

        for x in lines:
            ll = x.split();
            info.append(ll[-1]);
    
        Length = float(info[0]);
        Mass = float(info[1]);
        Points = float(info[2]);
        Stiffness = float(info[3]);
        Damping = float(info[4]);
        dt = float(info[5]);
        
        if show:
            print("Quantities in normalized units.\n")
            print("Length=%.2f" % Length)
            print("Mass=%.2f" % Mass)
            print("Number of points=%d" % Points)
            print("Stiffness=%.2f" % Stiffness)
            print("Damping=%.2f" % Damping)
            print("Time step size=%f" % dt)
            if (len(lines)>6):
                for ij in range(6,len(lines)):
                    print(lines[ij]);
        
        return Length, Mass, Points, Stiffness, Damping, dt
    else:
        print("ERROR: \n Folder/File: ../info/sheet"+str(sheetN)+"/generalInfo.txt \n does not exist!")
        return 0,0,0,0,0,0
    
    
def read_energy_file( sheetN ):
    # Reads energy file of (integer) sheetN.
    # The output of the function is the time, energy and instantaneous length.
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

        return simTime, energy, currLength
    else:
        print("ERROR: Folder/File does not exist!")
        return 0,0,0
    
def read_cpoint_file( sheetN, pointN ):
    # Reads the file of control point (integer) pointN from sheet
    # (integer) sheetN. 
    # The output of the function is time, x,y-position, x,y-velocity,
    # x,y-force. 
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
        
        return simTime, xpos, ypos, xvel, yvel, xforce, yforce
    else:
        print("ERROR: Folder/File of the requested control point does not exist!")
        return 0,0,0,0,0,0,0
    

def energy_plot( sheetN ):
    # Creates the figure for the total energy evolution.
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, energy, cLength = read_energy_file( sheetN );
    
    figTitle = "Length=%.2f, Mass=%.2f, Points=%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, P, S, D, dt)
        
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(simTime,energy)
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("energy",fontsize = fSize)
    plt.grid()
    plt.title(figTitle,fontsize = fSize)
    h.savefig("../info/FIGURES/energy_sheet"+str(sheetN)+".png")
    plt.close(h)
    
    return;
    
def length_plot( sheetN ):
    # Creates the figure for the evolution of the length of the sheet.
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, energy, cLength = read_energy_file( sheetN );
    
    figTitle = "Length=%.2f, Mass=%.2f, Points=%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, P, S, D, dt)
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(simTime,cLength)
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("length of sheet",fontsize = fSize)
    plt.grid()
    plt.title(figTitle,fontsize = fSize)
    h.savefig("../info/FIGURES/length_sheet"+str(sheetN)+".png")
    plt.close(h)
    
    return;
        

def point_plots( sheetN, pointN ):
    # The time series of position, velocity and force for a single control
    # point of one sheet are plotted. A phase space plot is also generated.
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, xpos, ypos, xvel, yvel, xforce, yforce = read_cpoint_file( sheetN, pointN );
    
    figTitle = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, pointN+1, P, S, D, dt)
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(simTime,xpos-xpos[0],label="streamwise")
    plt.plot(simTime,ypos-ypos[0],label="cross-stream")
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("displacement",fontsize = fSize)
    plt.grid();
    plt.legend();
    plt.title(figTitle,fontsize = fSize)
    h.savefig("../info/FIGURES/displacement_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
        
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(simTime,xvel,label="streamwise")
    plt.plot(simTime,yvel,label="cross-stream")
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("velocity",fontsize = fSize)
    plt.grid();
    plt.legend();
    plt.title(figTitle,fontsize = fSize)
    h.savefig("../info/FIGURES/velocity_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
        
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(simTime,xforce,label="streamwise")
    plt.plot(simTime,yforce,label="cross-stream")
    plt.xlabel("simulation time",fontsize = fSize)
    plt.ylabel("force",fontsize = fSize)
    plt.grid();
    plt.legend();
    plt.title(figTitle,fontsize = fSize)
    h.savefig("../info/FIGURES/force_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
        
    phtit = "Point %d from %d" % (pointN+1, P)
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
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
    plt.close(h)

    return;
    
def frequency_plot( sheetN, pointN ):
    # General frequency plot of the motion of control point 
    # (integer) pointN from sheet (integer) sheetN. 
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, xpos, ypos, xvel, yvel, xforce, yforce = read_cpoint_file( sheetN, pointN );
    
    figtit = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, pointN+1, P, S, D, dt)
        
    freq, spec = custom_FFT( xpos, dt );
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(freq, spec)
    plt.title(figtit, fontsize=fSize)
    plt.xlabel("frequency",fontsize=fSize)
    plt.ylabel("power",fontsize=fSize)
    plt.xlim(0,5)
    plt.legend()
    plt.grid();
    h.savefig("../info/FIGURES/fft_stream_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
    
    freq, spec = custom_FFT( ypos, dt );
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(freq, spec)
    plt.title(figtit, fontsize=fSize)
    plt.xlabel("frequency",fontsize=fSize)
    plt.ylabel("power",fontsize=fSize)
    plt.xlim(0,5)
    plt.legend()
    plt.grid();
    h.savefig("../info/FIGURES/fft_cross_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
    
    return;


def normalMode_frequency_plot( sheetN, pointN, mode, stretchRatio, align="y" ):
    # Frequency plot of the motion of control point 
    # (integer) pointN from sheet (integer) sheetN.
    # Information related to the normal mode under 
    # study is also printed on the figure and cmd.
    # The number of the mode should be input and 
    # the ratio of the initial stretch of the sheet.
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, xpos, ypos, xvel, yvel, xforce, yforce = read_cpoint_file( sheetN, pointN );
    
    figtit = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, pointN+1, P, S, D, dt)

    expFreq = np.sqrt(P*S*stretchRatio/M)*mode/(2*L)
#    expFreq = np.sqrt(P*S/M)*mode/(2*L)
    
    if (align=="x"):
        freq, spec = custom_FFT( ypos, dt );
    else:
        freq, spec = custom_FFT( xpos, dt );
    PeakInd = detect_peaks(spec);
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(freq, spec, label="data")
    plt.plot([expFreq, expFreq],[0, max(spec)], label="expected")
    plt.title(figtit, fontsize=fSize)
    plt.xlabel("frequency",fontsize=fSize)
    plt.ylabel("power",fontsize=fSize)
    plt.xlim(0,3*freq[PeakInd[0]])
    plt.legend()
    plt.grid();
    h.savefig("../info/FIGURES/fft_mode"+str(mode)+"_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
    
    print("Expected resonance frequency was %.3f and was found %.3f" % (expFreq, freq[PeakInd[0]]))
    return;
        
def impulse_frequency_analysis( sheetN, pointN, newL, g=10, align="y" ):
    # Perform frequency analysis of the motion of control point 
    # (integer) pointN from sheet (integer) sheetN. Refers to
    # the impulse test. The magnitude of gravity used must be 
    # given. 
    
    L, M, P, S, D, dt = read_sheet_info( sheetN, show=False );
    simTime, xpos, ypos, xvel, yvel, xforce, yforce = read_cpoint_file( sheetN, pointN );
    
    figtit = "Length=%.2f, Mass=%.2f, # Point=%d/%d \n Stiffness=%.1f, Damping=%.1f, dt=%.2e" % (L, M, pointN+1, P, S, D, dt)
    
    if (align=="x"):
        freq, spec = custom_FFT( ypos, dt );
    else:
        freq, spec = custom_FFT( xpos, dt );
    PeakInd = detect_peaks(spec)
    
    expFreq1 = (2.4048/(4*np.pi))*np.sqrt(g/newL);
    expFreq2 = (5.5021/(4*np.pi))*np.sqrt(g/newL);
    expFreq3 = (8.6537/(4*np.pi))*np.sqrt(g/newL);
    
    h = plt.figure(num=None, figsize=figSize, dpi=myDpi)
    plt.plot(freq, spec, label="data")
    plt.plot([expFreq1, expFreq1],[0, max(spec)], label="expected 1")
    plt.plot([expFreq2, expFreq2],[0, max(spec)], label="expected 2")
    plt.plot([expFreq3, expFreq3],[0, max(spec)], label="expected 3")
    plt.title(figtit, fontsize=fSize)
    plt.xlabel("frequency",fontsize=fSize)
    plt.ylabel("power",fontsize=fSize)
    plt.xlim(0,5)
    plt.legend()
    plt.grid();
    h.savefig("../info/FIGURES/fft_impulse_sheet"+str(sheetN)+"_cp"+str(pointN)+".png")
    plt.close(h)
    
    return;
    
def custom_FFT( data, stepT ):
    # Customize fft analysis.
    Nsize = data.size;
    fftOut = fft(data);
    freq = np.linspace(0.0, 1.0/(2.0*stepT), Nsize//2)
    
    return freq[1:], 2.0/Nsize * np.abs(fftOut[1:Nsize//2])
#    return freq, 2.0/Nsize * np.abs(fftOut[0:Nsize//2])


# Found it on: http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indices of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind
    
def _plot(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02*x.size, x.size*1.02-1)
        ymin, ymax = x[np.isfinite(x)].min(), x[np.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('Data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()   
    