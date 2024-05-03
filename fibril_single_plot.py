# fibril_hor_oscillation.py

'''Meant for plotting a sub field of view
with a fibril and its horizontal oscillations
(the cuts along the fibrils are already defined via crispex)
'''

# import modules
import sparsetools as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Cursor
import numpy as np
from sepid import *
import scipy as sci
from Tkinter import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.optimize import curve_fit
from datetime import date
from datetime import datetime
import scipy.ndimage as spnd
import scipy.signal as sgnl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import proj3d
from scipy.signal import argrelextrema as extrem
from scipy.integrate import cumtrapz

# Font properties
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# Directories
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'
invdir = '/home/seki2695/INV/stic/II/slabs_tseries/'
objdir = '/home/seki2695/OUTPUT/project2/class/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcsec
cad = 8
fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

plt.close("all")

#functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

class oscillation:
    def __init__(self, time, vpos_f, vlos_f, ncc_pos, ncc_los, pos_xmax, pos_ymax, pos_xmin, pos_ymin, pos_ncc_xmax, pos_ncc_ymax, los_xmax, los_ymax, los_xmin, los_ymin, los_ncc_xmax, los_ncc_ymax, osc_mid_coord, phase, osc_cat):
        self.time = time
        self.vpos_f = vpos_f
        self.vlos_f = vlos_f
        self.ncc_pos = ncc_pos
        self.ncc_los = ncc_los
        self.pos_xmax = pos_xmax
        self.pos_ymax = pos_ymax
        self.pos_xmin = pos_xmin
        self.pos_ymin = pos_ymin
        self.pos_ncc_xmax = pos_ncc_xmax
        self.pos_ncc_ymax = pos_ncc_ymax
        self.los_xmax = los_xmax
        self.los_ymax = los_ymax
        self.los_xmin = los_xmin
        self.los_ymin = los_ymin
        self.los_ncc_xmax = los_ncc_xmax
        self.los_ncc_ymax = los_ncc_ymax
        self.osc_mid_coord = osc_mid_coord
        self.phase = phase
        self.osc_cat = osc_cat


# boxcar smoothing
"""def smooth(y, w):
    N = y.shape[0]
    r = np.zeros(N)
    for i in range(N):
        if(i==0 or i==N-1):
            r[i] = y[i]
        elif(i>(w-1.)/2. and i<=N-(w+1.)/2.):
            r[i] = np.average(y[int(np.rint(i-w/2.)) : int(np.rint(i+w/2.-1))])
        else:
            r[i] = (y[i-1]+y[i]+y[i+1])/3.
    return r
"""

def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]

        
# Reading the observed data
pref = ['6173','8542','3950']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 =file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak =file_search(datadir,'crispex*'+pref[2]+'*.fcube')

cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])

nx = cube_cak.shape[4]
ny = cube_cak.shape[3]
nt = cube_cak.shape[0]
nw = cube_cak.shape[2]

fibdir = datadir+'big_cut/'

# to test the curve fitting
ti = np.arange(0, nt, 1.)

# plotting elements
tt = 30
rr = 40/1.5

tau_i = 23



# user contribution
for ii in range(len(file_search(fibdir,'crispex*3950*.csav'))):
    print(ii, file_search(fibdir,'crispex*3950*.csav')[ii])

    #ii = int(input("ENTER THE CUT INDEX: "))
    
    cut_file = (file_search(fibdir,'crispex*3950*.csav'))[ii]
    cut = restore(fibdir+cut_file)
    cube = cut.loop_slab[w_pos,:,:]
    
    # cut geometry
    x1 = cut.x_coords[0]
    x2 = cut.x_coords[1]
    y1 = cut.y_coords[0]
    y2 = cut.y_coords[1]
    r = np.sqrt((x1-x2)**2+(y1-y2)**2)
    sinth = np.abs(y2-y1)/r
    costh = np.abs(x2-x1)/r
    
    xx = cut.loop_slab.shape[2]
    
    # Extracting the inversion results
    inv_res = sp.model(invdir+
                       file_search(invdir, "*atmosout*"+cut_file[-11:-5]+"*nc")[0])
    
    temp_cube = inv_res.temp[0]
    vlos_cube = inv_res.vlos[0]
    vturb = inv_res.vturb[0]
    ltau = inv_res.ltau[0,0,0,:]
    ntau = ltau.size
    
    for ll in range(ntau):
        vlos_cube[:,:,ll] = vlos_cube[:,:,ll]
        - spnd.filters.gaussian_filter(inv_res.vlos[0,:,:,ll], [100,100], mode = 'constant')
        
    for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):
        
        #oo = int(input("NO. : "))
        osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
        coord = np.loadtxt(osc_fname)
        amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
        tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
        #arange = np.linspace(amin, amax, amax-amin+1)
        trange = np.linspace(tmin, tmax, tmax-tmin+1)
        
        # coords of the centere of oscillation
        # on the FOV map
        # sin = (y2-y1)/r = (yi-y1)/d
        d = np.mean(coord[:,0])
        yi = d*sinth + y1
        xi = d*sinth + x1
        osc_mid_coord = np.array([xi,yi])
        
        # interpolating the points in between the clicked fit
        y_interp = np.interp(trange, coord[:,1], coord[:,0])
        osc_fname_interp = outdir+"interp"+osc_fname[31:]
        #np.savetxt(osc_fname_interp, y_interp, fmt='%3.8f', encoding=None)
        
        # track max intensity based on the interp click_fit
        y_imax = np.zeros(tmax-tmin+1)
        dist = 1. # width of the interval in which the max intensity is taken
        for i in range(tmax-tmin+1):
            y_imax[i] = np.argmax(cube[i, int(y_interp[i]-dist):int(y_interp[i]+dist)])+int(y_interp[i])-dist
            #print( int(y_interp[i]-dist), int(y_interp[i]+dist))
            
        smth_y = smooth(y_imax, 10)
        osc_fname_smth = outdir+"smth"+osc_fname[31:]
        np.savetxt(osc_fname_smth, smth_y, fmt='%3.8f', encoding=None)
        
        # Calculation of the velocity of the oscillation in the plane of sky
        xp = trange*cad
        fxp = smth_y*res*725
        der = derivative(xp, fxp, dx = 1e-5)
        x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]
        
        # smooth fit to the oscillations
        Imax_smth = np.zeros(tmax-tmin+1)
        temp_smth = np.zeros((tmax-tmin+1,ntau))
        vlos_smth = np.zeros((tmax-tmin+1,ntau))
    
        for i in range(tmax-tmin+1):
            Imax_smth[i] = cube[i,int(np.rint(smth_y[i]))]
            for j in range(ntau):
                temp_smth[i,j] = temp_cube[i,int(np.rint(smth_y[i])),j]
                vlos_smth[i,j] = vlos_cube[i,int(np.rint(smth_y[i])),j]
                
        ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
        ttick_lab = ttick_pos*cad/60
    
        xtick_pos = np.arange(0,(np.round((xx)/rr)+1)*rr,rr)
        xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
        
        # plotting CALCULATIONS
        alph = 0.25
        sm_fact = 8
        sm_fact_los = 8
        per = vlos_smth.shape[0]
        brief_fact = f1x.shape[0]/per
        tt_osc = 25
        ttick_pos_osc = np.arange(0,(np.round((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,tt_osc)
        ttick_lab_osc = (ttick_pos_osc*8).astype('int')
        
        # plot Normalised Cross-Correlation
        ttick_ncc = np.arange(0,(np.round((vlos_smth.shape[0]*2)/tt_osc)+1)*tt_osc,tt_osc) - np.mean(np.arange(0,(np.round((vlos_smth.shape[0]*2)/tt_osc)+1)*tt_osc,tt_osc))
        #ttick_pos_ncc = np.linspace(np.min(ttick_ncc)+5, np.max(ttick_ncc)-5,9)
        #test = ttick_pos_ncc*cad
        #ttick_lab_ncc = test.astype("int")
        tt_osc = 25
        ttick_pos_ac = np.arange(-(np.round((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,(np.round((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,tt_osc)
        ttick_lab_ac = (ttick_pos_ac*8).astype('int')
        ttick_pos_ncc = np.arange(-(np.round((vlos_smth.shape[0])/(tt_osc*2))+1)*tt_osc,(np.round((vlos_smth.shape[0])/(tt_osc*2.))+1)*(tt_osc*2),tt_osc*2)
        ttick_lab_ncc = (ttick_pos_ncc*8).astype('int')
        x_brief = np.zeros(per)
        fx_brief = np.zeros(per)
        f1x_brief = np.zeros(per)
        
        for iii in range(per):
            x_brief[iii] = x[iii*brief_fact]
            fx_brief[iii] = fx[iii*brief_fact]
            f1x_brief[iii] = f1x[iii*brief_fact]
            
        vlos_f = smooth(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5), sm_fact_los) # LOS
        rlos_f = integral(x_brief,vlos_f, dx = 1e-5)
        vpos_f = smooth(f1x_brief - np.mean(f1x_brief), sm_fact) # POS
        rpos_f = smooth(fx_brief - np.mean(fx_brief), sm_fact) # POS
        t = x_brief/cad
        
        vpos_f_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f)*len(vpos_f))
        vpos_ff_norm = (vpos_f-np.mean(vpos_f)) / (np.std(vpos_f))
        vlos_f_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f))
        vlos_ff_norm = (vlos_f-np.mean(vlos_f)) / (np.std(vlos_f)*len(vlos_f))
        ncc = sgnl.correlate(vpos_f_norm, vlos_f_norm, 'full')
        ncc_pos = sgnl.correlate(vpos_f_norm, vpos_ff_norm, 'full') # POS
        ncc_pos_half = ncc_pos[len(ncc_pos)/2:]
        ncc_los = sgnl.correlate(vlos_ff_norm, vlos_f_norm, 'full') # LOS
        ncc_los_half = ncc_los[len(ncc_los)/2:]
        lag = np.arange(-len(vlos_f_norm)+1, len(vpos_f_norm))
        lag_half = lag[len(lag)/2:]
        
        
        # wave properties calculation
        pos_ncc_xmax = lag_half[extrem(ncc_pos_half, np.greater)[0]]
        pos_ncc_ymax = ncc_pos_half[extrem(ncc_pos_half, np.greater)]
        los_ncc_xmax = lag_half[extrem(ncc_los_half, np.greater)[0]]
        los_ncc_ymax = ncc_los_half[extrem(ncc_los_half, np.greater)]
        
        pos_ymax = (vpos_f)[extrem(vpos_f, np.greater, order = 11)[0]]
        pos_xmax = x_brief[extrem(vpos_f, np.greater, order = 11)[0]]/8.
        
        pos_ymin = (vpos_f)[extrem(vpos_f, np.less, order = 11)[0]]
        pos_xmin = x_brief[extrem(vpos_f, np.less, order = 11)[0]]/8.
        
        los_extremax = lag[extrem(ncc_los, np.greater)]
        los_periods_max = los_extremax[np.where(los_extremax>0)]
        los_ymax = (vlos_f)[extrem(vlos_f, np.greater, order = 11)[0]]
        los_xmax = x_brief[extrem(vlos_f, np.greater, order = 11)[0]]/8.
        
        los_extremin = lag[extrem(ncc_los, np.less)]
        los_periods_min = los_extremin[np.where(los_extremin<0)]
        los_ymin = (vlos_f)[extrem(vlos_f, np.less, order = 11)[0]]
        los_xmin = x_brief[extrem(vlos_f, np.less, order = 11)[0]]/8.
        
        # plotting settings
        fig = plt.figure(figsize = (5,4.5))
        gs = gridspec.GridSpec(3,1)
        ax0 = plt.subplot(gs[0,0])
        ax1 = plt.subplot(gs[1,0])
        ax2 = plt.subplot(gs[2,0])
        #ax0.set_title("no. "+ str(oo)+" - "+osc_fname[31:-4])
        ax0.set_xticks(ttick_pos_osc)
        #ax0.set_yticks(ttick_pos_osc)
        ax0.set_xticklabels(ttick_lab_osc, fontdict = font)
        ax0.set_xlabel("t [s]", fontdict = font)
        ax0.tick_params(axis = 'y', labelsize = 8)
        ax0.tick_params(axis = 'x', labelsize = 8)
        ax0.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax0.yaxis.set_minor_locator(AutoMinorLocator(2))
        if (oo==5 or oo==6):
            ax0.set_ylabel(r"$v$ [km s$^{-1}$]", fontdict = font)
        ax0.set_xlim(np.min(x_brief)/8., np.max(x_brief)/8.)
        if (oo==6):
            tit = 'No. 1'
        elif (oo==12):
            tit = 'No. 2'
        elif (oo==5):
            tit = 'No. 3'
        elif (oo==7):
            tit = 'No. 4'
        else:
            tit = 'Fr. ' + str(ii).zfill(2)+str(oo).zfill(2)
        ax0.set_title(tit, fontdict = font)
        
        # plot LOS and POS
        ax0.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        
        ax0.plot(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i])*1e-5, color = 'orangered', linestyle = "--", alpha = alph) # vLOS
        ax0.plot(smooth(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i])*1e-5,sm_fact_los), label = r"LOS", color = 'orangered') # smoothed vLOS
        
        ax0.plot(x_brief/cad, f1x_brief-np.mean(f1x_brief), linestyle = "--", alpha = alph, color = "darkgray") # vPOS
        ax0.plot(x_brief/cad, smooth(f1x_brief-np.mean(f1x_brief),sm_fact), label = r"POS", color = 'darkgray') # smooth vPOS
        
        #ax0.plot(pos_xmax, pos_ymax, "|", color = "black", alpha = 0.9)
        #ax0.plot(los_xmax, los_ymax, "|", color = "orangered", alpha = 0.9)
        #ax0.plot(pos_xmin, pos_ymin, "|", color = "black", alpha = 0.9)
        #ax0.plot(los_xmin, los_ymin, "|", color = "orangered", alpha = 0.9)
        
        # Auto-correlation plots
        ax1.set_xticks(ttick_pos_ac)
        ax1.set_xticklabels(ttick_lab_ac, fontdict = font)
        ax1.set_ylim(-1,1)
        ax1.set_xlim(0, np.max(ttick_ncc))
        ax1.tick_params(axis = 'y', labelsize = 8)
        ax1.set_xlabel(r"time lag [s]", fontdict = font)
        if (oo==5 or oo==6):
            ax1.set_ylabel(r"Auto-$\hat{C}$", fontdict = font)
        ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
        
        ax1.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        ax1.plot(lag, ncc_pos, color = "gray", alpha = 0.75)
        ax1.plot(lag, ncc_los, color = "orangered", alpha = 0.75)
        
        # Cross correlation plot
        ax2.set_xticks(ttick_pos_ncc)
        ax2.set_xticklabels(ttick_lab_ncc, fontdict = font)
        ax2.set_ylim(-1,1)
        ax2.set_xlim(np.min(ttick_ncc), np.max(ttick_ncc))
        ax2.tick_params(axis = 'y', labelsize = 8)
        ax2.set_xlabel(r"time lag [s]", fontdict = font)
        if (oo==5 or oo==6):
            ax2.set_ylabel(r"$\hat{C}$", fontdict = font)
        ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
        
        ax2.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        ax2.axvline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
        ax2.plot(lag, ncc, color = "black", alpha = 0.75)

        if (oo==12):
            fig.legend(loc = "upper right", ncol = 3, fontsize = 8.)
            
        #plt.show()

        plt.tight_layout()
        #fc.savefig(outdir +osc_fname[31:-4]+ "no"+ str(oo)+".pdf", quality = 100)
        osc_fname_plot = outdir+"proto"+osc_fname[31:-4]+"_NCC-"+ "no"+ str(oo)+".pdf"
        #osc_fname_plot = outdir+"cut_profiles/osc"+str(ii).zfill(2)+str(oo).zfill(2)+".png"
        fig.savefig(osc_fname_plot, quality = 100)
        print("file saved to: "+osc_fname_plot)
        plt.close("all")

        #  motion vector 3D
        #------------------------
        if(0):
            y = vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5)
            x = f1x_brief - np.mean(f1x_brief)
            t = x_brief/cad
            
            tt_osc = 12.5
            ttick_pos_osc = np.arange(0,(np.int((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,tt_osc)
            ttick_lab_osc = (ttick_pos_osc*8).astype('int')
            
            inty = cumtrapz(y,t*cad,initial=0)
            intx = cumtrapz(x,t*cad,initial=0)
            intt = t
            
            f3 = plt.figure(figsize = (4,13.5))#max(t)*8./75.))
            
            # Total velocity 3d
            ax = f3.add_subplot(1, 1, 1, projection='3d')#, aspect = ratio)
            ax.azim = 40
            ax.elev = 10
            ax.set_zticks(ttick_pos_osc)
            ax.set_zticklabels(ttick_lab_osc, size = 10.)
            ax.zaxis.set_minor_locator(AutoMinorLocator(8))
            ax.set_zlabel(r"$t$ [s]", size= 12)
            ax.set_ylabel(r"[km]", size = 12)
            ax.set_zlim3d(0,np.max(intt))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(10)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(10)
            ax.set_xlabel(r"[km]", size = 12)
            
            ax.plot3D(np.zeros(len(intx))+np.min(intx), inty, intt, 'orangered', alpha = 0.5, label = r'$A_{\rm{LOS}}$(t)')
            ax.plot3D(intx, np.zeros(len(inty))+np.min(inty), intt, 'rosybrown', alpha = 0.5, label = r'$A_{\rm{POS}}$(t)')
            ax.plot3D(intx, inty, np.zeros(len(intt)), 'lightslategray', linestyle = '--')#, label = r'[$A^2_{\rm{POS}}$(t)+$A^2_{LOS}$(t)]$^{1/2}$'
            ax.plot3D(intx, inty, intt, 'black', label = r'$\mathit{\mathbf{A}}$(t)', alpha = 0.75, linewidth = 2.5)
            ax.set_title(tit)
            '''
            if(oo==6):
                ax.text(np.max(intx), np.min(inty), intt[-1], '(1)', size = 14, color = 'black')
                # make a markerstyle class instance and modify its transform prop
                mark = '<'
                t = mpl.markers.MarkerStyle(marker=mark)
                orig = 12
                t._transform = t.get_transform().rotate_deg(np.arctan((inty[orig+2]-inty[orig])/(intx[orig+2]-intx[orig])))
                ax.scatter(intx[orig], inty[orig], 0, marker = t, s = 50, color = 'lightslategray')
            if(oo==12):
                ax.text(np.max(intx), np.min(inty), intt[-1], '(2)', size = 14, color = 'black')
            if(oo==5):
                ax.text(np.max(intx), np.min(inty), intt[-1], '(3)', size = 14, color = 'black')
            if(oo==7):
                ax.text(np.max(intx), np.min(inty), intt[-1], '(4)', size = 14, color = 'black')
                
            if (oo==6):
            '''
            plt.legend(loc = 'upper right', ncol = 3, fontsize = 12., handlelength = 1.5)

            #plt.tight_layout()
            #plt.show()
            plt.subplots_adjust(left=0.04, bottom = 0, right=0.96, top = 1.)
            #stop()
            
            #f3.savefig(outdir +osc_fname[31:-4]+"-"+ "no"+ str(oo)+"motion.pdf", quality = 100)
            f3.savefig(outdir +"cut_profiles/3d"+str(ii).zfill(2)+str(oo).zfill(2)+".png", quality = 100)
            plt.close(f3)
            

        #  velocity vector 3D
        #------------------------
        if(0):
            #y = smooth(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5), sm_fact)
            #x = smooth(f1x_brief - np.mean(f1x_brief), sm_fact)
            y = vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5)
            x = f1x_brief - np.mean(f1x_brief)
            t = x_brief/cad
            
            tt_osc = 12.5
            ttick_pos_osc = np.arange(0,(np.int((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,tt_osc)
            ttick_lab_osc = (ttick_pos_osc*8).astype('int')
            
            inty = cumtrapz(y,t*cad,initial=0)
            intx = cumtrapz(x,t*cad,initial=0)
            intt = t
            
            f3 = plt.figure(figsize = (4,max(t)*8./75.))
            # Total velocity 3d
            ax = f3.add_subplot(1, 1, 1, projection='3d')#, aspect = ratio)
            ax.azim = 40
            ax.elev = 10
            ax.set_zticks(ttick_pos_osc)
            ax.set_zticklabels(ttick_lab_osc, size = 10.)
            ax.zaxis.set_minor_locator(AutoMinorLocator(8))
            ax.set_zlabel(r"$t$ [s]", size= 12)
            ax.set_ylabel(r"[km s$^{-1}$]", size = 12)
            ax.set_zlim3d(0,np.max(intt))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(10)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(10)
            #ax.set_ylim(-norm_fact, norm_fact)
            ax.set_xlabel(r"[kms$^{-1}$]", size = 12)
            #ax.set_xlim(-norm_fact, norm_fact)
            #ax.plot3D(x, y, t, 'black', label = r'$\mathit{\mathbf{v}}$(t)', alpha = 0.75, linewidth = 2.5)
            #ax.plot3D(np.zeros(len(x))+np.min(x), y, t, 'orangered', alpha = 0.5, label = r'$v_{\rm{LOS}}$(t)')
            #ax.plot3D(x, np.zeros(len(y))+np.min(y), t, 'rosybrown', alpha = 0.5, label = r'$v_{\rm{POS}}$(t)')
            
            ax.plot3D(np.zeros(len(intx))+np.min(intx), inty, intt, 'orangered', alpha = 0.5, label = r'$v_{\rm{LOS}}$(t)')
            ax.plot3D(intx, np.zeros(len(inty))+np.min(inty), intt, 'rosybrown', alpha = 0.5, label = r'$v_{\rm{POS}}$(t)')
            ax.plot3D(intx, inty, np.zeros(len(intt)), 'lightslategray', linestyle = '--')#, label = r'[$A^2_{\rm{POS}}$(t)+$A^2_{LOS}$(t)]$^{1/2}$'
            # make a markerstyle class instance and modify its transform prop
            ax.plot3D(intx, inty, intt, 'black', label = r'$\mathit{\mathbf{A}}$(t)', alpha = 0.75, linewidth = 2.5)
            #if(oo==6):
             #   ax.text(np.max(intx), np.min(inty), intt[-1], '(1)', size = 14, color = 'black')
            #if(oo==12):
             #   ax.text(np.max(intx), np.min(inty), intt[-1], '(2)', size = 14, color = 'black')
            #if(oo==5):
             #   ax.text(np.max(intx), np.min(inty), intt[-1], '(3)', size = 14, color = 'black')
            #if(oo==7):
             #   mark = '<'
              #  t = mpl.markers.MarkerStyle(marker=mark)
               # orig = 12
                #t._transform = t.get_transform().rotate_deg(np.arctan((inty[orig+2]-inty[orig])/(intx[orig+2]-intx[orig])))
                #ax.scatter(intx[orig], inty[orig], 0, marker = t, s = 50, color = 'lightslategray')
                #ax.scatter(intx, inty)
                #ax.text(np.max(intx), np.min(inty), intt[-1], '(4)', size = 14, color = 'black')
                
            #ax.plot3D(x, y, np.zeros(len(t)), 'darkgray', alpha = 0.5)
            #ax.set_aspect('equal', 'box')
            if (oo==6):
                plt.legend(loc = 'upper right', ncol = 3, fontsize = 12., handlelength = 1.5)
            plt.tight_layout()
            plt.show()
            plt.subplots_adjust(left=0.04, bottom = 0, right=0.96, top = 1.)
            stop()
            
            f3.savefig(outdir +osc_fname[31:-4]+"-"+ "no"+ str(oo)+".pdf", quality = 100)
            plt.close(f3)
            
            
            
            
    #stop()
    plt.close("all")
        
        
        
