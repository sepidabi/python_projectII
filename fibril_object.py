# import modules
import sparsetools as sp
import matplotlib.pyplot as plt
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
from scipy.signal import argrelextrema as extrem

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
#fr_select = [29,75,135,167] # frames to display
#w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

plt.close("all")

#functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

class oscillation:
    def __init__(self, time, vpos_f, vlos_f, pos_xmax, pos_ymax, pos_xmin, pos_ymin, pos_ncc_xmax, pos_ncc_ymax, los_xmax, los_ymax, los_xmin, los_ymin, los_ncc_xmax, los_ncc_ymax):
        self.time = time
        self.vpos_f = vpos_f
        self.vlos_f = vlos_f
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


# boxcar smoothing

def smooth(x,window_len=11,window='flat'):

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

# plotting settings
fig = plt.figure(figsize = (8,3))
gs = gridspec.GridSpec(1,1)
ax0 = plt.subplot(gs[0,0])

# GUI stuff
cursor = Cursor(ax0, horizOn=True, vertOn=True, useblit=True,
                color = 'black', linewidth = 0.75, alpha = 0.75)
def onclick(event):
    xx, tt = event.xdata, event.ydata
    global coord
    coord.append((xx, tt))
    ax0.plot(xx,tt, 'x', color = 'dodgerblue')
    print('t [s]=%1.2f, v_los [km/s]=%1.2f' %
          (np.round(xx, decimals = 2),
           np.round(tt, decimals = 2)))
    fig.canvas.draw() #redraw the figure

# storing oscillation properties in object
if (1):

    # User contribution to the specify the oscillation
    for ii in range(len(file_search(objdir, "cut*osc0.obj"))):
        print(ii, file_search(objdir, "cut*osc0.obj")[ii][3:12])
    cc = input("Enter the number of the cut: ")
    
    for iii in range(cc,len(file_search(objdir, "*osc0.obj"))):
        cut_indx = file_search(objdir, "*osc0.obj")[iii][3:12]
        #oo = input('osc no. (0-'+str(len(file_search(objdir, 'cut'+cut_indx+'*osc*.obj'))-1)+'): ')
        for jjj in range(len(file_search(objdir, 'cut'+cut_indx+'*osc*.obj'))):
            osc_indx = str(jjj)
            
            osc_obj = (file_search(objdir, "cut*"+cut_indx+"_osc"+osc_indx+".obj"))[0]
            
            osc = load_obj(objdir+osc_obj)
            
            vlos_f = osc.vlos_f
            vpos_f = osc.vpos_f
            x = osc.time
            pos_xmax = osc.pos_xmax
            los_xmax = osc.los_xmax
            pos_ymax = osc.pos_ymax
            los_ymax = osc.los_ymax
            pos_xmin = osc.pos_xmin
            los_xmin = osc.los_xmin
            pos_ymin = osc.pos_ymin
            los_ymin = osc.los_ymin
            
            # Plotting
            ax0.cla()
            ax0.tick_params(axis = 'y', labelsize = 8)
            ax0.tick_params(axis = 'x', labelsize = 8)
            ax0.xaxis.set_minor_locator(AutoMinorLocator(8))
            ax0.set_ylabel(r"$v$ [km s$^{-1}$]", fontdict = font)
            ax0.set_xlim(0,np.max(x))
            
            ax0.axhline(0, linestyle = ":", color = "gray", linewidth = 1, alpha = 0.5)
            
            coord = []
            
            for i in range(4):
                cid = fig.canvas.mpl_connect('button_press_event', onclick)
                plt.show()
                #cir = fig.canvas.mpl_connect('<Button-3>', other_function)
                
                if i==0:
                    ax0.plot(x, vlos_f, color = 'orangered') # smoothed vLOS
                    ax0.plot(x,vpos_f, color = 'darkgray', alpha = 0.5) # smooth vPOS
                    print("Mark the MAX on v_los!")
                    ax0.plot(los_xmax, los_ymax, "|", color = "orangered", alpha = 0.9)
                    plt.legend()
                    stop()
                    if(coord!=[]):
                        osc.los_xmax = np.array(coord)[:,0]
                        osc.los_ymax = np.array(coord)[:,1]
                        coord = []

                if i==1:
                    ax0.cla()
                    ax0.tick_params(axis = 'y', labelsize = 8)
                    ax0.tick_params(axis = 'x', labelsize = 8)
                    ax0.xaxis.set_minor_locator(AutoMinorLocator(8))
                    ax0.set_ylabel(r"$v$ [km s$^{-1}$]", fontdict = font)
                    ax0.set_xlim(0,np.max(x))
                    ax0.plot(x, vlos_f, label = r"LOS", color = 'orangered') # smoothed vLOS
                    ax0.plot(x,vpos_f, color = 'darkgray', alpha = 0.5) # smooth vPOS
                    print("Mark the MIN on v_los!")
                    ax0.plot(los_xmin, los_ymin, "|", color = "orangered", alpha = 0.9)
                    plt.legend()
                    stop()
                    if(coord!=[]):
                        osc.los_xmin = np.array(coord)[:,0]
                        osc.los_ymin = np.array(coord)[:,1]
                        coord = []
                        
                if i==2:
                    ax0.cla()
                    ax0.tick_params(axis = 'y', labelsize = 8)
                    ax0.tick_params(axis = 'x', labelsize = 8)
                    ax0.xaxis.set_minor_locator(AutoMinorLocator(8))
                    ax0.set_ylabel(r"$v$ [km s$^{-1}$]", fontdict = font)
                    ax0.set_xlim(0,np.max(x))
                    ax0.plot(x, vlos_f, label = r"LOS", color = 'orangered', alpha = 0.5) # smoothed vLOS
                    ax0.plot(x,vpos_f, color = 'darkgray') # smooth vPOS
                    print("Mark the MAX on v_pos!")
                    ax0.plot(pos_xmax, pos_ymax, "|", color = "black", alpha = 0.9)
                    plt.legend()
                    stop()
                    if(coord!=[]):
                        osc.pos_xmax = np.array(coord)[:,0]
                        osc.pos_ymax = np.array(coord)[:,1]
                        coord = []
                        
                if i==3:
                    ax0.cla()
                    ax0.tick_params(axis = 'y', labelsize = 8)
                    ax0.tick_params(axis = 'x', labelsize = 8)
                    ax0.xaxis.set_minor_locator(AutoMinorLocator(8))
                    ax0.set_ylabel(r"$v$ [km s$^{-1}$]", fontdict = font)
                    ax0.set_xlim(0,np.max(x))
                    ax0.plot(x, vlos_f, label = r"LOS", color = 'orangered', alpha = 0.5) # smoothed vLOS
                    ax0.plot(x,vpos_f, label = r"POS", color = 'darkgray') # smooth vPOS
                    print("Mark the MIN on v_pos!")
                    ax0.plot(pos_xmin, pos_ymin, "|", color = "black", alpha = 0.9)
                    plt.legend()
                    stop()
                    if(coord!=[]):
                        osc.pos_xmin = np.array(coord)[:,0]
                        osc.pos_ymin = np.array(coord)[:,1]
                        coord = []

            print("Press r to save the object and \n Continue to next oscillation...")
            stop()
            
            save_obj(objdir+osc_obj, osc)

per_los = np.zeros(0)
per_los_m = np.zeros(0)
per_los_dom = np.zeros(0)
per_pos = np.zeros(0)
per_pos_m = np.zeros(0)
per_pos_dom = np.zeros(0)
amp_los = np.zeros(0)
amp_los_m = np.zeros(0)
amp_pos = np.zeros(0)
amp_pos_m = np.zeros(0)
A_los_one = np.zeros(0)
A_los = np.zeros(0)
A_los_m = np.zeros(0)
A_pos_one = np.zeros(0)
A_pos = np.zeros(0)
A_pos_m = np.zeros(0)
dphi = np.zeros(0)

plt.close("all")
# testing the LOS values
test_fig = plt.figure()
ax = test_fig.add_subplot(1, 1, 1)

A_los_test = np.zeros(0)
A_pos_test = np.zeros(0)
 
# Extract the oscillation objects
if (1):

    # User contribution to the specify the oscillation
    for iii in range(len(file_search(objdir, "*osc0.obj"))):
        cut_indx = file_search(objdir, "*osc0.obj")[iii][3:12]

        for jjj in range(len(file_search(objdir, 'cut'+cut_indx+'*osc*.obj'))):
            osc_indx = str(jjj)
            print(osc_indx)
            osc_obj = (file_search(objdir, "cut*"+cut_indx+"_osc"+osc_indx+".obj"))[0]
            osc = load_obj(objdir+osc_obj)

            # test
            if ((len(osc.los_ymin)==0 or len(osc.los_ymax)==0)
                or (len(osc.pos_ymin)==0 or len(osc.pos_ymax)==0)):
                print("doing nothing for now")
            else:
                A_los_test = np.append(A_los_test, [np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))])
                #ax.plot(osc.time, osc.vlos_f, color = "black", alpha = 0.01)
                ax.scatter(x = osc.los_xmax, y = np.abs(osc.los_ymax), s = 5, color = 'orangered', alpha = 0.25)
                ax.scatter(x = osc.los_xmin, y = np.abs(osc.los_ymin), s = 5, color = 'orangered', alpha = 0.25)
                ax.scatter(x = [np.mean(osc.los_xmax), np.mean(osc.los_xmin)],
                           y = [np.mean(np.abs(osc.los_ymax)), np.mean(np.abs(osc.los_ymin))],
                           color = 'red', alpha = 0.5, s = 10)
                A_pos_test = np.append(A_pos_test, [np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))])
                ax.scatter(x = osc.pos_xmax, y = np.abs(osc.pos_ymax), s = 5, color = 'cyan', alpha = 0.25)
                ax.scatter(x = osc.pos_xmin, y = np.abs(osc.pos_ymin), s = 5, color = 'cyan', alpha = 0.25)
                ax.scatter(x = [np.mean(osc.pos_xmax), np.mean(osc.pos_xmin)],
                           y = [np.mean(np.abs(osc.pos_ymax)), np.mean(np.abs(osc.pos_ymin))],
                           color = 'blue', alpha = 0.5, s = 10)
                

            # OSCILLATION VALIDITY CHECK
            if (len(osc.los_ncc_xmax)==0 or len(osc.pos_ncc_xmax)==0):
                print("doing nothing for now")
            else:

                #--------------------------------
                # LINE of SIGHT oscillation
                #--------------------------------
                # LOS Amplitude
                a = np.abs(osc.los_ymin)
                b = np.abs(osc.los_ymax)
                if (len(a)>len(b)):
                    A_l = np.array([a[0]])
                    for cc in range(len(b)):
                        A_l = np.append(A_l, (b[cc],a[cc+1]))
                if (len(b)>len(a)):
                    A_l = np.array([b[0]])
                    for cc in range(len(a)):
                        A_l = np.append(A_l, (a[cc], b[cc+1]))
                if (len(a)==len(b)):
                    A_l = np.zeros(0)
                    for cc in range(len(a)):
                        A_l = np.append(A_l, (a[cc],b[cc]))
                        
                if (len(A_l)==1):
                    A_los_one = A_l
                if (len(A_l)==2):
                    A_los_one = np.mean(A_l)
                if (len(A_l)>=3):
                    for aa in range(len(A_l)-2):
                        A_los_one = np.append(A_los_one, A_l[aa]+A_l[aa+1]+A_l[aa+2])
                A_los = np.append(A_los, A_los_one)
                A_los_m = np.append(A_los_m, np.mean(A_los_one))

                # LOS Period
                pe_l = np.sort(np.append(osc.los_xmax, osc.los_xmin))
                per_los_one = np.zeros(0)

                if (len(pe_l)==2):
                    per_los_one = np.append(per_los_one, 2*np.abs(osc.los_xmax-osc.los_xmin))
                if (len(pe_l)==3):
                    per_los_one = np.append(per_los_one, pe_l[-1] - pe_l[0])
                if (len(pe_l)>3):
                    for l in range (len(pe_l)-3):
                        per_los_one = np.append(per_los_one, (pe_l[l+3]-pe_l[l+1]+pe_l[l+2]-pe_l[l])/2.)
                if (len(pe_l)==1):
                    per_los_one = pe_l*2
                        
                per_los = np.append(per_los, per_los_one)
                per_los_dom = np.append(per_los_dom, per_los_one[np.argmin(np.abs(per_los_one-osc.los_ncc_xmax[0]))])
                per_los_m = np.append(per_los_m, np.mean(per_los_one))


                #-------------------------------
                # PLANE of SKY oscillation
                #--------------------------------
                # POS Amplitude
                a = np.abs(osc.pos_ymin)
                b = np.abs(osc.pos_ymax)
                if (len(a)>len(b)):
                    A_p = np.array([a[0]])
                    for cc in range(len(b)):
                        A_p = np.append(A_p, (b[cc],a[cc+1]))
                if (len(b)>len(a)):
                    A_p = np.array([b[0]])
                    for cc in range(len(a)):
                        A_p = np.append(A_p, (a[cc], b[cc+1]))
                if (len(a)==len(b)):
                    A_p = np.zeros(0)
                    for cc in range(len(a)):
                        A_p = np.append(A_p, (a[cc],b[cc]))
                if (len(A_p)==1):
                    A_pos_one = A_p
                if (len(A_p)==2):
                    A_pos_one = np.mean(A_p)
                if (len(A_p)>=3):
                    for aa in range(len(A_p)-2):
                        A_pos_one = np.append(A_pos_one, A_p[aa]+A_p[aa+1]+A_p[aa+2])
                A_pos = np.append(A_pos, A_pos_one)
                A_pos_m = np.append(A_pos_m, np.mean(A_pos_one))

                # POS period
                pe_p = np.sort(np.append(osc.pos_xmax, osc.pos_xmin))
                per_pos_one = np.zeros(0)

                if (len(pe_p)==2):
                    per_pos_one = np.append(per_pos_one, 2*np.abs(osc.pos_xmax-osc.pos_xmin))
                if (len(pe_p)==3):
                    per_pos_one = np.append(per_pos_one, pe_p[-1] - pe_p[0])
                if (len(pe_p)>3):
                    for m in range (len(pe_p)-3):
                        per_pos_one = np.append(per_pos_one, (pe_p[m+3]-pe_p[m+1]+pe_p[m+2]-pe_p[m])/2.)
                    
                per_pos = np.append(per_pos, per_pos_one) 
                per_pos_dom = np.append(per_pos_dom, per_pos_one[np.argmin(np.abs(per_pos_one-osc.pos_ncc_xmax[0]))])
                per_pos_m = np.append(per_pos_m, np.mean(per_pos_one))

    ##########
    # JOINT PLOTS
    ##########
    fontsize = 10
    xy = np.linspace(-0,1000,100)
            
    # Period
    x = per_pos_m
    y = per_los_m

    xmin, xmax = 80, 600
    ymin, ymax = 80, 600
                
    #plt.close('all')
    fig = plt.figure(figsize = (4,4))
    bins = 50
    gs = gridspec.GridSpec(ncols = 5,nrows = 5)

    ax_joint = fig.add_subplot(gs[1:5,0:4])
    ax_marg_x = fig.add_subplot(gs[0,0:4])
    ax_marg_y = fig.add_subplot(gs[1:5,4])

    angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
    ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1.)#, linestyle = '--')
    ax_joint.text(xmax - 60,ymax - 55,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

    ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
    ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
    ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)
    
    # Turn off ticks and labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    ax_marg_x.get_yaxis().set_visible(False)
    ax_marg_y.get_xaxis().set_visible(False)

    # horizontal hist
    ax_marg_x.spines['top'].set_visible(False)
    ax_marg_x.spines['right'].set_visible(False)
    ax_marg_x.spines['left'].set_visible(False)
    # vertical hist
    ax_marg_y.spines['top'].set_visible(False)
    ax_marg_y.spines['right'].set_visible(False)
    ax_marg_y.spines['bottom'].set_visible(False)
    # joint plot
    ax_joint.spines['top'].set_visible(False)
    ax_joint.spines['right'].set_visible(False)
    
    # Set labels on joint
    ax_joint.set_xlabel(r' $\overline{P}_{\rm{POS}}$ [s]', fontdict = font)
    ax_joint.set_ylabel(r' $\overline{P}_{\rm{LOS}}$ [s]', fontdict = font)
    ax_joint.tick_params(labelsize = fontsize)
    ax_marg_x.tick_params(labelsize = fontsize)
    ax_marg_y.tick_params(labelsize = fontsize)
    
    # Set ax limits on marginals
    ax_joint.set_xlim([xmin,xmax])
    ax_joint.set_ylim([ymin,ymax])
    ax_marg_y.set_ylim([ymin,ymax])
    ax_marg_x.set_xlim([xmin,xmax])

    pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
    pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
    ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
    
    plt.tight_layout()
    plt.subplots_adjust(left = 0.152,
                        bottom = 0.135,
                        right = 0.98,
                        top = 0.98,
                        wspace = 0.,
                        hspace = 0.
    )
    #plt.show()
    filename = outdir + 'period_comp.pdf'
    fig.savefig(filename, dpi = 1000)
    plt.close(fig)


    # Amplitude
    
    x = A_pos_test
    y = A_los_test

    xmin, xmax = 0, 6
    ymin, ymax = 0, 6
                
    #plt.close('all')
    fig = plt.figure(figsize = (4,4))
    bins = 50
    gs = gridspec.GridSpec(ncols = 5,nrows = 5)

    ax_joint = fig.add_subplot(gs[1:5,0:4])
    ax_marg_x = fig.add_subplot(gs[0,0:4])
    ax_marg_y = fig.add_subplot(gs[1:5,4])

    angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
    ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linewidth = 1., linestyle = '--')#, linestyle = '--')
    ax_joint.text(xmax - 1.25,ymax - 0.5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

    ax_joint.scatter(x, y, alpha = 0.25, color = 'orangered')
    ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = 100)
    ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = 25)
    
    # Turn off ticks and labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    ax_marg_x.get_yaxis().set_visible(False)
    ax_marg_y.get_xaxis().set_visible(False)

    # horizontal hist
    ax_marg_x.spines['top'].set_visible(False)
    ax_marg_x.spines['right'].set_visible(False)
    ax_marg_x.spines['left'].set_visible(False)
    # vertical hist
    ax_marg_y.spines['top'].set_visible(False)
    ax_marg_y.spines['right'].set_visible(False)
    ax_marg_y.spines['bottom'].set_visible(False)
    # joint plot
    ax_joint.spines['top'].set_visible(False)
    ax_joint.spines['right'].set_visible(False)
    
    # Set labels on joint
    ax_joint.set_xlabel(r' $\overline{A}_{\rm{POS}}$ [km s$^{-1}$]', fontdict = font)
    ax_joint.set_ylabel(r' $\overline{A}_{\rm{LOS}}$ [km s$^{-1}$]', fontdict = font)
    ax_joint.tick_params(labelsize = fontsize)
    ax_marg_x.tick_params(labelsize = fontsize)
    ax_marg_y.tick_params(labelsize = fontsize)
    
    # Set ax limits on marginals
    ax_joint.set_xlim([xmin,xmax])
    ax_joint.set_ylim([ymin,ymax])
    ax_marg_y.set_ylim([ymin,ymax])
    ax_marg_x.set_xlim([xmin,xmax])

    pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
    pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
    ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
    
    plt.tight_layout()
    plt.subplots_adjust(left = 0.155,
                        bottom = 0.135,
                        right = 0.98,
                        top = 0.98,
                        wspace = 0.,
                        hspace = 0.
    )
    #plt.show()
    filename = outdir + 'Amp_comp.pdf'
    fig.savefig(filename, dpi = 1000)
    plt.close('fig')

    plt.show()




        
