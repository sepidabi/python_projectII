# fibril_hor_oscillation.py

'''Meant for plotting a sub field of view
with a fibril and its horizontal oscillations
(the cuts along the fibrils are already defined via crispex)
'''

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
from skimage import exposure
from scipy.fftpack import fft2, ifft2


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

plt.close("all")

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcse
cad = 8
fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)

#functions
def test_func(x, a, b, c, d):
    return a * np.sin(b * x - c)+d

#boxcar smoothing
def smooth(y, w):
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

# user contribution
for ii in range(len(file_search(fibdir,'crispex*3950*.csav'))):
    print(ii, file_search(fibdir,'crispex*3950*.csav')[ii])
cc = input("Enter the number of the cut: ")

cut_file = (file_search(fibdir,'crispex*3950*.csav'))[cc]
cut = restore(fibdir+cut_file)
cube_raw = np.mean(cut.loop_slab[w_pos-3:w_pos+3,:,:],axis = 0).squeeze()*1e9

# background intensity from the wings
wg = 0.15
cube_bg = wg*np.mean(cut.loop_slab[:w_pos-6,:,:]+cut.loop_slab[w_pos+8:,:,:],axis=0).squeeze()*1e9
    
cube_trunc = cube_raw - cube_bg # removing the background intensity

# contrast correction
cube_med = exposure.rescale_intensity(sgnl.medfilt2d(cube_trunc,kernel_size = [3,1]), out_range=(0, 1.))

cube_sharp = exposure.rescale_intensity(sharpen(cube_trunc, sigma =[3,1], unsigma = [1,3]), out_range=(0, 1.))

cube = exposure.rescale_intensity(exposure.adjust_gamma(cube_sharp,0.1), out_range=(0, 1.))

cube_std = np.std(cube)
fact = 1.5
vmin, vmax = np.mean(cube) - fact*cube_std, np.mean(cube) + fact*cube_std

xx = cut.loop_slab.shape[2]

# to test the curve fitting
ti = np.arange(0, nt, 1.)

# plotting elements
tt = 30
ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

rr = 40/1.5
xtick_pos = np.arange(0,(np.round((xx)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

# Extracting the inversion results
inv_inc=raw_input("include inversion results [default=0 or 1]? ").strip() or 0
if(inv_inc):
    inv_res = sp.model(invdir+
                       file_search(invdir, "*atmosout*"+cut_file[-11:-5]+"*nc")[0])
    
    temp_cube = inv_res.temp[0]
    vlos_cube = inv_res.vlos[0]
    vturb = inv_res.vturb[0]
    ltau = inv_res.ltau[0,0,0,:]
    ntau = ltau.size
    
    for l in range(ntau):
        vlos_cube[:,:,l] = (vlos_cube[:,:,l]
        - spnd.filters.gaussian_filter(inv_res.vlos[0,:,:,l], [65,65], mode = 'constant')).squeeze()*1e-5
        
    tau_i = 23
    
# 3D view of the cuts
if(0):
    plt.close("all")
    thd = plt.figure()
    ax = thd.add_subplot(111, projection='3d')
    y = np.arange(0, cube.shape[0], 1)
    x = np.arange(0, cube.shape[1], 1)
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y,
                    cube,
                    cmap = "gray",
                    alpha = 0.75,
                    vmin=vmin,
                    vmax=vmax,
    )
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

    
# tracking mode
track = raw_input("Tracking on? [0 or default=1]: ").strip() or 1

if(track==1):
    ########################
    # CURVE FITTING
    ########################
    cut_fname = outdir+cut_file[-11:-5]+'.txt'

    yi = np.argmax(cube, axis = 1)
    
    #plt.close('all')
    f, ax1 = plt.subplots(figsize = (5.75,14))
    ax1.set_xlabel(r'length [pixel]', fontdict = font)
    ax1.set_ylabel(r'$t$ [fr]', fontdict = font)
    ax1.set_ylim(0,nt-1)
    
    ax1.imshow(cube,
               cmap = 'gray', origin = 'lower', aspect=10./6.,
               #interpolation = 'bilinear',
                           vmin=vmin,
                          vmax = vmax,
               #        aspect = 1.5,
    )
    for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):
        osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
        coords = np.loadtxt(osc_fname)
        ax1.plot(coords[:,0], coords[:,1], color = 'orange')
        [x,y] = coords[-1,:]
        ax1.text(x,y,osc_fname[-10:-4], color = 'orangered', fontweight='bold', rotation='vertical')

        ax1.plot()
    plt.tight_layout()
    
    coord = []
    
    # GUI stuff
    cursor = Cursor(ax1, horizOn=True, vertOn=True, useblit=True,
                    color = 'r', linewidth = 0.75, alpha = 0.75)
    
    def onclick(event):
        xx, tt = event.xdata, event.ydata
        global coord
        coord.append((xx, tt))
        ax1.plot(xx,tt, '+', color = 'red')
        print('a(t)=%1.2f, t=%1.2f' %
              (np.round(xx, decimals = 2),
               np.round(tt, decimals = 2)))
        f.canvas.draw() #redraw the figure
        
    c = 1
    while c==1:
        cid = f.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        

        stop()
        
        # save the clicked points
        now = datetime.now()
        osc_fname = outdir+"cut"+cut_file[-11:-5]+"_osc"+now.strftime("%y%m%d-%H%M%S")+".txt"        
        np.savetxt(osc_fname, coord, fmt='%3.8f', delimiter=' ', newline='\n', header='a(t) [pixel], t [fr]', footer='', comments='# ', encoding=None)
        coord = []

        c = raw_input("Still tracking? [0 or default=1]: ").strip() or 1
        
        if c==0:
            break

# intensity oscillation PLOTs
plt.close('all')
if(inv_inc==1):
    fig = plt.figure(figsize = (4,5.8))
    gs = gridspec.GridSpec(3,1)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[2,0])

else:
    fig = plt.figure(figsize = (8,4.5))
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])
    
for oo in range(len(file_search(outdir, "cut"+cut_file[-11:-5]+"*txt"))):

    osc_fname = outdir+file_search(outdir, "cut"+cut_file[-11:-5]+"*txt")[oo]
    coord = np.loadtxt(osc_fname)
    amin, amax = int(np.min(coord[:,0])), int(np.max(coord[:,0]))+1
    tmin, tmax = int(np.min(coord[:,1])), int(np.max(coord[:,1]))+1
    trange = np.linspace(tmin, tmax, tmax-tmin+1)
    
    # interpolating the points in between the clicked fit
    y_interp = np.interp(trange, coord[:,1], coord[:,0])
    osc_fname_interp = outdir+"interp"+osc_fname[31:]
    np.savetxt(osc_fname_interp, y_interp, fmt='%3.8f', encoding=None)
    
    # track max intensity based on the interp click_fit
    y_imax = np.zeros(tmax-tmin+1)
    dist = 1. # width of the interval in which the max intensity is taken
    for i in range(tmax-tmin+1):
        y_imax[i] = np.argmax(cube[i, int(y_interp[i]-dist):int(y_interp[i]+dist)])+int(y_interp[i])-dist

    smth_y = smooth(y_imax, 3)
    osc_fname_smth = outdir+"smth"+osc_fname[31:]
    np.savetxt(osc_fname_smth, smth_y, fmt='%3.8f', encoding=None)

    # Calculation of the velocity of the oscillation in the plane of sky
    xp = trange*cad
    fxp = smth_y*res*725
    der = derivative(xp, fxp, dx = 1e-5)
    x, fx, f1x = der[0]-np.min(der[0]), der[1], der[2]

    # Intensity
    #------------
    im_cak = ax1.imshow(np.transpose(cube),
               cmap = 'gray', origin = 'lower', #aspect=ratio
               #interpolation = 'bilinear',
               vmin=vmin,
               vmax = vmax,
               aspect = 0.6,
    )
    ax1.plot(coord[:,1], coord[:,0],color = 'wheat', alpha = 0.5, label = 'Manual fit', linestyle = '--', linewidth = 1)
    ax1.plot(trange, smth_y, color = 'orange', alpha = 0.5, label = r'smthd I$_{\rm{max}}$', linewidth = 2)
    [y,x] = coord[0,:]
    ax1.text(x,y,osc_fname[-10:-4], color = 'orange', fontdict = font)
    plt.show()
    #stop()

    if(inv_inc==1):
        # temp
        #--------
        tau = 23
        im_temp = ax2.imshow(np.transpose(temp_cube[:,:,tau]).squeeze()*1e-3,
                             vmin = 3.7, vmax = 4.8,
                             cmap = 'inferno', origin = 'lower', #aspect=ratio
                             #interpolation = 'bilinear',
                             #vmin=2.e-9,
                             #vmax = 3.5e-9,
                             aspect = 0.6,
                             )
                    
        # v_los
        #--------
        vlos_lim = np.round(np.max(np.abs(vlos_cube[:,:,tau])), decimals = 1)
        im_vlos = ax3.imshow(np.transpose(vlos_cube[:,:,tau]),
                             vmin = -vlos_lim, vmax = vlos_lim,
                             cmap = 'bwr', origin = 'lower', #aspect=ratio
                             #interpolation = 'bilinear',
                             #vmin=2.e-9,
                             #vmax = 3.5e-9,
                             aspect = 0.6,
                             )
        
           
        temp_smth = np.zeros((tmax-tmin+1,ntau))
        vlos_smth = np.zeros((tmax-tmin+1,ntau))

        for i in range(tmax-tmin+1):
            for j in range(ntau):
                temp_smth[i,j] = temp_cube[i,int(np.rint(smth_y[i])),j]
                vlos_smth[i,j] = vlos_cube[i,int(np.rint(smth_y[i])),j]

    # smooth fit to the oscillations
    Imax_smth = np.zeros(tmax-tmin+1)
    for i in range(tmax-tmin+1):
        Imax_smth[i] = cube[i,int(np.rint(smth_y[i]))]


    
    # plotting settings
    if (0):
        alph = 0.25
        sm_fact = 8
        per = vlos_smth.shape[0]
        brief_fact = f1x.shape[0]/per
        tt_osc = 10
        ttick_pos_osc = np.arange(0,(np.round((vlos_smth.shape[0])/tt_osc)+1)*tt_osc,tt_osc)
        ttick_lab_osc = ttick_pos_osc*cad
        pig = plt.figure(figsize = (8,4))
        gs = gridspec.GridSpec(1,1)
        ax = plt.subplot(gs[0,0])
        ax.set_title("no. "+ str(oo)+" - "+osc_fname[31:-4])
        ax.set_xticks(ttick_pos_osc)
        ax.set_xticklabels(ttick_lab_osc)
        ax.set_xlabel("t [s]")
        ax.xaxis.set_minor_locator(AutoMinorLocator(8))
        ax.set_ylabel(r"$v$ [km s$^{-1}$]")
        
        x_brief = np.zeros(per)
        f1x_brief = np.zeros(per)
        
        for iii in range(per):
            x_brief[iii] = x[iii*brief_fact]
            f1x_brief[iii] = f1x[iii*brief_fact]
        
    # plot of LOS vs POS
    #-------------------------
    if (0):
        ax.plot(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i])*1e-5, color = 'orangered', linestyle = "--", alpha = alph)
        ax.plot(smooth(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i])*1e-5,sm_fact), label = r"LOS", color = 'orangered')
        ax.plot(x_brief/cad, f1x_brief-np.mean(f1x_brief), linestyle = "--", alpha = alph, color = "darkgray")
        ax.plot(x_brief/cad, smooth(f1x_brief-np.mean(f1x_brief),sm_fact), label = r"POS", color = 'darkgray')
        plt.legend()
        plt.tight_layout()
        osc_fname_plot = osc_fname[:-4]+".pdf"
        pig.savefig(osc_fname_plot, quality = 100)
        plt.close(pig)# intensity oscillation PLOTs

    # velocity vector 3D
    #------------------------
    if(0):
        y = smooth(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5), sm_fact)
        x = smooth(f1x_brief - np.mean(f1x_brief), sm_fact)
        t = x_brief/cad
        ratio = len(t)/(6*np.max(np.abs([x,y])))
        
        f3 = plt.figure()#figsize = (4,3*ratio))
        # Total velocity 3d
        ax = f3.add_subplot(1, 1, 1, projection='3d')#, aspect = ratio)
        ax.azim = 40
        ax.elev = 10
        norm_fact =np.max([np.max(np.abs(vlos_smth[:,tau_i]*1e-5 - np.mean(vlos_smth[:,tau_i]*1e-5))), np.max(np.abs(f1x_brief - np.mean(f1x_brief)))])
        ax.set_zticks(ttick_pos_osc)
        ax.set_zticklabels(ttick_lab_osc)
        ax.zaxis.set_minor_locator(AutoMinorLocator(8))
        ax.set_zlabel(r"$t$ [s]")
        ax.set_ylabel(r"[km s$^{-1}$]")
        #ax.set_ylim(-norm_fact, norm_fact)
        ax.set_xlabel(r"[km s$^{-1}$]")
        #ax.set_xlim(-norm_fact, norm_fact)
        ax.plot3D(x, y, t, 'black', label = r'$\mathit{\mathbf{v}}$(t)', alpha = 0.75, linewidth = 2.5)
        ax.plot3D(np.zeros(len(x))+np.min(x), y, t, 'orangered', alpha = 0.5, label = r'$v_{\rm{LOS}}$(t)')
        ax.plot3D(x, np.zeros(len(y))+np.min(y), t, 'rosybrown', alpha = 0.5, label = r'$v_{\rm{POS}}$(t)')
        #ax.plot3D(x, y, np.zeros(len(t)), 'darkgray', alpha = 0.5)
        #ax.set_aspect('equal', 'box')
        plt.legend(loc = 'upper right')
        plt.tight_layout()
        plt.show()
        stop()

        f3.savefig(outdir +osc_fname[31:-4]+"-"+ "no"+ str(oo)+".pdf", quality = 100)
        plt.close(f3)
        
    # cross_correlation plot
    #----------------------------
    if(0):
        #plt.close("all")
        fc = plt.figure(figsize = (4,2))
        #fc.suptitle("no. "+ str(oo)+" - "+osc_fname[31:-4])
        ax = fc.add_subplot(1, 1, 1)#, aspect = 20.)
        ttick_ncc = np.arange(0,(np.round((vlos_smth.shape[0]*2)/tt_osc)+1)*tt_osc,tt_osc) - np.median(np.arange(0,(np.round((vlos_smth.shape[0]*2)/tt_osc)+1)*tt_osc,tt_osc))
        ttick_pos_ncc = np.linspace(np.min(ttick_ncc)+tt_osc, np.max(ttick_ncc)-tt_osc,9)
        test = ttick_pos_ncc*cad
        ttick_lab_ncc = test.astype("int")
        ax.set_xticks(ttick_pos_ncc)
        ax.set_xticklabels(ttick_lab_ncc, fontdict = font)
        x_norm = (x-np.mean(x)) / (np.std(x)*len(x))
        y_norm = (y-np.mean(y)) / (np.std(y))
        ncc = sgnl.correlate(x_norm, y_norm, 'full')
        lag = np.arange(-len(y_norm)+1, len(x_norm))
        ax.axhline(0, linestyle = "--", color = "black", linewidth = 0.75)
        ax.axvline(0, linestyle = "--", color = "black", linewidth = 0.75)
        ax.plot(lag, ncc, label = r"NCC", color = "red")
        ax.set_ylim(-1,1)
        ax.set_xlim(np.min(ttick_ncc), np.max(ttick_ncc))
        ax.tick_params(axis = 'y', labelsize = 8)
        ax.set_xlabel(r"lag", fontdict = font)
        ax.set_ylabel(r"NCC", fontdict = font)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.show()
        
        
        plt.tight_layout()
        fc.savefig(outdir +osc_fname[31:-4]+"_NCC-"+ "no"+ str(oo)+".pdf", quality = 100)
        plt.close(fc)

if(inv_inc==1):
    ax1.set_ylabel(r'length [arcsec]', fontdict = font)
    ax1.set_xticks(ttick_pos)
    ax1.set_xticklabels([])
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_xlabel('')
    ax1.set_yticks(xtick_pos)
    ax1.set_yticklabels(xtick_lab, fontdict = font)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_xlim(0,nt-1)    
    ax1.set_ylim(0,xx-1)
    
    ax2.set_ylabel(r'length [arcsec]', fontdict = font)
    ax2.set_xticks(ttick_pos)
    ax2.set_xticklabels([])
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_xlabel('')
    ax2.set_yticks(xtick_pos)
    ax2.set_yticklabels(xtick_lab, fontdict = font)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.set_xlim(0,nt-1)
    ax2.set_ylim(0,xx-1)
    
    ax3.set_ylabel(r'length [arcsec]', fontdict = font)
    ax3.set_xticks(ttick_pos)
    ax3.set_xticklabels(ttick_lab, fontdict = font)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xlabel('t [min]', fontdict = font)
    ax3.set_yticks(xtick_pos)
    ax3.set_yticklabels(xtick_lab, fontdict = font)
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.set_xlim(0,nt-1)
    ax3.set_ylim(0,xx-1)

    axins1 = inset_axes(ax1,
                        width="3%",  
                        height="90%",
                        loc='center left',
                        bbox_to_anchor=(1.01, 0., 1, 1),
                        bbox_transform=ax1.transAxes,
                        borderpad=0,
                        )
    
    
    axins2 = inset_axes(ax2,
                        width="3%",  
                        height="90%",
                        loc='center left',
                        bbox_to_anchor=(1.01, 0., 1, 1),
                        bbox_transform=ax2.transAxes,
                        borderpad=0,
                        )
    
    axins3 = inset_axes(ax3,
                        width="3%",  
                        height="90%",
                        loc='center left',
                        bbox_to_anchor=(1.01, 0., 1, 1),
                        bbox_transform=ax3.transAxes,
                        borderpad=0,
                        )
    
    cbar1 = fig.colorbar(im_cak, cax=axins1)
    cbar2 = fig.colorbar(im_temp, cax=axins2)
    cbar3 = fig.colorbar(im_vlos, cax=axins3)
    
    cbar1.set_label(r'I [$\times 10^{-6}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$ sr$^{-1}$]', size = 8.)
    cbar2.set_label(r'$v_{\rm{LOS}}$ [km s$^{-1}$]', size = 8.)
    cbar3.set_label(r'$T$ [kK]', size = 8.)
    
    cbar1.ax.tick_params(axis='y', labelsize = 8.)
    cbar2.ax.tick_params(axis='y', labelsize = 8.)
    cbar3.ax.tick_params(axis='y', labelsize = 8.)
    filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_non.pdf'

else:
    ax1.set_ylabel(r'length [arcsec]', fontdict = font)
    ax1.set_xticks(ttick_pos)
    ax1.set_xticklabels(ttick_lab, fontdict = font)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_xlabel('t [min]', fontdict = font)
    ax1.set_yticks(xtick_pos)
    ax1.set_yticklabels(xtick_lab, fontdict = font)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.set_xlim(0,nt-1)
    ax1.set_ylim(0,xx-1)
    filename = outdir + 'oscillation_curvefit_' +cut_file[-11:-5]+'_test.pdf'

gs.update(left=0.07,
          right=0.87,
          wspace=0.0,
          bottom=0.07,
          top=0.995,
          hspace = 0.05,
)

plt.show()
fig.savefig(filename, quality = 100)
print 'file saved to: ' + filename


