# fibril_roi.py

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

# Functions
def length(x1, x2, y1,y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

# Font properties
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# Directories
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'
invdir = '/home/seki2695/INV/stic/II/roi/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcse
cad = 8
fr_select = [29,75,135,167] # frames to display
w_pos = 13 # wavelength position to display (K2V, K3, K2R) = (10,13,16)
w_pos8 = 12 # 8542 line center

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

# extracting the fibril and cuts information
fibdir = datadir+'fr29/'
fib_file = (file_search(fibdir,'crispex*3950*.csav'))[0]
fib = restore(fibdir+fib_file)

cut1_file = (file_search(fibdir,'crispex*3950*.csav'))[1]
cut2_file = (file_search(fibdir,'crispex*3950*.csav'))[2]
cut3_file = (file_search(fibdir,'crispex*3950*.csav'))[3]

cut1 = restore(fibdir+cut1_file)
cut2 = restore(fibdir+cut2_file)
cut3 = restore(fibdir+cut3_file)

xmin = int(np.min([np.min([cut1.x_coords, cut2.x_coords,
                           cut3.x_coords]),np.min(fib.x_coords)]) - 5)
ymin = int(np.min([np.min([cut1.y_coords, cut2.y_coords,
                           cut3.y_coords]),np.min(fib.y_coords)]) - 5)
xmax = int(np.max([np.max([cut1.x_coords, cut2.x_coords,
                           cut3.x_coords]),np.max(fib.x_coords)]) + 5)
ymax = int(np.max([np.max([cut1.y_coords, cut2.y_coords,
                           cut3.y_coords]),np.max(fib.y_coords)]) + 5)
#stop()
# Exxtracting the inversion results
fr0_file = file_search(invdir, '*atmosout*fr29.nc')[0]
fr1_file = file_search(invdir, '*atmosout*fr75.nc')[0]
fr2_file = file_search(invdir, '*atmosout*fr135.nc')[0]
fr3_file = file_search(invdir, '*atmosout*fr167.nc')[0]

fr0_inv = sp.model(invdir+fr0_file)
fr1_inv = sp.model(invdir+fr1_file)
fr2_inv = sp.model(invdir+fr2_file)
fr3_inv = sp.model(invdir+fr3_file)

savedir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/w_calibrated/'
ck0  = mf.readfits(savedir+'roi_obs3950_fr29'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ck1  = mf.readfits(savedir+'roi_obs3950_fr75'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ck2  = mf.readfits(savedir+'roi_obs3950_fr135'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ck3  = mf.readfits(savedir+'roi_obs3950_fr167'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6

ca0  = mf.readfits(savedir+'roi_obs8542_fr29'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ca1  = mf.readfits(savedir+'roi_obs8542_fr75'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ca2  = mf.readfits(savedir+'roi_obs8542_fr135'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6
ca3  = mf.readfits(savedir+'roi_obs8542_fr167'+'.fits')[1,:, ymin:ymax, xmin:xmax].squeeze()*1e6

# chosen  factors for plotting
tau0, tau1, tau2 = 52, 26, 21
ltau = fr1_inv.ltau[0,0,0,:]
crap = np.load('/storage/users/sepid/INV/stic/fov/crap.npy')[:,:,:,1]

cmaps = ['gray', 'gray', 'gray', 'inferno', 'inferno' , 'inferno' , 'bwr', 'bwr', 'bwr']*4
vmin = [ 13, 10., 2.,
         6., 4.45, 4.8,
         -4.5, -6.7, -7.5]*4
vmax = [31, 16, 3.9,
        7.15, 5.1, 5.95,
        +4.5, +6.7, +7.5]*4

cb_ticks = [[14, 21, 28],
            [11, 13, 15],
            [2.5, 3, 3.5],
            [6, 6.5, 7],
            [4.5, 4.75, 5],
            [4.5, 5, 5.5],
            [-4, 0, 4],
            [-6, 0, 6],
            [-7, 0, 7]]
ann = [r'cont 4000 $\rm{\AA}$', r'Ca II 8542 $\rm{\AA}$', r'Ca II K',
       r'log $\tau$ = '+str(np.round(ltau[tau0], decimals = 1)),
       r'log $\tau$ = '+str(np.round(ltau[tau1], decimals = 1)),
       r'log $\tau$ = '+str(np.round(ltau[tau2], decimals = 1)),
       r'log $\tau$ = '+str(np.round(ltau[tau0], decimals = 1)),
       r'log $\tau$ = '+str(np.round(ltau[tau1], decimals = 1)),
       r'log $\tau$ = '+str(np.round(ltau[tau2], decimals = 1)),
]

maps = [ck0[-1,:,:],
        ca0[w_pos8,:,:],
        ck0[w_pos,:,:],
        fr0_inv.temp[0,:,:,tau0].squeeze()*1e-3,
        fr0_inv.temp[0,:,:,tau1].squeeze()*1e-3,
        fr0_inv.temp[0,:,:,tau2].squeeze()*1e-3,
        fr0_inv.vlos[0,:,:,tau0].squeeze()*1e-5,# - crap[ymin:ymax,xmin:xmax,tau0],
        fr0_inv.vlos[0,:,:,tau1].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau1],
        fr0_inv.vlos[0,:,:,tau2].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau2],
        ck1[-1,:,:],
        ca1[w_pos8,:,:],
        ck1[w_pos,:,:],
        fr1_inv.temp[0,:,:,tau0].squeeze()*1e-3,
        fr1_inv.temp[0,:,:,tau1].squeeze()*1e-3,
        fr1_inv.temp[0,:,:,tau2].squeeze()*1e-3,
        fr1_inv.vlos[0,:,:,tau0].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau0],
        fr1_inv.vlos[0,:,:,tau1].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau1],
        fr1_inv.vlos[0,:,:,tau2].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau2],
        ck2[-1,:,:],
        ca2[w_pos8,:,:],
        ck2[w_pos,:,:],
        fr2_inv.temp[0,:,:,tau0].squeeze()*1e-3,
        fr2_inv.temp[0,:,:,tau1].squeeze()*1e-3,
        fr2_inv.temp[0,:,:,tau2].squeeze()*1e-3,
        fr2_inv.vlos[0,:,:,tau0].squeeze()*1e-5,# - crap[ymin:ymax,xmin:xmax,tau0],
        fr2_inv.vlos[0,:,:,tau1].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau1],
        fr2_inv.vlos[0,:,:,tau2].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau2],
        ck3[-1,:,:],
        ca3[w_pos8,:,:],
        ck3[w_pos,:,:],
        fr3_inv.temp[0,:,:,tau0].squeeze()*1e-3,
        fr3_inv.temp[0,:,:,tau1].squeeze()*1e-3,
        fr3_inv.temp[0,:,:,tau2].squeeze()*1e-3,
        fr3_inv.vlos[0,:,:,tau0].squeeze()*1e-5,# - crap[ymin:ymax,xmin:xmax,tau0],
        fr3_inv.vlos[0,:,:,tau1].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau1],
        fr3_inv.vlos[0,:,:,tau2].squeeze()*1e-5 - crap[ymin:ymax,xmin:xmax,tau2],
]

# plot settings
xfig, yfig = 8,5.3
plt.close('all')
f = plt.figure(figsize = (xfig, yfig))

nrow, ncol = 4,9
axes = gridspec.GridSpec(nrow, ncol) # grid scale of the left col.
axes.update(left=0.045,
               right=0.995,
               wspace=0.1,
               bottom=0.06,
               top=0.91,
               hspace = 0.05
    )
    

# axis info
rr = 40/1.5
ytick_pos = np.arange(0,(np.round((ymax-ymin)/rr)+1)*rr,rr)
ytick_lab = np.intc(np.round(ytick_pos*res, decimals = 0))
xtick_pos = np.arange(0,(np.round((xmax-xmin)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))

cuts_file = file_search('/home/seki2695/CODE/py/II/', '*npy')
fr1 = np.zeros((3,2))
for j in range(ncol):
    for i in range(nrow):
        ax = plt.subplot(axes[i,j],adjustable = 'box')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks(xtick_pos)
        ax.set_yticks(ytick_pos)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        if(j==0):
            ax.set_yticklabels(ytick_lab, fontdict = font)
            if(i==3):
                ax.set_xlabel(r'x [arcsec]', fontdict = font, labelpad = 0)
                ax.set_ylabel(r'y [arcsec]', fontdict = font)
        if(i==3):
            ax.set_xticklabels(ytick_lab, fontdict = font)
        im = ax.imshow(maps[i*ncol+j], cmap = cmaps[i*ncol+j], origin = 'lower',
                  vmin = vmin[i*ncol+j],
                  vmax = vmax[i*ncol+j]
        )

        if(j==0):
            props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='none')
            ax.text(3, 77, str(np.around(fr_select[i]*8./60, decimals = 1))+' min', color = 'black', fontdict = font, bbox=props)

        if(i==0):
            axin = inset_axes(ax,
                               width="90%",  # width = 10% of parent_bbox width
                               height="5%",  # height : 50%
                               loc='lower center',
                               bbox_to_anchor=(0, 1.03, 1, 1),
                               bbox_transform=ax.transAxes,
                               borderpad=0,
            )

            cb_tick_labels = [str((cb_ticks[j])[0]),str((cb_ticks[j])[1]),str((cb_ticks[j])[2])]
            cb = plt.colorbar(im, cax=axin, orientation="horizontal", ticks = cb_ticks[j])
            cb.ax.set_xticklabels(cb_tick_labels, fontdict = font)
            cb.ax.xaxis.tick_top()
            cb.ax.xaxis.set_tick_params(pad = 0)
            if(j==1):
                cb.ax.set_title(r'I [$ \times $10$^{3}$ W m$^{-2}$ Hz$^{-1}$ sr$^{-1}$]', fontdict = font, pad = 18)
            if(j==4):
                cb.ax.set_title(r'T [kK]', fontdict = font, pad = 18)
            if(j==7):
                cb.ax.set_title(r'$v_{\rm{LOS}}$ [km s$^{-1}$]', fontdict = font, pad = 18)

            if(j>5):
                ax.text(1, 1, ann[j], fontdict = font, color = 'black')
            elif(j==3):
                ax.text(1, 1, ann[j], fontdict = font, color = 'black')
                ax.text(1, 1, ann[j], fontdict = font, color = 'white')                
            else:
                ax.text(1, 1, ann[j], fontdict = font, color = 'white')
        #L1 = length(cut1.x_coords[-1], cut1.x_coords[0], cut1.y_coords[-1], cut1.y_coords[0])
        #L2 = length(cut2.x_coords[-1], cut2.x_coords[0], cut2.y_coords[-1], cut2.y_coords[0])
        #L3 = length(cut3.x_coords[-1], cut3.x_coords[0], cut3.y_coords[-1], cut3.y_coords[0])
        #for l in range(len(cuts_file)):
          #  fr1 = np.load(cuts_file[l])
           # cut_loc = np.load()
        x01 = cut3.x_coords[0] + np.load(cuts_file[0])[i]*np.abs((cut3.x_coords[0]-cut3.x_coords[-1]))
        y01 = cut3.y_coords[0] + np.load(cuts_file[0])[i]*np.abs((cut3.y_coords[0]-cut3.y_coords[-1]))
        x02 = cut2.x_coords[0] + np.load(cuts_file[1])[i]*np.abs((cut2.x_coords[0]-cut2.x_coords[-1]))
        y02 = cut2.y_coords[0] + np.load(cuts_file[1])[i]*np.abs((cut2.y_coords[0]-cut2.y_coords[-1]))
        x03 = cut1.x_coords[0] + np.load(cuts_file[2])[i]*np.abs((cut1.x_coords[0]-cut1.x_coords[-1]))
        y03 = cut1.y_coords[0] + np.load(cuts_file[2])[i]*np.abs((cut1.y_coords[0]-cut1.y_coords[-1]))

        r1 = np.sqrt((x02-x01)**2+(y02-y01)**2)
        r2 = np.sqrt((x02-x03)**2+(y02-y03)**2)
        #print(np.mean([r1/(8*dphi.shape
          #                 (71-68)),r2/(8*(79-71))])*res*725)
        #print(np.mean([r1/(8*(123-119)),r2/(8*(129-123))])*res*725)

        if(j==2):
            #ax.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'orangered', alpha = 0.5)
            ax.plot(cut1.x_loop_pts[4:-1]-xmin, cut1.y_loop_pts[4:-1]-ymin, color = 'orangered', linestyle = '--', alpha = 0.65,linewidth = 1)
            ax.plot(cut2.x_loop_pts[4:-1]-xmin, cut2.y_loop_pts[4:-1]-ymin, color = 'orangered', linestyle = '--', alpha = 0.65,linewidth = 1)
            ax.plot(cut3.x_loop_pts[4:-1]-xmin, cut3.y_loop_pts[4:-1]-ymin, color = 'orangered', linestyle = '--', alpha = 0.65,linewidth = 1)
            if(i==0):
                ax.text(cut1.x_coords[0]-xmin-4, cut1.y_coords[0]-ymin+6, '(3)', fontdict=font, color = 'orangered')
                ax.text(cut2.x_coords[0]-xmin-4, cut2.y_coords[0]-ymin+6, '(2)', fontdict=font, color = 'orangered')
                ax.text(cut3.x_coords[0]-xmin-4, cut3.y_coords[0]-ymin+6, '(1)', fontdict=font, color = 'orangered')
                #else:
                #ax.scatter(x = np.array([x01, x02, x03])-xmin, y = np.array([y01, y02, y03])-ymin, color = 'orangered', s = 5, alpha = 0.7)
                #xint = np.linspace(np.int(np.round(x01)),np.int(np.round(x03)),np.int(np.round(x03))-np.int(np.round(x01))+1)
                #yint = np.interp(xint, np.array([x01, x02, x03]), np.array([y01, y02, y03]))
                #ax.plot(xint-xmin, yint-ymin, alpha = 0.5)

        #if(j==8):
         #   if(i==3):
          #      print('')
            #else:
                #ax.scatter(x = np.array([x01, x02, x03])-xmin, y = np.array([y01, y02, y03])-ymin, color = 'black', s = 5, alpha = 0.5)
        
        if(i==0 and j==5):
            #ax.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'white', alpha = 0.5)
            ax.plot(cut1.x_loop_pts[4:-1]-xmin, cut1.y_loop_pts[4:-1]-ymin, color = 'white', linestyle = '--', alpha = 0.5,linewidth = 1)
            ax.plot(cut2.x_loop_pts[4:-1]-xmin, cut2.y_loop_pts[4:-1]-ymin, color = 'white', linestyle = '--', alpha = 0.5,linewidth = 1)
            ax.plot(cut3.x_loop_pts[4:-1]-xmin, cut3.y_loop_pts[4:-1]-ymin, color = 'white', linestyle = '--', alpha = 0.5,linewidth = 1)

        if(i==0 and j==8):
            #ax.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'black', alpha = 0.5)
            ax.plot(cut1.x_loop_pts[4:-1]-xmin, cut1.y_loop_pts[4:-1]-ymin, color = 'black', linestyle = '--', alpha = 0.5,linewidth = 1)
            ax.plot(cut2.x_loop_pts[4:-1]-xmin, cut2.y_loop_pts[4:-1]-ymin, color = 'black', linestyle = '--', alpha = 0.5,linewidth = 1)
            ax.plot(cut3.x_loop_pts[4:-1]-xmin, cut3.y_loop_pts[4:-1]-ymin, color = 'black', linestyle = '--', alpha = 0.5,linewidth = 1)

plt.show()

fname = outdir + 'fibril_roi.pdf'
f.savefig(fname, quality = 100)
