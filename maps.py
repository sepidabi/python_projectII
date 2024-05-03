# A brief script to generate cool images
# for the thesis
# from chromospheric observation at SST
# written by Sepideh Kianfar

import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns
import scipy.ndimage as spnd
import matplotlib
from matplotlib.gridspec import GridSpec

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }


datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'
fibdir = datadir+'big_cut/'

w_cak = 11
w_fe = 7
w_ca8 = 9
fr = 59

# data extract
res = 0.0375 # CHROMIS pixel size in arcsec
pref = ['6173','8542','3950']#,'6563']
pref_name = ['fe','ca8','cak']#,'ha']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])[fr]
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])[fr]
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])[fr]
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fr]
cak_cont = cube_cak[0,-1,:,:]
fe_cont = cube_fe[0,-1,:,:]

# RoI coords
x1, x2, y1, y2 = 425, 900, 400, 1110

# photosphere LPtot
LPtot = np.zeros((cube_fe.shape[2], cube_fe.shape[3]))
for ii in range(cube_fe.shape[1]-1):
    LPtot += np.sqrt(cube_fe[1,ii,:,:]**2 + cube_fe[2,ii,:,:]**2)/fe_cont
LPtot = LPtot/cube_fe.shape[1]

# photosphere CPtot
CPtot = unsharp(0.5*(np.mean(cube_fe[3,0:7,:,:], axis=0) - np.mean(cube_fe[3,7:-1,:,:], axis=0))/fe_cont, alpha = 0.5, sigma = 2.)

# CaK wavelength-summed and unsharped
wi, wf = 10, 16
cak_sum_un = unsharp(np.mean(cube_cak[0, wi:wf+1, :, :],axis=0), alpha = 0.5, sigma = 10.)

# CaK Core
cak_wc = 13
cak_core = unsharp(cube_cak[0, cak_wc, :, :], alpha = 0.5, sigma = 5.)

# Ca 8542 Core
ca8_wc = 11
ca8_core = unsharp(cube_ca8[0,ca8_wc,:,:], alpha = 0.5, sigma = 5.)

# maps
plt.close('all')

f = plt.figure(figsize = [8.,5.9])
gs = gridspec.GridSpec(4,4) # grid scale of the 1st col.
gs.update(left=0.06,
          right=0.98,
          wspace=0.08,
          bottom=0.07,
          top=1.,
          hspace = 0.01,
)

text = ['', r'a) continuum 4000 $\mathrm{\AA}$', r'b) Ca II K', r'c) Fe I 6173 LP', r'd) Fe I 6173 CP', r'e) Ca II 8542 $\mathrm{\AA}$ core', r'f) Ca II K core']
textcolor = ['','white', 'white', 'white', 'black', 'white', 'white']

xx, yy = cube_fe.shape[3], cube_fe.shape[2]

# axis info
sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact)+1)*sc_fact,sc_fact)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact)+1)*sc_fact,sc_fact)
xtick_lab = np.round(xtick_pos*res).astype(int)

ax1 = plt.subplot(gs[0:2,0:2],adjustable = 'box')
ax1.imshow(cak_cont, origin = 'lower', cmap = 'gray')
alp = 1.
lstyle = '--'
lw = 0.75
lcolor = 'white'
ax1.plot([x1,x1],[y1,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax1.plot([x1,x2],[y1,y1], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax1.plot([x2,x2],[y1,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax1.plot([x1,x2],[y2,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax1.set_xticks(xtick_pos)
ax1.set_xlim(0,xx)
ax1.set_ylim(0,yy)
ax1.set_yticks(ytick_pos)
ax1.set_ylabel(r'y [arcsec]', fontdict=font)
ax1.set_yticklabels(ytick_lab, fontdict = font)
ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
ax1.yaxis.set_minor_locator(AutoMinorLocator(10))
ax1.tick_params(which='minor', length=2)
ax1.set_xticks(xtick_pos)
ax1.set_xticklabels(xtick_lab, fontdict = font)
ax1.text(10, 10 , text[1], fontdict = font, color = textcolor[1])

ax2 = plt.subplot(gs[0:2,2:],adjustable = 'box')
ax2.imshow(cak_sum_un, origin = 'lower', cmap = 'gray')
# ROI coords
xmax = 520
xmin = 456
ymax = 885
ymin = 797
ax2.plot([x1,x1],[y1,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax2.plot([x1,x2],[y1,y1], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax2.plot([x2,x2],[y1,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax2.plot([x1,x2],[y2,y2], color = lcolor, linestyle = lstyle, linewidth = lw, alpha = alp)
ax2.plot([xmin,xmin],[ymin,ymax], color = "red", linewidth = 0.5, alpha = alp)
ax2.plot([xmin,xmax],[ymin,ymin], color = "red", linewidth = 0.5, alpha = alp)
ax2.plot([xmax,xmax],[ymin,ymax], color = "red", linewidth = 0.5, alpha = alp)
ax2.plot([xmin,xmax],[ymax,ymax], color = "red", linewidth = 0.5, alpha = alp)
ax2.set_xticks(xtick_pos)
ax2.set_xlim(0,xx)
ax2.set_ylim(0,yy)
ax2.set_yticks(ytick_pos)
ax2.set_ylabel('')
ax2.set_yticklabels([])
ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
ax2.yaxis.set_minor_locator(AutoMinorLocator(10))
ax2.tick_params(which='minor', length=2)
ax2.set_xticks(xtick_pos)
ax2.set_xticklabels(xtick_lab, fontdict = font)
ax2.text(10, 10 , text[2], fontdict = font, color = textcolor[2])

xx, yy = x2-x1, y2-y1

# axis info
sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact))*sc_fact,sc_fact/2.)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact))*sc_fact,sc_fact/2.)
xtick_lab = np.round(xtick_pos*res).astype(int)

ax3 = plt.subplot(gs[2:,0],adjustable = 'box')
ax3.imshow(LPtot[y1:y2,x1:x2], origin = 'lower', cmap = 'gray')
ax3.set_xticks(xtick_pos)
ax3.set_xlim(0,xx)
ax3.set_ylim(0,yy)
ax3.set_yticks(ytick_pos)
ax3.set_ylabel(r'y [arcsec]', fontdict=font)
ax3.set_yticklabels(ytick_lab, fontdict = font)
ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
ax3.tick_params(which='minor', length=2)
ax3.set_xticks(xtick_pos)
ax3.set_xticklabels(xtick_lab, fontdict = font)
ax3.set_xlabel(r'x [arcsec]', fontdict=font)
ax3.text(10, 10 , text[3], fontdict = font, color = textcolor[3])

ax4 = plt.subplot(gs[2:,1],adjustable = 'box')
ax4.imshow(CPtot[y1:y2,x1:x2], origin = 'lower', cmap = 'gray')
ax4.set_xticks(xtick_pos)
ax4.set_xlim(0,xx)
ax4.set_ylim(0,yy)
ax4.set_yticks(ytick_pos)
ax4.set_ylabel('')
ax4.set_yticklabels([])
ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
ax4.yaxis.set_minor_locator(AutoMinorLocator(5))
ax4.tick_params(which='minor', length=2)
ax4.set_xticks(xtick_pos)
ax4.set_xticklabels(xtick_lab, fontdict = font)
ax4.set_xlabel(r'x [arcsec]', fontdict=font)
ax4.text(10, 10 , text[4], fontdict = font, color = textcolor[4])

ax5 = plt.subplot(gs[2:,2],adjustable = 'box')
ax5.imshow(ca8_core[y1:y2,x1:x2], origin = 'lower', cmap = 'gray')
ax5.set_xticks(xtick_pos)
ax5.set_xlim(0,xx)
ax5.set_ylim(0,yy)
ax5.set_yticks(ytick_pos)
ax5.set_ylabel('')
ax5.set_yticklabels([])
ax5.xaxis.set_minor_locator(AutoMinorLocator(5))
ax5.yaxis.set_minor_locator(AutoMinorLocator(5))
ax5.tick_params(which='minor', length=2)
ax5.set_xticks(xtick_pos)
ax5.set_xticklabels(xtick_lab, fontdict = font)
ax5.set_xlabel(r'x [arcsec]', fontdict=font)
ax5.text(10, 10 , text[5], fontdict = font, color = textcolor[5])

ax6 = plt.subplot(gs[2:,3],adjustable = 'box')
ax6.imshow(cak_core[y1:y2,x1:x2], origin = 'lower', cmap = 'gray')

for ii in range(len(file_search(fibdir,'crispex*3950*.csav'))):
    #print("processing cut "+ str(ii) + " out of " + str(len(file_search(fibdir,'crispex*3950*.csav')))+"...")

    cut_file = (file_search(fibdir,'crispex*3950*.csav'))[ii]
    print(cut_file)
    cut = restore(fibdir+cut_file)
    xcut = cut.x_coords - x1
    ycut = cut.y_coords - y1
    
    if (cut_file == "crispex_3950_2018-07-22T08:23:58_scans=0-212_time-corrected_2021Oct04_165713.csav"):
        #plt.arrow(np.mean(xcut)+25, np.mean(ycut)+20, -20, 0, width =3., color = "red" )
        ax6.plot(xcut,ycut, color = "white", alpha = 0.5, linewidth = 1, linestyle = '--')
        ax6.plot(xcut[0], ycut[0], 'go', color = 'white', markersize = 2.5, alpha = 0.75)
    else:
        ax6.plot(xcut,ycut, color = "red", alpha = 0.5, linewidth = 1)
    plt.show()
    #stop()
    
ax6.set_xticks(xtick_pos)
ax6.set_xlim(0,xx)
ax6.set_ylim(0,yy)
ax6.set_yticks(ytick_pos)
ax6.set_ylabel('')
ax6.set_yticklabels([])
ax6.xaxis.set_minor_locator(AutoMinorLocator(5))
ax6.yaxis.set_minor_locator(AutoMinorLocator(5))
ax6.tick_params(which='minor', length=2)
ax6.set_xticks(xtick_pos)
ax6.set_xticklabels(xtick_lab, fontdict = font)
ax6.set_xlabel(r'x [arcsec]', fontdict=font)
ax6.text(10, 10 , text[6], fontdict = font, color = textcolor[6])

# clipping the image
#sigma_factor = 4.
#cak_min = np.mean(cak_scan) - 2.5*np.std(cak_scan)
#cak_max = np.mean(cak_scan) + 6*np.std(cak_scan)


#plt.tight_layout()
plt.show()
filename = 'maps_cuts.pdf'
plt.savefig(outdir + filename, dpi = 1000)

