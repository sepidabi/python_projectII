# fibril_hor_oscillation.py

'''Meant for plotting a sub field of view
with a fibril and its horizontal oscillations
(the cuts along the fibrils are already defined via crispex)
'''

# import modules
import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# Font properties
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }


# Directories
datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
outdir = '/home/seki2695/OUTPUT/project2/'

# Difinitions
res = 0.0375 # CHROMIS pixel size in arcse
cad = 8

# Reading the observed data
pref = ['6173','8542','3950']
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 =file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak =file_search(datadir,'crispex*'+pref[2]+'*.fcube')

cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])

# getting the fibril and cuts information
fibdir = datadir+'fr29/'
fib_file = (file_search(fibdir,'crispex*3950*.csav'))[0]
fib = restore(fibdir+fib_file)
cut1_file = (file_search(datadir,'crispex*3950*.csav'))[1]
cut2_file = (file_search(datadir,'crispex*3950*.csav'))[2]
#cut3_file = (file_search(datadir,'crispex*3950*.csav'))[3]
cut4_file = (file_search(datadir,'crispex*3950*.csav'))[4]

cut1 = restore(datadir+cut1_file)
cut2 = restore(datadir+cut2_file)
#cut3 = restore(datadir+cut3_file)
cut4 = restore(datadir+cut4_file)

print(np.sqrt((cut1.x_coords[0]-cut1.x_coords[1])**2 + (cut1.y_coords[0]-cut1.y_coords[1])**2),
      cut1.loop_size)
print(np.sqrt((cut2.x_coords[0]-cut2.x_coords[1])**2 + (cut2.y_coords[0]-cut2.y_coords[1])**2),
      cut2.loop_size)
print(np.sqrt((cut4.x_coords[0]-cut4.x_coords[1])**2 + (cut4.y_coords[0]-cut4.y_coords[1])**2),
      cut4.loop_size)

xmin = int(np.min([np.min([cut1.x_coords, cut2.x_coords,
                           #cut3.x_coords,
                           cut4.x_coords]),np.min(fib.x_coords)]) - 5)
ymin = int(np.min([np.min([cut1.y_coords, cut2.y_coords,
                           #cut3.y_coords,
                           cut4.y_coords]),np.min(fib.y_coords)]) - 5)
xmax = int(np.max([np.max([cut1.x_coords, cut2.x_coords,
                           #cut3.x_coords,
                           cut4.x_coords]),np.max(fib.x_coords)]) + 5)
ymax = int(np.max([np.max([cut1.y_coords, cut2.y_coords,
                           #cut3.y_coords,
                           cut4.y_coords]),np.max(fib.y_coords)]) + 5)

# PLOT the left pannel
xfig_left = 3.5
xfig_right = 4
yfig = 5.
plt.close('all')
f = plt.figure(figsize = (xfig_left, yfig))

gs1 = gridspec.GridSpec(2, 2) # grid scale of the 1st col.
gs1.update(left=0.15,
           right=0.98,
           wspace=0.,
           bottom=0.1,
           top=0.98,
           hspace = 0.
)

# axis info
rr = 40/1.5
ytick_pos = np.arange(0,(np.round((ymax-ymin)/rr)+1)*rr,rr)
ytick_lab = np.intc(np.round(ytick_pos*res, decimals = 0))
xtick_pos = np.arange(0,(np.round((xmax-xmin)/rr)+1)*rr,rr)
xtick_lab = np.intc(np.round(xtick_pos*res, decimals = 0))
xtick_lab1 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut1.x_coords[0]-cut1.x_coords[1])**2 +
                                               (cut1.y_coords[0]-cut1.y_coords[1])**2)/cut1.loop_size),
                       decimals = 0))
xtick_lab2 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut2.x_coords[0]-cut2.x_coords[1])**2 +
                                               (cut2.y_coords[0]-cut2.y_coords[1])**2)/cut2.loop_size),
                       decimals = 0))
xtick_lab4 = np.intc(np.round(xtick_pos*res*(np.sqrt((cut4.x_coords[0]-cut4.x_coords[1])**2 +
                                               (cut4.y_coords[0]-cut4.y_coords[1])**2)/cut4.loop_size),
                       decimals = 0))

tt = 30
ttick_pos = np.arange(0,(np.round((213)/tt)+1)*tt,tt)
ttick_lab = ttick_pos*cad/60

ax11 = plt.subplot(gs1[0,0],adjustable = 'box')#,aspect = 'equal')
#ax11.set_xlabel(r'x [arcsec]', fontdict = font)
ax11.set_ylabel(r'y [arcsec]', fontdict = font)
ax11.xaxis.set_minor_locator(AutoMinorLocator(5))
ax11.yaxis.set_minor_locator(AutoMinorLocator(5))
ax11.tick_params(which='minor', length=2)
ax11.set_xticks(xtick_pos)
ax11.set_xticklabels([])
ax11.set_yticks(ytick_pos)
ax11.set_yticklabels(ytick_lab, fontdict = font)

ax12 = plt.subplot(gs1[1,0],adjustable = 'box')#,aspect = 'equal')
ax12.set_xlabel(r'x [arcsec]', fontdict = font)
ax12.set_ylabel(r'y [arcsec]', fontdict = font)
ax12.xaxis.set_minor_locator(AutoMinorLocator(5))
ax12.yaxis.set_minor_locator(AutoMinorLocator(5))
ax12.tick_params(which='minor', length=2)
ax12.set_xticks(xtick_pos)
ax12.set_xticklabels(xtick_lab, fontdict = font)
ax12.set_yticks(ytick_pos)
ax12.set_yticklabels(ytick_lab, fontdict = font)

ax13 = plt.subplot(gs1[0,1],adjustable = 'box')#,aspect = 'equal')
#ax13.set_xlabel(r'x [arcsec]', fontdict = font)
#ax13.set_ylabel(r'y [arcsec]', fontdict = font)
ax13.xaxis.set_minor_locator(AutoMinorLocator(5))
ax13.yaxis.set_minor_locator(AutoMinorLocator(5))
ax13.tick_params(which='minor', length=2)
ax13.set_xticks(xtick_pos)
ax13.set_xticklabels([])
ax13.set_yticks(ytick_pos)
ax13.set_yticklabels([])

ax14 = plt.subplot(gs1[1,1],adjustable = 'box')#,aspect = 'equal')
ax14.set_xlabel(r'x [arcsec]', fontdict = font)
#ax14.set_ylabel(r'y [arcsec]', fontdict = font)
ax14.xaxis.set_minor_locator(AutoMinorLocator(5))
ax14.yaxis.set_minor_locator(AutoMinorLocator(5))
ax14.tick_params(which='minor', length=2)
ax14.set_xticks(xtick_pos)
ax14.set_xticklabels(xtick_lab, fontdict = font)
ax14.set_yticks(ytick_pos)
ax14.set_yticklabels([])

fr_select = [29,75,135,167]

ax11.imshow(cube_cak[fr_select[0],0,13,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower')
#ax11.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
#ax11.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
#ax11.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
#plt.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
#ax11.plot(cut4.x_coords-xmin, cut4.y_coords-ymin, color = 'red')

ax12.imshow(cube_cak[fr_select[1],0,13,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower')
#ax12.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
#ax12.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
#ax12.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
#plt.plot(cut3.x_coords-xmin, cut3.y_coords-ymin, color = 'red')
#ax12.plot(cut4.x_coords-xmin, cut4.y_coords-ymin, color = 'red')

ax13.imshow(cube_cak[fr_select[2],0,13,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower')
#ax13.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
#ax13.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
#ax13.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')

ax14.imshow(cube_cak[fr_select[3],0,13,ymin:ymax,xmin:xmax], cmap = 'gray', origin = 'lower')
ax14.plot(fib.x_coords-xmin, fib.y_coords-ymin, color = 'red')
ax14.plot(cut1.x_coords-xmin, cut1.y_coords-ymin, color = 'red')
ax14.plot(cut2.x_coords-xmin, cut2.y_coords-ymin, color = 'red')
ax14.plot(cut4.x_coords-xmin, cut4.y_coords-ymin, color = 'red')
ax14.text(cut1.x_coords[0]-xmin-4, cut1.y_coords[0]-ymin-4, '(3)', fontdict=font, color = 'red')
ax14.text(cut2.x_coords[0]-xmin-4, cut2.y_coords[0]-ymin-4, '(2)', fontdict=font, color = 'red')
ax14.text(cut4.x_coords[0]-xmin-4, cut4.y_coords[0]-ymin-4, '(1)', fontdict=font, color = 'red')

plt.show()
filename = outdir+'oscilation_left.pdf'
f.savefig(filename, quality = 100)
#f.tight_layout()
print 'file saved to: ' + filename

#stop()

# PLOT right panel
plt.close('all')
f = plt.figure(figsize = (xfig_right,yfig))

gs = gridspec.GridSpec(1, 3) # grid scale of the 1st col.
gs1.update(left=0.15,
           right=1.,
           wspace=0.,
           bottom=0.08,
           top=0.98,
           hspace = 0.
)


ax2 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
ax2.set_xlabel(r'length [arcsec]', fontdict = font)
ax2.set_yticks(ttick_pos)
ax2.set_yticklabels(ttick_lab, fontdict = font)
ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
ax2.set_ylabel(r'$t$ [m]', fontdict = font)
ax2.set_xticks(xtick_pos)
ax2.set_xticklabels(xtick_lab4, fontdict = font)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.set_ylim(0,cube_cak.shape[0]-1)

ax3 = plt.subplot(gs[0,1],adjustable = 'box')#,aspect = 'equal')
ax3.set_yticks([])
ax3.set_xlabel(r'length [arcsec]', fontdict = font)
ax3.set_yticks(ttick_pos)
ax3.set_yticklabels([])
ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
ax3.set_xticklabels(xtick_lab2, fontdict = font)
ax3.set_xticks(xtick_pos)
ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
ax3.set_ylim(0,cube_cak.shape[0]-1)

ax4 = plt.subplot(gs[0,2],adjustable = 'box')#,aspect = 'equal')
ax4.set_yticks([])
ax4.set_xlabel(r'length [arcsec]', fontdict = font)
ax4.set_yticks(ttick_pos)
ax4.set_yticklabels([])
ax4.yaxis.set_minor_locator(AutoMinorLocator(4))
ax4.set_xticklabels(xtick_lab1, fontdict = font)
ax4.set_xticks(xtick_pos)
ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
ax4.set_ylim(0,cube_cak.shape[0]-1)

ratio =0.8
ax2.imshow(cut4.loop_slab[14,:,:], cmap = 'gray', origin = 'lower', aspect=ratio)
ax2.text(2,205,'(1)', color = 'red')

ax3.imshow(cut2.loop_slab[14,:,:], cmap = 'gray', origin = 'lower', aspect=ratio)
ax3.text(2,205,'(2)', color = 'red')

ax4.imshow(cut1.loop_slab[14,:,:], cmap = 'gray', origin = 'lower', aspect=ratio)
ax4.text(2,205,'(3)', color = 'red')

plt.show()

#print('YOU NEED TO ADJUST THE MARGINS OF THE FIGURE MANUALLY')
#stop()
filename = outdir+'oscilation_right.pdf'
f.savefig(filename, quality = 100)
f.tight_layout()
print 'file saved to: ' + filename
