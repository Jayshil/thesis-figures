import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import ConnectionPatch
from matplotlib import cm
import matplotlib
import os
import plotstyles

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
for i in range(11):
    c1 = cm.plasma(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10

# Sector 8, time 500
fname = os.getcwd() + '/Data/TESS/tess2019036055935-s0008-1-3-0136-s/tess2019036055935-s0008-1-3-0136-s_ffic.fits'
# Sector 35, time 500
#fname = os.getcwd() + '/mastDownload/TESS/tess2021043234908-s0035-1-4-0205-s/tess2021043234908-s0035-1-4-0205-s_ffic.fits'
# Sector 45, time 500
#fname = os.getcwd() + '/mastDownload/TESS/tess2021312154857-s0045-3-3-0216-s/tess2021312154857-s0045-3-3-0216-s_ffic.fits'

hdul = fits.open(fname)

# WCS and data
wcs = WCS(hdul[1].header)
cal = hdul[1].data
hdr = hdul[1].header
mid_time = hdr['BJDREFI'] + (hdr['TSTOP'] + hdr['TSTART'])/2

print(wcs)

print(wcs.pixel_to_world_values(1000, 1000))
print(wcs.world_to_pixel_values(148.18516666666667, 6.2161027777777775))

pix_x, pix_y = wcs.world_to_pixel_values(148.18516666666667, 6.2161027777777775)
print(pix_x, pix_y)

print(mid_time)

# Plot
plt.figure(figsize = (6.,6.))
axs = plt.subplot(111, projection = wcs)
axs.imshow(cal, vmin = np.percentile(cal,4), vmax = np.percentile(cal, 98), origin = "lower", cmap='plasma')
axs.set_xlabel('RA', fontsize=16, fontfamily='serif')
axs.set_ylabel('Dec', fontsize=16, fontfamily='serif')
axs.tick_params(labelfontfamily='serif')
axs.set_title("TESS FFI for Sector 8, Camera 1, CCD 3\n at %f BJD" % mid_time, fontfamily='serif')
axs.grid()

plt.setp(axs.get_xticklabels(), fontsize=15)
plt.setp(axs.get_yticklabels(), fontsize=15)

## Inset-axis for zoom-in on TESS bandpass
axins = axs.inset_axes([0.05, 0.05, 0.30, 0.30], visible=True)
axins.imshow(cal, vmin = np.percentile(cal,4), vmax = np.percentile(cal, 98), origin = "lower", cmap='plasma')
axins.set_xlim(pix_x-50, pix_x+100)
axins.set_ylim(pix_y-50, pix_y+100)

clr = chex[0]
axins.vlines(x=pix_x, ymin=pix_y-20, ymax=pix_y-10, lw=2.5, colors=clr)
axins.vlines(x=pix_x, ymin=pix_y+10, ymax=pix_y+20, lw=2.5, colors=clr)
axins.hlines(y=pix_y, xmin=pix_x-20, xmax=pix_x-10, lw=2.5, colors=clr)
axins.hlines(y=pix_y, xmin=pix_x+10, xmax=pix_x+20, lw=2.5, colors=clr)

axins.set_xticklabels([])
axins.set_yticklabels([])

#axins.scatter(pix_x, pix_y, color='r', marker='o')

axs.indicate_inset_zoom(axins, edgecolor="black", lw=2., alpha=1.)

cp1 = ConnectionPatch(xyA=(pix_x-50, pix_y-50), xyB=(0, 0), axesA=axs, axesB=axins,
                      coordsA="data", coordsB="axes fraction", lw=1.9)#, ls=":")
cp2 = ConnectionPatch(xyA=(pix_x+100, pix_y+100), xyB=(1, 1), axesA=axs, axesB=axins,
                      coordsA="data", coordsB="axes fraction", lw=1.9)#, ls=":")

plt.tight_layout()
#plt.show()
plt.savefig('tess_ffi.pdf')#, dpi=500)