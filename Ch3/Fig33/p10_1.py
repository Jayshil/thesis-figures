import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import cm
from matplotlib.patches import Wedge
import matplotlib
import os

# ------------------------------------------------------------------------------
#
#                   To plot CHEOPS Subarray and Imagette
#
# ------------------------------------------------------------------------------

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
for i in range(11):
    c1 = cm.plasma(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10


# ------------------------------
# 
#           Subarray
#
# ------------------------------

hdul = fits.open(os.getcwd() + '/Ch3/Fig33/CH_PR110047_TG000401_TU2021-07-11T04-05-52_SCI_CAL_SubArray_V0300.fits')
cal = hdul[1].data

iter = 51
cal_iter = cal[iter, :, :]
up_lim, lo_lim = np.nanmedian(cal_iter) + 4000*np.nanmedian(np.abs(np.diff(cal_iter))), np.nanmedian(cal_iter) - 4000*np.nanmedian(np.abs(np.diff(cal_iter)))

fig, axs = plt.subplots(figsize=(5.,5.))
axs.imshow(cal_iter, origin = "lower", cmap='plasma', vmin = lo_lim, vmax = up_lim)

plt.xticks(ticks=np.array([]), labels=np.array([]))
plt.yticks(ticks=np.array([]), labels=np.array([]))
plt.axis('off')
plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig33/kelt20_subarray.png', dpi=250)

# ------------------------------
# 
#           Imagette
#
# ------------------------------

hdul = fits.open('Ch3/Fig33/CH_PR110047_TG000401_TU2021-07-11T04-05-40_SCI_RAW_Imagette_V0300.fits')
cal = hdul[1].data

iter = 51
cal_iter = cal[iter, :, :]
up_lim, lo_lim = np.nanmedian(cal_iter) + 350*np.nanmedian(np.abs(np.diff(cal_iter))), np.nanmedian(cal_iter) - 350*np.nanmedian(np.abs(np.diff(cal_iter)))

fig, axs = plt.subplots(figsize=(5.,5.))
axs.imshow(cal_iter, origin = "lower", cmap='plasma', vmin = lo_lim, vmax = up_lim)

plt.xticks(ticks=np.array([]), labels=np.array([]))
plt.yticks(ticks=np.array([]), labels=np.array([]))
plt.axis('off')
plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig33/kelt20_imagette.png', dpi=250)