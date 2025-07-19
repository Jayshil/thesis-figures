import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib
import os
import plotstyles

# ------------------------------------------------------------------------------
#
#                       To plot CHEOPS subarray
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

hdul = fits.open(os.getcwd() + '/Ch3/Fig31/CH_PR100008_TG000801_TU2021-02-28T20-16-28_SCI_CAL_SubArray_V0200.fits')

cal = hdul[1].data
tim = hdul[2].data['BJD_TIME']

iter = 51

cal_iter = cal[iter, :, :]

up_lim, lo_lim = np.nanmedian(cal_iter) + 1000*np.nanmedian(np.abs(np.diff(cal_iter))), np.nanmedian(cal_iter) - 1000*np.nanmedian(np.abs(np.diff(cal_iter)))

fig, axs = plt.subplots(figsize=(6.,6.))
axs.imshow(cal_iter, origin = "lower", cmap='plasma', vmin = lo_lim, vmax = up_lim)
axs.set_xlabel('Column number', fontsize=16)
axs.set_ylabel('Row number', fontsize=16)
axs.set_title("CHEOPS calibrated subarray data at\n %f BJD" % tim[iter])

plt.setp(axs.get_xticklabels(), fontsize=14)
plt.setp(axs.get_yticklabels(), fontsize=14)

plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig31/cheops_subarray.pdf')#, dpi=500)