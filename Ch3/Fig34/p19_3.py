import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib import cm
import seaborn as sns
import os
import plotstyles
import h5py

# ------------------------------------------------------------------------------
#
#                   To plot MIRI spectral light curves
#
# ------------------------------------------------------------------------------

p1 = os.getcwd() + '/Ch3/Fig34/1_Light_Curves'

# --------------- Loading the data
f1 = h5py.File(p1 + '/sparta.h5')

time = np.asarray( f1['time'] )
wave = np.asarray( f1['wavelength'] )
fl, fle = np.asarray( f1['flux'] ), np.asarray( f1['err'] )
msk = np.asarray( f1['mask_white'] )

mask = np.ones(len(msk), dtype=bool)
mask[msk == 1] = False

time, fl, fle = time[mask], fl[:,mask], fle[:,mask]
times_hours = (time + 2400000.5 - 2459915.1207842459) * 24
#times_hours = (time - np.mean(time)) * 24

med_fl = np.nanmedian(fl, axis=1)

fl = fl / med_fl[:,None]
fle = fle / med_fl[:,None]

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
cmap = sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
for i in range(11):
    c1 = cmap(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10

fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, sharex=True, figsize=(25/2.3, 7/2.3))

axs[0].errorbar(times_hours, fl[2,:], yerr=fle[2,:], fmt='.', color=chex[1])#, mew=0.3, zorder=10)#, mec='k')
axs[1].errorbar(times_hours, fl[6,:], yerr=fle[6,:], fmt='.', color=chex[3])
axs[2].errorbar(times_hours, fl[10,:], yerr=fle[10,:], fmt='.', color=chex[9])

axs[0].text(-10., 0.98, str(wave[2]) + r' $\mu$m', color=chex[1])
axs[1].text(-10., 0.98, str(wave[6]) + r' $\mu$m', color=chex[3])
axs[2].text(-10., 0.98, str(wave[10]) + r' $\mu$m', color=chex[9])

axs[2].set_xlim([np.min(times_hours), np.max(times_hours)])
axs[2].set_ylim([0.9675, 1.012])

axs[0].set_ylabel('Relative flux')

fig.supxlabel('Time since mid-transit [hr]', y=0.1)

sns.despine()
#axs[1].spines.left.set_visible(False)
#axs[2].spines.left.set_visible(False)
#axs[1].yaxis.set_visible(False)
#axs[2].yaxis.set_visible(False)

plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig34/w43_miri_lightcurves.pdf', transparent=True)

print(chex[1], chex[3], chex[9])