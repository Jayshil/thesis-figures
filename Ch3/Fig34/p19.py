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
#                   To plot MIRI 2D spectral timeseries
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
print(wave)

fl = fl / med_fl[:,None]
fle = fle / med_fl[:,None]

print(fl.T.shape)

fl = fl[0:12,:]
wave = wave[0:12]

print(fl.T.shape)

# ---------------- Plotting data
# construct cmap
my_cmap = 'plasma'#sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)

fig, axs = plt.subplots(figsize=(15/1.5,5/1.5))
im = axs.imshow(fl.T, aspect='auto', cmap=my_cmap)#, cmap='plasma')#, interpolation='none')
im.set_clim([0.98, 1.003])

#axs.axvline(2, color='k')

ticks = np.arange(0, len(wave), 2)
ticklabels = ["{:.2f}".format(wave[i]) for i in ticks]
axs.set_xticks(ticks)
axs.set_xticklabels(ticklabels, fontsize=14)

axs.xaxis.tick_top()
axs.xaxis.set_label_position('top')

# Y axis:
ticks = np.arange(0, len(time), 1685)
ticklabels = ["{:.1f}".format(times_hours[i]) for i in ticks]
axs.set_yticks(ticks)
axs.set_yticklabels(ticklabels, fontsize=14)

axs.set_xlabel('Wavelength [$\mu$m]', fontsize = 15, labelpad=10)
axs.set_ylabel('Time since mid transit [hr]', fontsize = 15)

# Colorbar:
divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="3%", pad=0.05)

cbar = fig.colorbar(im, shrink = 0.08, cax=cax)
cbar.ax.tick_params(labelsize=14)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel('Relative flux', rotation=270, fontsize=15, labelpad=20)

#axs.set_xlim([5,10.5])

plt.tight_layout()

#plt.show()
plt.savefig('Ch3/Fig34/wasp43_2d_miri_top_label.pdf')