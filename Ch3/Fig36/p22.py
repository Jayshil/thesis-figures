import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from glob import glob
import plotstyles

# ------------------------------------------------------------------------------
#
#                          To plot MIRI 2D spectrum
#
# ------------------------------------------------------------------------------

# --------------- Loading the data
## 2D Specrum data
corrected_2d_spec_ints = np.load('Ch3/Fig36/Corrected_data.npy')
corrected_2d_spec = np.nanmedian(corrected_2d_spec_ints, axis=0)

## Trace positions and wavelengths
#tr1 = pickle.load(open(p1 + '/stark/Outputs/StelSpec/Traces_seg007_V4.pkl', 'rb'))
xpos = np.arange(corrected_2d_spec.shape[0])
jw = fits.open(glob('Ch3/Fig36/*mirimage_calints.fits')[0])
wav_soln = jw['WAVELENGTH'].data
wav_obs = np.zeros(len(xpos))
for i in range(len(xpos)):
    wav_obs[i] = wav_soln[xpos[i], 35]

print(corrected_2d_spec.shape, wav_obs.shape)

#plt.plot(corrected_2d_spec[:,1500], 'k-')
#plt.show()

corrected_2d_spec = np.flip(corrected_2d_spec.T, axis=1)
wav_obs = np.flip(wav_obs)

# ---------------- Plotting data
# construct cmap
my_cmap = 'plasma'#sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)

fig, axs = plt.subplots(figsize=(15/1.5,5/1.5))
im = axs.imshow(corrected_2d_spec, aspect='auto', cmap=my_cmap, interpolation='none')#, cmap='plasma')#, interpolation='none')
im.set_clim([0.4e2, 5e2])

ticks = np.arange(21, 408, 48)
ticklabels = ["{:.2f}".format(wav_obs[i]) for i in ticks]
axs.set_xticks(ticks)
axs.set_xticklabels(ticklabels, fontsize=14)

#axs.xaxis.tick_top()
#axs.xaxis.set_label_position('top')

axs.set_xlabel('Wavelength [$\mu$m]', fontsize = 15, labelpad=10)
axs.set_ylabel('Spatial direction', fontsize = 15)

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
plt.savefig('Ch3/Fig36/wasp43_2d_miri_spec.pdf')