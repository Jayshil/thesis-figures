import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Wedge
import plotstyles
import os

# ------------------------------------------------------------------------------
#
#                         To plot CHEOPS full array
#
# ------------------------------------------------------------------------------

from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

fname = os.getcwd() + '/Ch3/Fig33/CH_PR110047_TG000401_TU2021-07-11T04-04-50_SCI_RAW_FullArray_V0300.fits'

hdul = fits.open(fname)
cal = hdul[1].data

# C: 250-310; R: 800-860

# Plot
plt.figure(figsize = (4.,4.))
axs = plt.subplot(111)#, projection = wcs)
axs.imshow(cal, vmin = np.percentile(cal,2), vmax = np.percentile(cal, 98), origin = "lower", cmap='plasma')

# Adding a sub-array aperture
sub_array = Wedge((280., 828.), 100., 0, 360, width=0.5, color='white')
axs.add_patch(sub_array)
imagette = Wedge((280., 828.), 30., 0, 360, width=0.5, color='cyan')#, ls='--')
axs.add_patch(imagette)

plt.tight_layout()
plt.show()
plt.savefig('Ch3/Fig33/kelt20_full_array.pdf')#, dpi=250)