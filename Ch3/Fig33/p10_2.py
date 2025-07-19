import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import os
import plotstyles

# ------------------------------------------------------------------------------
#
#                   To plot CHEOPS DRP light curve
#
# ------------------------------------------------------------------------------

hdul = fits.open('Ch3/Fig33/CH_PR110047_TG000401_TU2021-07-11T04-05-52_SCI_COR_Lightcurve-DEFAULT_V0300.fits')

tab = Table(hdul[1].data)

tim = (tab['BJD_TIME'] - tab['BJD_TIME'][0])*24
fl, fle = tab['FLUX'] / np.nanmedian(tab['FLUX']), tab['FLUXERR'] / np.nanmedian(tab['FLUX'])

fig, ax = plt.subplots(figsize=(16/2, 9/2))
ax.errorbar(tim, fl, yerr=fle, fmt='.', c='k')

ax.set_ylim([0.984, 1.0015])
ax.set_xlim([tim[0], tim[-1]])

ax.set_yticks(ticks=np.array([1.00, 0.995, 0.990, 0.985]),\
              labels=np.array(['1.000', '0.995', '0.990', '0.985']))
              #labels=np.array([1.0, 0.998, 0.996, 0.994, 0.992, 0.990, 0.988, 0.986, 0.984]))

ax.set_xlabel('Time since the beginning (hr)', fontsize=24)
ax.set_ylabel('Relative flux', fontsize=24)
plt.setp(ax.get_xticklabels(), fontsize=20)
plt.setp(ax.get_yticklabels(), fontsize=20)

plt.tight_layout()

#plt.show()
plt.savefig('Ch3/Fig33/kelt20_drplc.pdf')#, dpi=250)