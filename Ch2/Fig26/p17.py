import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as con
from matplotlib import cm
import matplotlib
import astropy.units as u
import plotstyles
import os

# ------------------------------------------------------------------------------
#
#   To plot emmision spectra of planet, star, and occultation depth spectra 
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

#star, planet = '55Cnc', '55Cnce'
star, planet = 'TOI-561', 'TOI-561b'
#star, planet = 'TOI-1807', 'TOI-1807b'

# Occultation depth spectrum
wav_fpfs, fpfs = np.loadtxt(os.getcwd() + '/Ch2/Fig26/Zilinskas/' + planet + '_FpFs.dat', usecols=(0,1), unpack=True)

# Planetary spectrum
wav_fp, fp = np.loadtxt(os.getcwd() + '/Ch2/Fig26/Zilinskas/' + planet + '_Fp.dat', usecols=(0,1), unpack=True)
fp = ( fp * 1e-6 * u.erg / u.s / u.cm**2 / u.Hz ).to(u.W / u.m**2 / u.Hz)

# Stellar spectrum
wav_st, st = np.loadtxt(os.getcwd() + '/Ch2/Fig26/Zilinskas/' + star + '_stel_spec.dat', usecols=(0,1), unpack=True)
wav_st = ( wav_st * u.cm ).to(u.micron)
st = ( st * u.erg / u.s / u.cm**3 ).to(u.W / u.m**2 / u.Hz, equivalencies=u.spectral_density(wav=wav_st))

# First let's plot planetary spectrum (Fp/F*)
fig, axs = plt.subplots(figsize=(16/2, 9/2))
axs.plot(wav_fpfs, fpfs, color=chex[2], lw=1.5)
axs.set_xlim([np.min(wav_fp), np.max(wav_fp)])

axs.set_xscale('log')

axs.set_xticks(ticks=np.array([1., 10.]),\
               labels=np.array(['1', '10']))

axs.set_xlabel('Wavelength [$\mu$m]')
axs.set_ylabel(r'F$_p$/F$_\star$ [ppm]')

plt.tight_layout()
#plt.show()
plt.savefig(os.getcwd() + '/Ch2/Fig26/pl_fpfs.pdf')#,dpi=500)"""

# Stellar spectrum
fig, axs = plt.subplots(figsize=(16/2.5, 9/2.5))
axs.plot(wav_st, st, color='crimson', lw=1.5)
axs.set_xlim([np.min(wav_fp), np.max(wav_fp)])

axs.set_xscale('log')

axs.set_xticks(ticks=np.array([1., 10.]),\
               labels=np.array(['1', '10']))

axs.set_xlabel(r'Wavelength [$\mu$m]')
axs.set_ylabel(r'Stellar Flux [W m$^{-2}$ Hz$^{-1}$]')

plt.tight_layout()

#plt.show()
plt.savefig(os.getcwd() + '/Ch2/Fig26/st_fl.pdf')#,dpi=500)

# Planet's spectrum
fig, axs = plt.subplots(figsize=(16/2.5, 9/2.5))

axs.plot(wav_fp, fp, color=chex[7], lw=1.5)

# Inset axes 1
axins = axs.inset_axes([0.7, 0.35, 0.25, 0.55], xlim=(6,11), ylim=(2e-9,5.2e-9), xticklabels=[], yticklabels=[])#, visible=True)
axins.plot(wav_fp, fp, color=chex[7], lw=1.5)

#axins.set_xticks(ticks=np.array([6., 7., 8., 9., 10.]),\
#               labels=np.array(['', '7', '', '', '10']))

axins1 = axs.inset_axes([0.27, 0.1, 0.15, 0.5], xlim=(0.55,0.65), ylim=(0.25e-8,1.1e-8), xticklabels=[], yticklabels=[])#, visible=True)
axins1.plot(wav_fp, fp, color=chex[7], lw=1.5)

axs.set_xscale('log')
axs.set_xticks(ticks=np.array([1., 10.]),\
               labels=np.array(['1', '10']))

axs.set_xlabel(r'Wavelength [$\mu$m]')
axs.set_ylabel(r'Planetary Flux [W m$^{-2}$ Hz$^{-1}$]')

axs.set_xlim([np.min(wav_fp), np.max(wav_fp)])
axs.indicate_inset_zoom(axins, edgecolor="black")
axs.indicate_inset_zoom(axins1, edgecolor="black")

plt.tight_layout()
#plt.show()
plt.savefig(os.getcwd() + '/Ch2/Fig26/pl_fl.pdf')#,dpi=500)