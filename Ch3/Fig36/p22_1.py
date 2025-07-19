import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as con
from astropy.io import fits
from glob import glob
import pickle
import plotstyles

# ------------------------------------------------------------------------------
#
#                   To plot MIRI 1D spectrum
#
# ------------------------------------------------------------------------------

# Stellar properties
rst, dist = 0.939, 213.982
teff = 5485

# Computing solid angle to the star
solid_angle = (np.pi * ((rst * con.R_sun)**2)) / (((dist * u.pc).to(u.m))**2)
solid_angle_1au = (np.pi * ((rst * con.R_sun)**2)) / (((1. * u.au).to(u.m))**2)


## Trace positions and wavelengths
#tr1 = pickle.load(open(p1 + '/stark/Outputs/StelSpec/Traces_seg007_V4.pkl', 'rb'))
xpos = np.arange(150, 390, 1)
jw = fits.open(glob('Ch3/Fig36/*mirimage_calints.fits')[0])
wav_soln = jw['WAVELENGTH'].data
wav_obs = np.zeros(len(xpos))
for i in range(len(xpos)):
    wav_obs[i] = wav_soln[xpos[i], 35]

# Loading the dataset
spec_cube = pickle.load(open('Ch3/Fig36/Spectrum_cube_seg.pkl', 'rb'))
spectra, vars = spec_cube['spectra'], spec_cube['variance']
obs_spec = np.nanmedian(spectra, axis=0) * u.Jy / solid_angle# / u.sr
obs_spec = obs_spec.to(u.W / u.m**2 / u.micron, equivalencies=u.spectral_density(wav=wav_obs*u.micron))

# Planck blackbody energy density for the star
def planck_func_at_photosphere(lam, temp):
    """
    Given the wavelength and temperature
    this function will compute the specific
    intensity using the Planck's law
    """
    coeff1 = (2*con.h*con.c*con.c)/(lam**5)
    expo = np.exp((con.h*con.c)/(lam*con.k_B*temp)) - 1
    planck = (coeff1/expo).to(u.W / u.m**2 / u.micron)
    return planck
fl_planck_star = planck_func_at_photosphere(wav_obs*u.micron, teff*u.K)     # Blackbody energy density of W/m2/micron

fig, axs = plt.subplots(figsize=(15/1.5, 6/1.5))
axs.plot(wav_obs[:-2], obs_spec[:-2], color='mediumpurple')

axs.set_xlabel('Wavelength [$\mu$m]')
axs.set_ylabel('Intensity [W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]')

axs.set_yscale('log')

plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig36/stel_spec_w43_miri.pdf')