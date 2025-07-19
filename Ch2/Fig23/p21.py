import numpy as np
import matplotlib.pyplot as plt
from pytransit import OblateStarModel
from pytransit.contamination import TabulatedFilter
import astropy.constants as con
import astropy.units as u
import plotstyles
import batman
import juliet
import os

# ------------------------------------------------------------------------------
#
#   To plot transit light curves with batman and pytransit (gravity darkening) 
#
# ------------------------------------------------------------------------------

def symmetric_transit(tim, per, t0, rprs, inc, ar, u1, u2, ecc, w):
    pars = batman.TransitParams()
    pars.per = per
    pars.t0 = t0
    pars.rp = rprs
    pars.inc = inc
    pars.a = ar
    pars.ecc = ecc
    pars.w = w
    pars.limb_dark = 'quadratic'
    pars.u = [u1, u2]

    m = batman.TransitModel(params=pars, t=tim)
    return m.light_curve(params=pars)

def gravity_darkened_transit(tim, filter, per, t0, rprs, inc, ar, u1, u2, ecc, w, mst, rst, tpole, phi, lamp, vsini):
    model = OblateStarModel(filters=filter, rstar=rst, model='husser2013', sres=100, pres=8, tmin=5000, tmax=10000)
    # Converting angles to radians
    phi, lamp, inc, w = np.radians(phi), np.radians(lamp), np.radians(inc), np.radians(w)
    # Computing rho and rper from the existing data
    ## rper
    rper = 2*np.pi*rst*con.R_sun/(vsini*(u.km/u.s))
    rper = rper * np.sin( (np.pi/2) - phi )
    rper = rper.to(u.day)
    rper = rper.value
    ## rho
    term1 = con.G * mst*con.M_sun * rper*u.day * rper*u.day/(2 * np.pi * np.pi * ((rst * con.R_sun)**3))
    term1 = term1.decompose()
    fstar = 1/(1 + term1)
    rho = 3 * np.pi/(2 * con.G * fstar * rper*u.day * rper*u.day)
    rho = rho.to(u.g/u.cm**3)
    rho = rho.value
    ## LDCs
    ldcs = np.array([u1, u2])
    ## rprs to numpy.ndarray (as `pytransit` only takes array for k)
    rprs = np.array([rprs])

    # Now computing models
    ## Oblate star transit model
    model.rstar = rst*con.R_sun.value
    model.set_data(tim)
    flux_tra = model.evaluate_ps(k=rprs, rho=rho, rperiod=rper, tpole=tpole, phi=phi, beta=0.199, ldc=ldcs,\
        t0=t0, p=per, a=ar, i=inc, l=lamp, e=ecc, w=w)
    
    return flux_tra

# Planetary parameters
per, rprs = 2.14877381, 0.07884
inc, rho = 88.45, 0.2966*1e3
q1, q2 = 0.234, 0.405
ecc, w = 0.00034, -16
rst, mst = 2.082, 1.900
vsini, tpole = 101.7, 7490
phi, lamp = 90-55.5, -69.2
ar = 4.1676
t14 = 4.226/24
u1, u2 = juliet.utils.reverse_ld_coeffs('quadratic', q1, q2)

wav_c, band_c = np.loadtxt(os.getcwd() + '/Data/cheops_response_fun.txt', usecols=(0,1), unpack=True)
wav_c, band_c = wav_c/10, band_c/np.max(band_c)
### Saving bandpass so that `pytransit` can understand it!
filt = TabulatedFilter('CHEOPS', wav_c, band_c)

times = np.linspace(-1.2*t14, 1.2*t14, 10000)

batman_tra = symmetric_transit(tim=times, per=per, t0=0, rprs=rprs, inc=inc, ar=ar, u1=u1, u2=u2, ecc=ecc, w=w)
pytransit_tra = gravity_darkened_transit(tim=times, filter=filt, per=per, t0=0, rprs=rprs, inc=inc, ar=ar,\
                                         u1=u1, u2=u2, ecc=ecc, w=w, mst=mst, rst=rst, tpole=tpole, phi=phi,\
                                         lamp=lamp, vsini=vsini)

fig, axs = plt.subplots(figsize=(15/2.5, 9/2.5))
axs.plot(times/per, pytransit_tra, color='crimson', zorder=10)
axs.plot(times/per, batman_tra, color='dodgerblue', lw=1, zorder=5)

axs.set_xlabel('Orbital phase')
axs.set_ylabel('Relative flux')

axs.text(-0.0379, 0.9995, 'Gravity darkened transit', color='crimson', fontsize=14)
axs.text(-0.02775, 0.9988, 'Symmetric transit', color='dodgerblue', fontsize=14)

axs.set_xlim([-0.075, 0.075])

axs.set_yticks(ticks=np.array([0.992, 0.994, 0.996, 0.998, 1.000]))
axs.set_xticks(ticks=np.array([-0.06, -0.03, 0., 0.03, 0.06]))

plt.tight_layout()
#plt.show()
plt.savefig('Ch2/Fig23/pytransit.pdf')