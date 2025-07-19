import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
import batman
from exotoolbox.utils import reverse_ld_coeffs, tdur
import juliet
import plotstyles

# ------------------------------------------------------------------------------
#
#                       To plot phase curve models
#
# ------------------------------------------------------------------------------

# Planetary parameters (Patel & Espinoza 2022)
per, tc = 0.8134749, 2458555.80567
rprs, bb = 0.1615, 0.698
aR = 4.72
inc1 = np.rad2deg(np.arccos(bb / aR))
q1, q2 = 0.3106615413, 0.5494710423
u1, u2 = reverse_ld_coeffs('quadratic', q1, q2)

params = batman.TransitParams()
params.t0 = tc                       #time of inferior conjunction
params.per = per                     #orbital period
params.rp = rprs                     #planet radius (in units of stellar radii)
params.a = aR                        #semi-major axis (in units of stellar radii)
params.inc = inc1                    #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [u1, u2]                  #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

tecl = tc + (per/2)

t14 = tdur(per=per, ar=aR, rprs=rprs, bb=bb)

phases = np.linspace(-0.6, 0.6, 10000)
tim = tc + (phases*per)
#tim = np.linspace(tc-per/2-t14, tc+per/2+t14, 10000)
#phases = juliet.utils.get_phases(t=tim, P=per, t0=tc, phmin=1.)

fp, c1, d1 = 0.0057608005, 0.0021444080, -0.0004145816
c2, d2 = 0.0001746626, -0.0000055102

def CowanPC_model(times, te, per, E, C1, D1, C2, D2):
    omega_t = 2 * np.pi * (times - te) / per
    pc = E + ( C1 * (np.cos( omega_t ) - 1.) ) + ( D1 * np.sin( omega_t ) ) +\
             ( C2 * (np.cos( 2*omega_t ) - 1.) ) + ( D2 * np.sin( 2*omega_t ) )
    return pc

m = batman.TransitModel(params, tim)    #initializes model
flux_tra = m.light_curve(params)          #calculates light curve

params.fp = 80e-6
params.t_secondary = tc + (per/2)
m = batman.TransitModel(params, tim, transittype='secondary')    #initializes model
flux_ecl = m.light_curve(params)          #calculates light curve
flux_ecl = (flux_ecl - 1.) / 80e-6

pc = CowanPC_model(times=tim, te=tecl, per=per, E=fp, C1=c1, D1=d1, C2=c2, D2=d2)

fl_tot = flux_tra + (pc*flux_ecl)


fig, axs = plt.subplots(figsize=(16/3, 9/3))
axs.plot(phases, fl_tot, 'k-', zorder=10)
axs.axhline(1., color='k', ls='--', lw=0.7, zorder=5)
axs.set_xlim([-0.6, 0.6])

axs.set_xlabel('Orbital phase', labelpad=10)
axs.set_ylabel('Relative flux')

axs.xaxis.tick_top()
axs.xaxis.set_label_position('top')

plt.tight_layout()

#plt.show()
plt.savefig('Ch2/Fig28/full_pc.pdf')

fig, axs = plt.subplots(figsize=(16/3, 9/3))
axs.plot(phases, fl_tot, 'k-', zorder=10)
axs.axhline(1., color='k', ls='--', lw=0.7, zorder=5)
axs.set_xlim([-0.6, 0.6])
axs.set_ylim([0.997, 1.007])

axs.set_xlabel('Orbital phase', labelpad=10)
axs.set_ylabel('Relative flux')

plt.tight_layout()

#plt.show()
plt.savefig('Ch2/Fig28/full_pc_zoom_in.pdf')