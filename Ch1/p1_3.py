import numpy as np
import matplotlib.pyplot as plt
import batman
from exotoolbox.utils import reverse_ld_coeffs
from juliet.utils import get_phases
import plotstyles

# ------------------------------------------------------------------------------
#
#                        To plot transit light curve
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

tim = np.linspace(tc-0.15, tc+0.15, 10000)

m = batman.TransitModel(params, tim)    #initializes model
flux = m.light_curve(params)          #calculates light curve

fig, ax = plt.subplots(figsize=(16/2, 9/2))
ax.plot((tim-tc)*24, flux, color='k')

ax.set_yticks(ticks=np.array([1.0, 0.995, 0.990, 0.985, 0.980, 0.975, 0.970]),\
              labels=np.array(['1.0', '', '0.99', '', '0.98', '', '0.97']))
ax.set_xticks(ticks=np.array([-1.0, -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.0]),\
              labels=np.array(['$-$1', '', '$-$0.5', '', '0', '', '0.5', '', '1']))

ax.set_xlabel('Time since mid transit [hr]', fontsize=24)
ax.set_ylabel('Relative flux', fontsize=24)
plt.setp(ax.get_xticklabels(), fontsize=20)
plt.setp(ax.get_yticklabels(), fontsize=20)

ax.set_xlim([-1., 1.])
ax.set_ylim([0.965, 1.003])
plt.tight_layout()
#plt.show()
plt.savefig('Ch1/transitW431_tim.pdf')#, dpi=250)