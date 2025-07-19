import numpy as np
import matplotlib.pyplot as plt
import batman
from juliet.utils import get_phases
from exotoolbox.utils import reverse_ld_coeffs
from matplotlib import cm
import matplotlib
import plotstyles

# ------------------------------------------------------------------------------
#
#          To plot transit light curves with different limb darkening laws
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

# Planetary parameters (Patel & Espinoza 2022)
per, tc = 0.8134749, 2458555.80567
rprs, bb = 0.1615, 0.698
aR = 4.72
inc1 = np.rad2deg(np.arccos(bb / aR))

tim = np.linspace(tc-0.15, tc+0.15, 10000)

def transit_light_curve(q1, q2):
    u1, u2 = reverse_ld_coeffs('quadratic', q1, q2)
    #u1, u2 = q1, q2
    print(u1, u2)
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
    m = batman.TransitModel(params, tim) #initializes model
    flux = m.light_curve(params)         #calculates light curve
    return flux


fl1 = transit_light_curve(0., 0.)
fl2 = transit_light_curve(0.15, 0.35)
fl3 = transit_light_curve(0.5, 0.6)

phs = get_phases(tim, per, tc)
fig, ax = plt.subplots(figsize=(16/2, 9/2))
ax.plot(phs, fl1, ls='-', color=chex[2])#, label='No limb darkening')
ax.plot(phs, fl2, ls='-', color=chex[5])#, label='Quadratic limb darkening 1')
ax.plot(phs, fl3, ls='-', color=chex[8])#, label='Quadratic limb darkening 2')

ax.set_yticks(ticks=np.array([1.0, 0.995, 0.990, 0.985, 0.980, 0.975, 0.970]),\
              labels=np.array(['1.0', '', '0.99', '', '0.98', '', '0.97']))
ax.set_xlabel('Orbital phase', fontsize=24)
ax.set_ylabel('Relative flux', fontsize=24)
plt.setp(ax.get_xticklabels(), fontsize=20)
plt.setp(ax.get_yticklabels(), fontsize=20)

ax.set_xlim([-0.05, 0.05])
ax.set_ylim([0.970, 1.003])

ax.text(-0.0143, 0.998, 'No limb darkening', color=chex[2], fontsize=14)#, fontweight='bold')
ax.text(-0.0209, 0.995, 'Quadratic limb darkening 1', color=chex[5], fontsize=14)#, fontweight='bold')
ax.text(-0.0209, 0.992, 'Quadratic limb darkening 2', color=chex[8], fontsize=14)

#ax.legend(loc='best', framealpha=0., fontsize=14)

plt.tight_layout()
#plt.show()
plt.savefig('Ch2/Fig22/transitLD.pdf')#, dpi=250)