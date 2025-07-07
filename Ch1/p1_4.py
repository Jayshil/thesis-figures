import numpy as np
import matplotlib.pyplot as plt
import radvel
import plotstyles

# ------------------------------------------------------------------------------
#
#                        To plot radial velocity
#
# ------------------------------------------------------------------------------

# Planetary parameters (Patel & Espinoza 2022)
per, tc = 0.8134749, 2458555.80567
rprs, bb = 0.1615, 0.698
aR = 4.72
inc1 = np.rad2deg(np.arccos(bb / aR))
q1, q2 = 0.3106615413, 0.5494710423
K = 551.7

def radvel_rv_model(tim):
    anybasis_params = radvel.Parameters(1,basis='per tc e w k', planet_letters={1: 'b'})    # initialize Parameters object

    anybasis_params['per1'] = radvel.Parameter(value=per)      # period of 1st planet
    anybasis_params['tc1'] = radvel.Parameter(value=tc)     # time of inferior conjunction of 1st planet
    anybasis_params['e1'] = radvel.Parameter(value=0.)          # eccentricity of 1st planet
    anybasis_params['w1'] = radvel.Parameter(value=np.pi/2.)      # argument of periastron of the star's orbit for 1st planet
    anybasis_params['k1'] = radvel.Parameter(value=K)          # velocity semi-amplitude for 1st p

    rv_radvel = radvel.model.RVModel(anybasis_params).__call__(tim)
    return rv_radvel

tim = np.linspace(tc-1.5*per, tc+1.5*per, 1000)

fig, axs = plt.subplots(figsize=(16/2, 9/2))
axs.plot(tim-tc, radvel_rv_model(tim), 'k-')

axs.set_xlabel('Time [d]', fontsize=24)
axs.set_ylabel('Radial velocity [m/s]', fontsize=24)

axs.set_ylim([-600., 600.])

plt.setp(axs.get_xticklabels(), fontsize=20)
plt.setp(axs.get_yticklabels(), fontsize=20)


plt.tight_layout()

#plt.show()
plt.savefig('Ch1/RVW431_tim.pdf')#, dpi=250)