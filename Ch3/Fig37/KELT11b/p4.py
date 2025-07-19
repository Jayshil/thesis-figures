import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as con
import seaborn as sns
import matplotlib
import os
import juliet
import pickle
from glob import glob
import plotstyles

# ------------------------------------------------------------------------------
#
#                        To plot posterior distributions
#
# ------------------------------------------------------------------------------

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
my_cmap = matplotlib.cm.plasma#sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
for i in range(11):
    c1 = my_cmap(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10

G = con.G.value

# Loading the data
fname = glob(os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis/*.pkl')[0]
post = pickle.load(open(fname, 'rb'))
post1 = post['posterior_samples']

tim = np.loadtxt(os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis/lc.dat', usecols=0, unpack=True)
tref = tim[0]

# Planetary parameters
per = 4.7360990
tc, rho = post1['t0_p1'], post1['rho']
u1, u2 = juliet.utils.reverse_ld_coeffs('quadratic', post1['q1_CHEOPS'], post1['q2_CHEOPS'])
rprs, bb = post1['p_p1'], post1['b_p1']
ar = ((rho * G * ((per * 24. * 3600.)**2)) / (3. * np.pi))**(1. / 3.)
inc = np.rad2deg( np.arccos(bb / ar) )

clr = chex[7]

fig, axs = plt.subplots(nrows=3, ncols=2, sharey=True, figsize=(10/3, 16/3))

axs[0,0].hist((tc-tref)*24, bins=100, histtype='step', color=clr, lw=1.5);
axs[0,1].hist(ar, bins=100, histtype='step', color=clr, lw=1.5);

axs[1,0].hist((rprs**2) * 1e6, bins=100, histtype='step', color=clr, lw=1.5);
axs[1,1].hist(inc, bins=100, histtype='step', color=clr, lw=1.5);

axs[2,0].hist(u1, bins=100, histtype='step', color=clr, lw=1.5);
axs[2,1].hist(u2, bins=100, histtype='step', color=clr, lw=1.5);

axs[0,0].yaxis.set_major_formatter(plt.NullFormatter())
axs[1,0].yaxis.set_major_formatter(plt.NullFormatter())
axs[2,0].yaxis.set_major_formatter(plt.NullFormatter())

#fig.supylabel('Density', x=0.05)

plt.tight_layout()

sns.despine()

#plt.show()
plt.savefig('Ch3/Fig37/KELT11b/Figures/posteriors.pdf')