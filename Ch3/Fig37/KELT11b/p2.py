import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import matplotlib
import os
import plotstyles

# ------------------------------------------------------------------------------
#
#                     To plot raw KELT-11 b transit light curve
#
# ------------------------------------------------------------------------------

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
my_cmap = sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
for i in range(11):
    c1 = my_cmap(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10


# Plotting the raw data
tim, fl, fle = np.loadtxt(os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis/lc.dat', usecols=(0,1,2), unpack=True)
tim = ( tim - tim[0] ) * 24

# Raw data
fig, axs = plt.subplots(figsize=(16/3, 9/3))
axs.errorbar(tim, (fl-1.)*1e4, yerr=fle*1e4, fmt='.', color='dodgerblue')

axs.set_xlim([ np.min(tim), np.max(tim) ])

axs.set_xlabel('Time since beginning [hr]')
#axs.set_ylabel(r'Relative flux [$\times$ 10$^2$ ppm]')
fig.supylabel(r'Relative flux [$\times$ 10$^2$ ppm]', x=0.07)

sns.despine()

plt.tight_layout()

#plt.show()
plt.savefig('Ch3/Fig37/KELT11b/Figures/kelt11b_raw.pdf')