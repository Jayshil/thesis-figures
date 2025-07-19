import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import plotstyles
import matplotlib
import os

# ------------------------------------------------------------------------------
#
#             To plot spectra of rock vapour atmospheres
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

names = np.array(['TOI-1807b', '55Cnce', 'TOI-561b'])
lbls = np.array(['TOI-1807 b', '55 Cnc e', 'TOI-561 b'])
clrs = np.array([chex[8], chex[5], chex[2]])

fig, axs = plt.subplots(figsize=(16/1.5, 9/1.5))
axins = axs.inset_axes([0.1, 0.42, 0.25, 0.55])#, visible=True)
axins.set_yscale('log')

for i in range(len(names)):
    # For spectra:
    wav0, spec0 = np.loadtxt(os.getcwd() + '/Ch5/Fig55/Data/Fp_F00_' + names[i] + '.dat', usecols=(0,1), unpack=True)
    wav6, spec6 = np.loadtxt(os.getcwd() + '/Ch5/Fig55/Data/Fp_F60_' + names[i] + '.dat', usecols=(0,1), unpack=True)
    # Plotting
    axs.plot(wav0, spec0, color=clrs[i], ls='-', zorder=10, label=lbls[i])
    axs.plot(wav6, spec6, color='k', ls='-', zorder=1)

    # Inset axis
    p0, t0 = np.loadtxt(os.getcwd() + '/Ch5/Fig55/Data/TP_F00_' + names[i] + '.txt', usecols=(0,1), unpack=True)
    p6, t6 = np.loadtxt(os.getcwd() + '/Ch5/Fig55/Data/TP_F60_' + names[i] + '.txt', usecols=(0,1), unpack=True)
    ## Converting pressures to bar
    p0, p6 = p0 * 1e-6, p6 * 1e-6

    axins.plot(t0, p0, color=clrs[i], ls='-', zorder=10)
    axins.plot(t6, p6, color=clrs[i], ls='--', zorder=1)

    axins.set_xlabel('Temperature [K]', fontsize=14)
    axins.set_ylabel('Pressure [bar]', fontsize=14)

    axins.invert_yaxis()

    plt.setp(axins.get_xticklabels(), fontsize=12)
    plt.setp(axins.get_yticklabels(), fontsize=12)


axs.set_xlabel(r'Wavelength [$\mu$m]', fontsize=18)
axs.set_ylabel('Occultation depth [ppm]', fontsize=18)
plt.setp(axs.get_xticklabels(), fontsize=16)
plt.setp(axs.get_yticklabels(), fontsize=16)

#axs.legend(loc='lower right', fontsize=15)

axs.text(14, 90, 'TOI-1807 b', color=chex[8], fontsize=16)#, fontweight='bold')
axs.text(14, 60, '55 Cnc e', color=chex[5], fontsize=16)#, fontweight='bold')
axs.text(14, 30, 'TOI-561 b', color=chex[2], fontsize=16)#, fontweight='bold')

axs.text(8.8, 477, 'SiO', color='k', fontsize=16)#, fontweight='bold')
axs.text(6.8, 440, r'SiO$_2$', color='k', fontsize=16, rotation=90)#, fontweight='bold')

x = np.linspace(7.55,12.95,1000)
y1 = np.ones(1000) * 0
y2 = np.ones(1000) * 800
poly = axs.fill_between(x=x, y1=y1, y2=y2, color='none', alpha=0.5)#cmap='plasma')
verts = np.vstack([p.vertices for p in poly.get_paths()])
gradient = plt.imshow(np.linspace(0, -1, 256).reshape(-1, 1), cmap='Greens', aspect='auto',
                      extent=[verts[:, 0].min(), verts[:, 0].max(), verts[:, 1].min(), verts[:, 1].max()], alpha=0.7)
gradient.set_clip_path(poly.get_paths()[0], transform=plt.gca().transData)


x = np.linspace(6.75,7.55,1000)
y1 = np.ones(1000) * 0
y2 = np.ones(1000) * 800
poly = axs.fill_between(x=x, y1=y1, y2=y2, color='none', alpha=0.5)#cmap='plasma')
verts = np.vstack([p.vertices for p in poly.get_paths()])
gradient = plt.imshow(np.linspace(0, -1, 256).reshape(-1, 1), cmap='Purples', aspect='auto',
                      extent=[verts[:, 0].min(), verts[:, 0].max(), verts[:, 1].min(), verts[:, 1].max()], alpha=0.7)
gradient.set_clip_path(poly.get_paths()[0], transform=plt.gca().transData)

axs.set_xscale('log')

axs.set_xticks(ticks=np.array([1., 10.]),\
               labels=np.array(['1', '10']))
axs.set_xlim([np.min(wav0), np.max(wav0)])
axs.set_ylim([-10, 500])

plt.tight_layout()

#plt.show()
plt.savefig('Ch5/Fig55/rock_vapour.pdf')