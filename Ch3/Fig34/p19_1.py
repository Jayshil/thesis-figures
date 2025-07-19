import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as con
import astropy.units as u
import matplotlib
import seaborn as sns
import os
import plotstyles
import h5py

# ------------------------------------------------------------------------------
#
#                  To plot MIRI transit and occ dep spectra
#
# ------------------------------------------------------------------------------


p1 = os.getcwd() + '/Ch3/Fig34/2_Planetary_Spectra'
p1_mod = os.getcwd() + '/Ch3/Fig34/3_Models'

# --------------- Loading the data
f1 = h5py.File(p1 + '/sparta.h5')

## Transit depths
tra_dep, tra_dep_perr, tra_dep_nerr = np.asarray( f1['transitDepth'] ), np.asarray( f1['transitDepth_errorPos'] ), np.asarray( f1['transitDepth_errorNeg'] )

## Planetary flux ratios
fpfs, fpfs_perr, fpfs_nerr = np.asarray( f1['fp_fs'] ), np.asarray( f1['fp_fs_errorPos'] ), np.asarray( f1['fp_fs_errorNeg'] )

## Phases and wavelengths
phases, waves = np.asarray( f1['phase'] ), np.asarray( f1['wavelength'] )

# Defining colors from colormap (will define 10 colors -- and will choose 2 out of them)
chex = np.array([])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)
cmap = sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
for i in range(fpfs.shape[1]):
    c1 = cmap(norm(100*i))     # Use cm.viridis for viridis colormap and so on...
    c2 = matplotlib.colors.rgb2hex(c1, keep_alpha=True)
    chex = np.hstack((chex, c2))
# 0 1 2 3 4 5 6 7 8 9 10


# --------------------------------------------------------
# ------------------- Planck Function --------------------
# --------------------------------------------------------
def planck_func(lam, temp):
    """
    Given the wavelength and temperature
    this function will compute the specific
    intensity using the Planck's law
    """
    coeff1 = (2*con.h*con.c*con.c)/(lam**5)
    expo = np.exp((con.h*con.c)/(lam*con.k_B*temp)) - 1
    planck = (coeff1/expo).to(u.W / u.m**2 / u.micron)
    return planck

def bb_occ_dep(rprs, temp_pl, temp_st):
    lams = np.linspace(1.,20.,10000)*u.micron
    pl_bb = planck_func(lam=lams, temp=temp_pl * u.K)
    st_bb = planck_func(lam=lams, temp=temp_st * u.K)
    return lams, ( ( rprs**2 ) * 1e6 * pl_bb / st_bb ).decompose()

# --------------------------------------------------------
# --------------------------------------------------------

lam_00, fpfs_00 = bb_occ_dep(rprs=0.1615, temp_pl=850, temp_st=4400)
lam_25, fpfs_25 = bb_occ_dep(rprs=0.1615, temp_pl=1250, temp_st=4400)
lam_50, fpfs_50 = bb_occ_dep(rprs=0.1615, temp_pl=1550, temp_st=4400)
lam_75, fpfs_75 = bb_occ_dep(rprs=0.1615, temp_pl=1150, temp_st=4400)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(2*16/2.7, 9/2.7))

for i in range(len(tra_dep)):
    yers = np.array([tra_dep_nerr[i]/1e4, tra_dep_nerr[i]/1e4]).reshape((2,1))
    axs[0].errorbar(waves[i], tra_dep[i]/1e4, yerr=yers, color=chex[i], fmt='.', mfc='white', elinewidth=2, capthick=2, capsize=3, zorder=10)
axs[0].axhline(np.mean(tra_dep/1e4), color='gray', lw=0.7, zorder=5)

axs[0].set_xlim([5, 10.5])
axs[0].set_ylim([2.44, 2.54])

axs[0].text(5.5, 2.455, 'Band averaged\ntransit depth', color='gray')

axs[0].set_yticks(ticks=np.array([2.44, 2.48, 2.52]))

#axs.set_xlabel('Wavelength [$\mu$m]')
axs[0].set_ylabel('Transit depth [%]')

# --------------------------------------------------------
# --------------------------------------------------------


# Phase 0 (night side)
for i in range(fpfs.shape[1]):
    yers = np.array([ fpfs_nerr[0,i]/1e3, fpfs_perr[0,i]/1e3 ]).reshape((2,1))
    axs[1].errorbar(waves[i], fpfs[0,i]/1e3, yerr=yers, color=chex[i], fmt='.', mfc='white', elinewidth=2, capthick=2, capsize=3, zorder=10)
#axs.plot(lam_00.value, fpfs_00/1e3, color='gray', lw=0.7, zorder=5)

# Phase 0.25 (evening side)
for i in range(fpfs.shape[1]):
    yers = np.array([ fpfs_nerr[1,i]/1e3, fpfs_perr[1,i]/1e3 ]).reshape((2,1))
    axs[1].errorbar(waves[i], fpfs[1,i]/1e3, yerr=yers, color=chex[i], fmt='.', mfc='white', elinewidth=2, capthick=2, capsize=3, zorder=10)
#axs.plot(lam_25.value, fpfs_25/1e3, color='gray', lw=0.7, zorder=5)

# Phase 0.5 (day side)
for i in range(fpfs.shape[1]):
    yers = np.array([ fpfs_nerr[2,i]/1e3, fpfs_perr[2,i]/1e3 ]).reshape((2,1))
    axs[1].errorbar(waves[i], fpfs[2,i]/1e3, yerr=yers, color=chex[i], fmt='.', mfc='white', elinewidth=2, capthick=2, capsize=3, zorder=10)
#axs.plot(lam_50.value, fpfs_50/1e3, color='gray', lw=0.7, zorder=5)

# Phase 0.75 (morning side)
for i in range(fpfs.shape[1]):
    yers = np.array([ fpfs_nerr[3,i]/1e3, fpfs_perr[3,i]/1e3 ]).reshape((2,1))
    axs[1].errorbar(waves[i], fpfs[3,i]/1e3, yerr=yers, color=chex[i], fmt='.', mfc='white', elinewidth=2, capthick=2, capsize=3, zorder=10)
#axs.plot(lam_75.value, fpfs_75/1e3, color='gray', lw=0.7, zorder=5)

clrs = np.array([ 'gray', 'orangered', 'darkorange', 'cornflowerblue' ])
for i in range(len(phases)):
    wav, mod = np.loadtxt(p1_mod + '/models_phase' + str(phases[i]) + '.txt', usecols=(0,1), unpack=True)
    idx = np.argsort(wav)
    axs[1].plot(wav[idx], mod[idx]/1e3, color=clrs[i], lw=0.8, zorder=5)

axs[1].set_xlim([5, 10.5])
axs[1].set_ylim([0, 7.85])

#axs.text(5.5, 8000., 'Black body at T = 1550 K', color='gray')

#axs.set_yticks(ticks=np.array([4000, 5000, 6000, 7000, 8000]))

#axs.set_xlabel('Wavelength [$\mu$m]')
axs[1].set_ylabel(r'F$_p$/F$_\star$ [$\times$ 10$^3$ ppm]')
fig.supxlabel('Wavelength [$\mu$m]', y=0.1)

sns.despine()

plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig34/w43_miri_spectra.pdf')