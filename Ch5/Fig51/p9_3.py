import numpy as np
import matplotlib.pyplot as plt
from astropy import table, units as u
import matplotlib
import matplotlib.cm as cm
import plotstyles

# ------------------------------------------------------------------------------
#
#               To plot Radius-Period diagram for close-in planets
#
# ------------------------------------------------------------------------------

tbl = table.Table.read('Ch5/Fig51/PSCompPars_2025.05.06_06.09.28.votable')

# Assumed Bond Albedo
ab = 0.

# Removing all planets with no period and not enough precision on the radius
max_prec = 0.3
idx_rad = np.isfinite(tbl["pl_radj"])
idx_per = np.isfinite(tbl["pl_orbper"])
idx_a1 = np.isfinite(tbl['pl_orbsmax'])

idx_rad_prec = ((np.abs(tbl["pl_radjerr1"]) < tbl["pl_radj"]*max_prec)
                & (np.abs(tbl["pl_radjerr2"]) < tbl["pl_radj"]*max_prec))

idx_per_prec = ((np.abs(tbl["pl_orbpererr1"]) < tbl["pl_orbper"]*max_prec)
                & (np.abs(tbl["pl_orbpererr2"]) < tbl["pl_orbper"]*max_prec))

idx_ok = idx_rad & idx_per & idx_rad_prec & idx_per_prec & idx_a1

# Period of the planet (x-axis)
per = np.asarray( tbl['pl_orbper'][idx_ok] ) * u.d
per_err_up = np.abs(np.asarray( tbl['pl_orbpererr1'][idx_ok] )) * u.d
per_err_lo = np.abs(np.asarray( tbl['pl_orbpererr2'][idx_ok] )) * u.d

# Radius of the planet
rp = np.asarray( tbl['pl_rade'][idx_ok] ) * u.R_earth
rp_err_up = np.abs(np.asarray( tbl['pl_radeerr1'][idx_ok] )) * u.R_earth
rp_err_lo = np.abs(np.asarray( tbl['pl_radeerr2'][idx_ok] )) * u.R_earth

# Stellar Radii
rst = np.array(tbl['st_rad'][idx_ok]) * u.R_sun
rst_err_up = np.abs(np.array(tbl['st_raderr1'][idx_ok])) * u.R_sun
rst_err_lw = np.abs(np.array(tbl['st_raderr2'][idx_ok])) * u.R_sun

# Semi-major axis
a1 = np.array(tbl['pl_orbsmax'][idx_ok]) * u.au
a1_err_up = np.abs(np.array(tbl['pl_orbsmaxerr1'][idx_ok])) * u.au
a1_err_lw = np.abs(np.array(tbl['pl_orbsmaxerr2'][idx_ok])) * u.au

# Stellar Effective Temperature
teff = np.array(tbl['st_teff'][idx_ok]) * u.K
teff_err_up = np.abs(np.array(tbl['st_tefferr1'][idx_ok])) * u.K
teff_err_lw = np.abs(np.array(tbl['st_tefferr2'][idx_ok])) * u.K

# Equilibrium temperature calculation
tm1 = 1 - ab
tm2 = (rst**2) * (teff**4) / (4 * a1 * a1)
temp_eq = ((tm1 * tm2)**0.25).to(u.K)

len(np.where(np.isnan(temp_eq))[0])

#idx1 = (per < 1 * u.d) & ( rp>10 * u.R_earth )
#print(np.asarray(tbl['pl_name'][idx_ok][idx1]))


#per, rp, teq = per.value, rp.value, temp_eq.value
#per_err, rp_err = np.sqrt((per_err_lo.value)**2 + (per_err_up.value)**2),\
#    np.sqrt((rp_err_lo.value)**2 + (rp_err_up.value)**2)
#"""
#usps = ['55 Cnc e', 'TOI-561 b', 'CoRoT-7 b', 'K2-141 b']#, 'LHS 3844 b']
usps = ['55 Cnc e', 'TOI-561 b', 'HD 189733 b', 'K2-141 b', 'WASP-189 b', 'WASP-18 b']#, 'LHS 3844 b']
arrowprops = dict(
    arrowstyle="->",
    connectionstyle="angle,angleA=0,angleB=90,rad=0",
    lw=1.5)

teq_lbls = np.arange(0,3500,500)

# Figure
fig, axs = plt.subplots(figsize=(16/1.5, 9/1.5))
axs.scatter(per.value, rp.value, c=temp_eq.value, s=0)

norm = matplotlib.colors.Normalize(vmin=100., vmax=3000., clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap='plasma')
vcolor = np.array([(mapper.to_rgba(v)) for v in temp_eq.value])

for x, xerr1, xerr2, y, yerr1, yerr2, clr in zip(per.value, per_err_lo.value, per_err_up.value, rp.value, rp_err_lo.value, rp_err_up.value, vcolor):
    axs.errorbar(x, y, xerr=np.vstack((xerr1, xerr2)), yerr=np.vstack((yerr1, yerr2)), fmt='o', c=clr, mew=0.7, zorder=10, mec='k')

axs.errorbar(np.array([87.969, 224.701, 365.256, 686.980, 4332.589, 10755.699, 30685.400, 60189.018]), np.array([0.383, 0.949, 1., 0.532, 11.209, 9.449, 4.007, 3.883]), color='k', fmt='s', mfc='white', mew=1.5)

for i in range(len(usps)):
    per_usp, rp_usp = per[tbl['pl_name'][idx_ok] == usps[i]].value, rp[tbl['pl_name'][idx_ok] == usps[i]].value
    if usps[i] == 'CoRoT-7 b':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(per_usp-0.085, rp_usp+0.65), arrowprops=arrowprops, zorder=100, fontsize=14)
    elif usps[i] == 'LHS 3844 b':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(per_usp-0.005, rp_usp+0.8), arrowprops=arrowprops, zorder=100, fontsize=14)
    elif usps[i] == '55 Cnc e':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(per_usp-0.1, rp_usp+0.85), arrowprops=arrowprops, zorder=100, fontsize=14)
    elif usps[i] == 'HD 189733 b':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(1., 25.), arrowprops=arrowprops, zorder=100, fontsize=14)
    elif usps[i] == 'WASP-18 b':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(0.7, 8.), arrowprops=arrowprops, zorder=100, fontsize=14)
    elif usps[i] == 'WASP-189 b':
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(1.55, 7.4), arrowprops=arrowprops, zorder=100, fontsize=14)
    else:
        axs.annotate(text=usps[i], xy=(per_usp, rp_usp), xytext=(per_usp-0.1, rp_usp+1.), arrowprops=arrowprops, zorder=100, fontsize=14)

mapper.set_array([])
cbar = plt.colorbar(mapper)
cbar.set_label('Equilibrium Temperature [K]', fontsize=20, rotation=270, labelpad=30)
cbar.set_ticklabels(fontsize=16, ticklabels=teq_lbls)

axs.set_xlabel("Period [d]", fontsize=20)
axs.set_ylabel("Radius of the planet [$R_\oplus$]", fontsize=20)

#axs.set_xlim([1e-1, 1e3])
#axs.set_xlim([1e-1, 7e4])
#axs.set_ylim([0.25, 30.])

axs.set_xlim([0.1, 3.])
axs.set_ylim([0.4, 30.])

axs.set_xscale('log')
axs.set_yscale('log')

axs.set_xticks(ticks=np.array([0.1, 1.]),\
               labels=np.array(['0.1', '1']))
#axs.set_xticks(ticks=np.array([0.1, 1., 10, 1e2, 1e3, 1e4]),\
#               labels=np.array(['0.1', '1', '10', r'10$^2$', r'10$^3$', r'10$^4$']))
axs.set_yticks(ticks=np.array([1., 10.]),\
               labels=np.array(['1', '10']))

plt.xticks(fontsize=16);
plt.yticks(fontsize=16);

plt.tight_layout()

#plt.show()
plt.savefig('PhDTalk/population_close_in_planets.png', dpi=500)