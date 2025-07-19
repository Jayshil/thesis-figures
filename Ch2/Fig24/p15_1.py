import numpy as np
import matplotlib.pyplot as plt
from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS import physical_constants as cst
from petitRADTRANS.physics import temperature_profile_function_guillot_global
import astropy.constants as con
import astropy.units as u
import plotstyles

# ------------------------------------------------------------------------------
#
#                  To plot models of transmission spectra 
#
# ------------------------------------------------------------------------------

co2_aubun = [1e-4, 1e-5]


# -------------------------------------------------------
#
#                 Physical quantities
#
# -------------------------------------------------------

infrared_mean_opacity = 0.01
gamma = 0.4
intrinsic_temperature = 200
equilibrium_temperature = 1166

## Planetary and stellar parameters
mp, rp = 0.281 * u.M_jup, 1.279 * u.R_jup
mst, rst = 0.913 * u.M_sun, 0.939 * u.R_sun

dep_white = ( 0.1457 ** 2 ) * 1e2

g = ( con.G *  mp / ( rp ** 2) ).to(u.cm / u.s**2)

planet_radius = rp.value * cst.r_jup_mean       # This is in cm
reference_gravity = 10 ** np.log10(g.value)     # This is in CGS units
reference_pressure = 1

# -------------------------------------------------------
# -------------------------------------------------------
#
#                 Radtrans object
#
# -------------------------------------------------------

radtrans = Radtrans(
    pressures=np.logspace(-6, 2, 100),
    line_species=[
        'H2O',
        'CO-NatAbund',
        'CO2',
        'Na',
        'K',
        'SO2'
    ],
    rayleigh_species=['H2', 'He'],
    gas_continuum_contributors=['H2-H2', 'H2-He'],
    wavelength_boundaries=[0.4, 7]
)


# -------------------------------------------------------
# -------------------------------------------------------
#
#                 Temperature profile
#
# -------------------------------------------------------

pressures = radtrans.pressures * 1e-6 # cgs to bar
temperatures = temperature_profile_function_guillot_global(
    pressures=pressures,
    infrared_mean_opacity=infrared_mean_opacity,
    gamma=gamma,
    gravities=reference_gravity,
    intrinsic_temperature=intrinsic_temperature,
    equilibrium_temperature=equilibrium_temperature
)


# -------------------------------------------------------
# -------------------------------------------------------
#
#                      Abundaces
#
# -------------------------------------------------------

mass_fractions_all = []

for i in range(len(co2_aubun)):
    mass_fractions = {
        'H2': 0.74 * np.ones(temperatures.size),
        'He': 0.24 * np.ones(temperatures.size),
        'H2O': 1e-3 * np.ones(temperatures.size),
        'CO-NatAbund': 1e-2 * np.ones(temperatures.size),
        'CO2': co2_aubun[i] * np.ones(temperatures.size),
        'SO2': 5e-5 * np.ones(temperatures.size),
        'Na': 1e-6 * np.ones(temperatures.size),
        'K': 1e-6 * np.ones(temperatures.size)
    }
    mass_fractions_all.append(mass_fractions)

mean_molar_masses = 2.33 * np.ones(temperatures.size)

# -------------------------------------------------------
# -------------------------------------------------------
#
#                 And transit depth
#
# -------------------------------------------------------

transit_depths_all = []

for i in range(len(mass_fractions_all)):
    wavelengths, transit_radii, _ = radtrans.calculate_transit_radii(
        temperatures=temperatures,
        mass_fractions=mass_fractions_all[i],
        mean_molar_masses=mean_molar_masses,
        reference_gravity=reference_gravity,
        planet_radius=planet_radius,
        reference_pressure=reference_pressure
    )

    transit_depth = ( ( transit_radii * u.cm / rst ).decompose() )**2 * 1e2

    transit_depths_all.append( transit_depth )


# -------------------------------------------------------
# -------------------------------------------------------
#
#                 Figures
#
# -------------------------------------------------------
colors = ['deeppink', 'dodgerblue']
lbls = ['1e-4', '1e-5']

fig, ax = plt.subplots(figsize = (15/1.5, 5/1.5))

for i in range(len(transit_depths_all)):
    ax.plot(wavelengths * 1e4, transit_depths_all[i], color=colors[i], label=lbls[i])

#ax.axhline(dep_white, ls='--', color='k')
ax.set_xscale('log')
ax.set_xlabel(r'Wavelength [$\mu$m]')
ax.set_ylabel(r'Transit depth [%]')

ax.set_xticks(ticks=np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4., 5., 6., 7., 8.]),\
              labels=np.array(['', '0.5', '', '0.7', '', '', '1', '', '3', '', '5', '', '7', '']))

ax.spines[['right', 'top']].set_visible(False)

ax.text(4.2, 2.02, r'CO$_2$ mass fraction = 10$^{-4}$', color='deeppink')
ax.text(4.2, 1.98, r'CO$_2$ mass fraction = 10$^{-5}$', color='dodgerblue')

#ax.set_facecolor('transparent')

plt.tight_layout()
#plt.show()
plt.savefig('Ch2/Fig24/trans_spec_model.pdf')#, dpi=500, transparent=True)