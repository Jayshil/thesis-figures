import numpy as np
import matplotlib.pyplot as plt
import os
import juliet
from exotoolbox import plots
import matplotlib.gridspec as gd
from astropy.stats import mad_std
from scipy.signal import medfilt
from astropy.io import fits
from astropy.table import Table

import multiprocessing
multiprocessing.set_start_method('fork')

# ------------------------------------------------------------------------------
#
#                    To analyse KELT-11 b transit light curve
#
# ------------------------------------------------------------------------------

def standard(x):
    return ( x - np.nanmedian(x) ) / mad_std(x)

pout = os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis'

# Loading the data
hdul = fits.open(os.getcwd() + '/Ch3/Fig37/KELT11b/CH_PR300024_TG000101_TU2020-03-09T14-50-41_SCI_COR_Lightcurve-DEFAULT_V0300.fits')
tab = Table.read(hdul[1])

times, flux, flux_err = np.asarray( tab['BJD_TIME'] ), np.asarray( tab['FLUX'] ), np.asarray( tab['FLUXERR'] )
background, roll_angle = np.asarray( tab['BACKGROUND'] ), np.asarray( tab['ROLL_ANGLE'] )
cenx, ceny = np.asarray( tab['CENTROID_X'] ), np.asarray( tab['CENTROID_Y'] )
conta, smear = np.asarray( tab['CONTA_LC'] ), np.asarray( tab['SMEARING_LC'] )

# Masking high background points
idx_highbg = np.ones(len(background), dtype=bool)
idx_highbg[background > 4e5] = False

times, flux, flux_err = times[idx_highbg], flux[idx_highbg], flux_err[idx_highbg]
background, roll_angle = background[idx_highbg], roll_angle[idx_highbg]
cenx, ceny = cenx[idx_highbg], ceny[idx_highbg]
conta, smear = conta[idx_highbg], smear[idx_highbg]

med_flux = np.nanmedian(flux)
flux, flux_err = flux / med_flux, flux_err / med_flux

# Masking other outliers
d = abs(medfilt(flux-1, 11)+1-flux)
mad = d.mean()
ok = d < 5*mad

print('\nRejected {} points more than {:0.1f} x MAD = {:0.0f} '
        'ppm from the median'.format(sum(~ok),5,1e6*mad*5))

times, flux, flux_err = times[ok], flux[ok], flux_err[ok]
background, roll_angle = background[ok], roll_angle[ok]
cenx, ceny = cenx[ok], ceny[ok]
conta, smear = conta[ok], smear[ok]

roll = np.radians(roll_angle)

lin_vecs = np.vstack([ np.sin(roll), np.sin(2*roll), np.sin(3*roll), np.sin(4*roll), np.sin(5*roll),\
                       np.cos(roll), np.cos(2*roll), np.cos(3*roll), np.cos(4*roll), np.cos(5*roll),\
                       standard(background), standard(cenx), standard(ceny), standard(conta), standard(smear) ])
lin_vecs = np.transpose( lin_vecs )

# ---------------------------------------------------------------------------
#
#          Preparing the data so that juliet can understand
#
# ---------------------------------------------------------------------------

instrument = 'CHEOPS'
tim, fl, fle = {}, {}, {}
lin_pars, gp_pars = {}, {}

tim[instrument], fl[instrument], fle[instrument] = times, flux, flux_err
lin_pars[instrument], gp_pars[instrument] = lin_vecs, times

# ---------------------------------------------------------------------------
#
#                               Priors
#
# ---------------------------------------------------------------------------

per, per_err = 4.7360990, np.sqrt(0.0000290**2 + 0.0000270**2)
tc, tc_err = 2458553.81381, 0.00033
rho, rho_err = 0.104*1e3, 0.003*1e3

cycle = round((tim[instrument][0] - tc)/per)
tc1, tc1_err = tc + (cycle*per), np.sqrt(tc_err**2 + (cycle*per_err)**2)

par_P = ['P_p1', 't0_p1', 'p_p1', 'b_p1', 'q1_' + instrument, 'q2_' + instrument , 'ecc_p1', 'omega_p1', 'rho']
dist_P = ['fixed', 'normal', 'uniform', 'uniform', 'uniform', 'uniform', 'fixed', 'fixed', 'normal']
hyper_P = [per, [tc1, tc1_err], [0., 1.], [0., 1.], [0., 1.], [0., 1.], 0., 90., [rho, rho_err]]

par_ins = ['mdilution_' + instrument, 'mflux_' + instrument, 'sigma_w_' + instrument]
dist_ins = ['fixed', 'normal', 'loguniform']
hyper_ins = [1., [0., 0.1], [0.1, 1e4]]

par_lin, dist_lin, hyper_lin = [], [], []
for j in range(lin_pars[instrument].shape[1]):
    par_lin.append('theta' + str(j) + '_' + instrument)
    dist_lin.append('uniform')
    hyper_lin.append([-1., 1.])

par_gp = ['GP_S0_' + instrument, 'GP_omega0_' + instrument, 'GP_Q_' + instrument]
dist_gp = ['loguniform', 'loguniform', 'fixed']
hyper_gp = [[np.exp(-30), np.exp(0)], [np.exp(-2.3), np.exp(8)], 1/np.sqrt(2)]

par_tot = par_P + par_ins + par_lin + par_gp
dist_tot = dist_P + dist_ins + dist_lin + dist_gp
hyper_tot = hyper_P + hyper_ins + hyper_lin + hyper_gp

priors_tot = juliet.utils.generate_priors(par_tot, dist_tot, hyper_tot)

# ---------------------------------------------------------------------------
#
#                               Fitting
#
# ---------------------------------------------------------------------------

# And fitting
dataset = juliet.load(priors=priors_tot, t_lc=tim, y_lc=fl, yerr_lc=fle,\
                      GP_regressors_lc=gp_pars, linear_regressors_lc=lin_pars, out_folder=pout)
res = dataset.fit(sampler = 'dynamic_dynesty', nthreads=8)

# ---------------------------------------------------------------------------
#
#                               Plotting
#
# ---------------------------------------------------------------------------


def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
    data = np.ma.masked_invalid(np.ma.copy(data))
    
    # bin data into multiple bin sizes
    npts = data.size
    if maxnbins is None:
        maxnbins = npts / 10.
    binsz = np.arange(1, maxnbins + binstep, step=binstep, dtype=int)
    nbins = np.zeros(binsz.size, dtype=int)
    rms = np.zeros(binsz.size)
    rmserr = np.zeros(binsz.size)
    for i in range(binsz.size):
        nbins[i] = int(np.floor(data.size / binsz[i]))
        bindata = np.ma.zeros(nbins[i], dtype=float)
        # bin data
        # ADDED INTEGER CONVERSION, mh 01/21/12
        for j in range(nbins[i]):
            bindata[j] = np.ma.mean(data[j * binsz[i]:(j + 1) * binsz[i]])
        # get rms
        rms[i] = np.sqrt(np.ma.mean(bindata ** 2))
        rmserr[i] = rms[i] / np.sqrt(2. * nbins[i])
    # expected for white noise (WINN 2008, PONT 2006)
    stderr = (np.ma.std(data) / np.sqrt(binsz)) * np.sqrt(nbins / (nbins - 1.))
    if isrmserr is True:
        return rms, stderr, binsz, rmserr
    else:
        return rms, stderr, binsz


model = res.lc.evaluate(instrument)

# Let's make sure that it works:
fig = plt.figure(figsize=(16/1.5,9/1.5))
gs = gd.GridSpec(2,1, height_ratios=[2,1])

# Top panel
ax1 = plt.subplot(gs[0])
ax1.errorbar(tim[instrument], fl[instrument], yerr=fle[instrument], fmt='.')#, alpha=0.3)
ax1.plot(tim[instrument], model, c='k', zorder=100)
ax1.set_ylabel('Relative Flux')
ax1.set_xlim(np.min(tim[instrument]), np.max(tim[instrument]))
ax1.xaxis.set_major_formatter(plt.NullFormatter())

# Bottom panel
ax2 = plt.subplot(gs[1])
ax2.errorbar(tim[instrument], (fl[instrument]-model)*1e6, yerr=fle[instrument]*1e6, fmt='.')#, alpha=0.3)
ax2.axhline(y=0.0, c='black', ls='--', zorder=100)
ax2.set_ylabel('Residuals (ppm)')
ax2.set_xlabel('Time (BJD)')
ax2.set_xlim(np.min(tim[instrument]), np.max(tim[instrument]))
#plt.show()
plt.savefig(pout + '/full_model_' + instrument + '.png')

residuals = fl[instrument] - model
rms, stderr, binsz = computeRMS(residuals, binstep=1)
normfactor = 1e-6

fig = plt.figure(figsize=(8,6))
plt.plot(binsz, rms / normfactor, color='black', lw=1.5,
                label='Fit RMS', zorder=3)
plt.plot(binsz, stderr / normfactor, color='red', ls='-', lw=2,
                label=r'Std. Err. ($1/\sqrt{N}$)', zorder=1)
plt.xlim(0.95, binsz[-1] * 2)
plt.ylim(stderr[-1] / normfactor / 2., stderr[0] / normfactor * 2.)
plt.xlabel("Bin Size (N frames)", fontsize=14)
plt.ylabel("RMS (ppm)", fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
#plt.show()
plt.savefig(pout + '/alan_deviation_' + instrument + '.png')

plots.corner_plot(pout)