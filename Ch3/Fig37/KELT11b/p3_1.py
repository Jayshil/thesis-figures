import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
from utils import planetary_model, linear_model, gp_model
import matplotlib
import seaborn as sns
import juliet
import plotstyles
import os

# ------------------------------------------------------------------------------
#
#                To plot Stellar, instrumental, and planetary models
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

# ------------------------------------------------------------------------------------------
#
#                               Loading the juliet
#
# ------------------------------------------------------------------------------------------

data = juliet.load(input_folder=os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis/')
res = data.fit(sampler = 'dynamic_dynesty')

model = res.lc.evaluate('CHEOPS')

# ------------------------------------------------------------------------------------------
#
#                               Generating the posteriors
#
# ------------------------------------------------------------------------------------------

# Loading the posteriors
post1 = res.lc.posteriors#['posterior_samples']


# Loading the posteriors
N = int(len(post1['t0_p1'])/2)
per = 4.7360990
tc, rprs = np.random.choice(post1['t0_p1'], size=N, replace=False) , np.random.choice(post1['p_p1'], size=N, replace=False)
b, rho = np.random.choice(post1['b_p1'], size=N, replace=False), np.random.choice(post1['rho'], size=N, replace=False)
q1, q2 = np.random.choice(post1['q1_CHEOPS'], size=N, replace=False), np.random.choice(post1['q2_CHEOPS'], size=N, replace=False)
s0, w0 = np.random.choice(post1['GP_S0_CHEOPS'], size=N, replace=False) , np.random.choice(post1['GP_omega0_CHEOPS'], size=N, replace=False)
mflx, sigw = np.random.choice(post1['mflux_CHEOPS'], size=N, replace=False), np.random.choice(post1['sigma_w_CHEOPS'], size=N, replace=False)*1e-6

# ------------------------------------------------------------------------------------------
#
#                             Loading the light curve data
#
# ------------------------------------------------------------------------------------------

# Loding the light curve data
times, fluxes, fl_errs = np.loadtxt(os.getcwd() + '/Ch3/Fig37/KELT11b/Analysis/lc.dat', usecols=(0,1,2), unpack=True)
times_highres = np.linspace(np.min(times), np.max(times), 1000)

times_hr = ( times - times[0] ) * 24
times_highres_hr = ( times_highres - times[0] ) * 24

# ------------------------------------------------------------------------------------------
#
#                               Generating the models
#
# ------------------------------------------------------------------------------------------

# Generating the models
print('>>>> --- Planetary model (high cadence)')
fl_samps, fl_med, fl_up, fl_lo = planetary_model(tim=times_highres, per=per, tc=tc, rprs=rprs,\
                                                 b=b, q1=q1, q2=q2, rho=rho)

print('>>>> --- Linear model')
lin_samples, lin_med, lin_up, lin_lo = linear_model(lin_pars=data.lm_lc_arguments['CHEOPS'],\
                                                    parameters=post1, N=N)

print('>>>> --- Planetary model (original cadence)')
fl_samps_tim, _, _, _ = planetary_model(tim=times, per=per, tc=tc, rprs=rprs,\
                                       b=b, q1=q1, q2=q2, rho=rho)

resids_samps = np.zeros(fl_samps_tim.shape)
res_errs_samps = np.zeros(fl_samps_tim.shape)
for i in range(resids_samps.shape[1]):
    resids_samps[:,i] = fluxes - ( fl_samps_tim[:,i] / (1 + mflx[i]) ) - lin_samples[:,i]
    res_errs_samps[:,i] = np.sqrt(fl_errs**2 + sigw[i]**2)

print('>>>> --- GP model')
gp_samples, gp_med, gp_up, gp_lo = gp_model(times=times, residuals_samps=resids_samps, res_errs_samps=res_errs_samps, s0=s0, w0=w0)


# ------------------------------------------------------------------------------------------
#
#                           Plotting the models: quantile models
#
# ------------------------------------------------------------------------------------------

fig = plt.figure(figsize=(15/2.7, 18/2.7))
gs = gd.GridSpec(4, 1, height_ratios=[3,1,1,3], hspace=0.1)

ax00 = plt.subplot(gs[0])
ax00.errorbar(times_hr, (fluxes-1)*1e4, yerr=fl_errs*1e4, fmt='.', color='dodgerblue', zorder=1)
ax00.plot(times_hr, (model-1)*1e4, color='navy', lw=2., zorder=10)
ax00.xaxis.set_major_formatter(plt.NullFormatter())

ax0 = plt.subplot(gs[1])
ax0.plot(times_hr, lin_med*1e4, lw=1., color=chex[1])
#ax0.plot(times_hr, (lin_med+lin_up)*1e4, lw=0.5, color=chex[1])
#ax0.plot(times_hr, (lin_med-lin_lo)*1e4, lw=0.5, color=chex[1])
ax0.fill_between(times_hr, y1=(lin_med-lin_lo)*1e4, y2=(lin_med+lin_up)*1e4, color=chex[1], alpha=0.3)

ax0.xaxis.set_major_formatter(plt.NullFormatter())

ax1 = plt.subplot(gs[2])
ax1.plot(times_hr, gp_med*1e4, lw=1.5, color=chex[5])
#ax1.plot(times_hr, (gp_med+gp_up)*1e4, lw=0.3, color=chex[5])
#ax1.plot(times_hr, (gp_med-gp_lo)*1e4, lw=0.3, color=chex[5])
ax1.fill_between(times_hr, y1=(gp_med-gp_lo)*1e4, y2=(gp_med+gp_up)*1e4, color=chex[5], alpha=0.3)
ax1.xaxis.set_major_formatter(plt.NullFormatter())

ax2 = plt.subplot(gs[3])
ax2.plot(times_highres_hr, (fl_med-1)*1e4, lw=1.5, color=chex[8])
#ax2.plot(times_highres_hr, (fl_med+fl_up-1)*1e4, lw=0.3, color=chex[8])
#ax2.plot(times_highres_hr, (fl_med-fl_lo-1)*1e4, lw=0.3, color=chex[8])
ax2.fill_between(times_highres_hr, y1=(fl_med-fl_lo-1)*1e4, y2=(fl_med+fl_up-1)*1e4, color=chex[8], alpha=0.3)
ax2.set_xlabel('Time since the beginning [hr]')

fig.supylabel(r'Relative flux [$\times$ 10$^2$ ppm]', x=0.0001)

sns.despine()

plt.tight_layout()
plt.show()
#plt.savefig('Ch3/Fig37/KELT11b/Figures/kel11b_models.pdf')