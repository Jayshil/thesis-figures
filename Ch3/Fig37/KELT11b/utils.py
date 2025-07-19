import numpy as np
import batman
import astropy.constants as con
from exotoolbox.utils import get_quantiles
from celerite.terms import SHOTerm
from celerite import GP
from tqdm import tqdm
import juliet

G = con.G.value

def planetary_model(tim, per, tc, rprs, b, q1, q2, rho):
    # First, compute u1, u2 from q1, q2
    u1, u2 = juliet.utils.reverse_ld_coeffs('quadratic', q1=q1, q2=q2)
    ar = ((rho * G * ((per * 24. * 3600.)**2)) / (3. * np.pi))**(1. / 3.)
    inc = np.rad2deg( np.arccos(b / ar) )

    # Running batman recursively to store all samples
    fl_samples = np.zeros(( len(tim), len(tc) ))
    for i in tqdm(range(len(tc))):
        params = batman.TransitParams()
        params.t0 = tc[i]                       #time of inferior conjunction
        params.per = per                        #orbital period
        params.rp = rprs[i]                     #planet radius (in units of stellar radii)
        params.a = ar[i]                        #semi-major axis (in units of stellar radii)
        params.inc = inc[i]                     #orbital inclination (in degrees)
        params.ecc = 0.                         #eccentricity
        params.w = 90.                          #longitude of periastron (in degrees)
        params.u = [u1[i], u2[i]]               #limb darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"          #limb darkening model

        m1 = batman.TransitModel(params, tim)    #initializes model
        fl_samples[:,i] = m1.light_curve(params) #calculates light curve

    fl_med, fl_up, fl_lo = np.zeros(len(tim)), np.zeros(len(tim)), np.zeros(len(tim))
    for i in tqdm(range(len(tim))):
        qua = get_quantiles(fl_samples[i,:])
        fl_med[i], fl_up[i], fl_lo[i] = qua[0], qua[1]-qua[0], qua[0]-qua[2]

    return fl_samples, fl_med, fl_up, fl_lo

def linear_model(lin_pars, parameters, N):
    # Number of linear vectors
    nos = lin_pars.shape[1]

    # First selecting random posteriors for thetai_CHEOPS
    theta_random = np.zeros((nos, N))
    for i in tqdm(range(nos)):
        theta_random[i,:] = np.random.choice(parameters['theta' + str(i) + '_CHEOPS'], size=N, replace=False)

    # Computing random samples
    lin_samples = np.zeros((lin_pars.shape[0], N))
    for i in tqdm(range(N)):
        lin_samples[:,i] = lin_pars @ theta_random[:,i]

    # Computing quantiles
    lin_med, lin_up, lin_lo = np.zeros(lin_pars.shape[0]), np.zeros(lin_pars.shape[0]), np.zeros(lin_pars.shape[0])
    for i in tqdm(range(len(lin_med))):
        qua = get_quantiles(lin_samples[i,:])
        lin_med[i], lin_up[i], lin_lo[i] = qua[0], qua[1]-qua[0], qua[0]-qua[2]

    return lin_samples, lin_med, lin_up, lin_lo

def gp_model(times, residuals_samps, res_errs_samps, s0, w0):

    gp_samples = np.zeros( (len(times), len(s0)) )
    for i in tqdm(range(len(s0))):
        # Building a GP model
        kernel = SHOTerm(log_S0=np.log(s0[i]), log_omega0=np.log(w0[i]), log_Q=np.log(1/np.sqrt(2)))
        gp = GP(kernel=kernel, mean=0.)
        gp.compute(times, yerr=res_errs_samps[:,i], check_sorted=False)

        gp_samples[:,i] = gp.predict(y=residuals_samps[:,i], t=times, return_cov=False)

    # Computing quantiles
    gp_med, gp_up, gp_lo = np.zeros(len(times)), np.zeros(len(times)), np.zeros(len(times))
    for i in tqdm(range(len(times))):
        qua = get_quantiles(gp_samples[i,:])
        gp_med[i], gp_up[i], gp_lo[i] = qua[0], qua[1]-qua[0], qua[0]-qua[2]

    return gp_samples, gp_med, gp_up, gp_lo