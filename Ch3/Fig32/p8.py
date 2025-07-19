import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import plotstyles
import os

# ------------------------------------------------------------------------------
#
#                   To plot TESS and CHEOPS response functions
#
# ------------------------------------------------------------------------------


wavt, trt = np.loadtxt(os.getcwd() + '/Data/tess_response_fun.txt', usecols=(0,1), unpack=True)
wavc, trc = np.loadtxt(os.getcwd() + '/Data/cheops_response_fun.txt', usecols=(0,1), unpack=True)

wavt = ( wavt * u.AA ).to(u.micron)
wavc = ( wavc * u.AA ).to(u.micron)

fig, axs = plt.subplots(figsize=(16/2.5,9/2.5))

axs.plot(wavt, trt/np.max(trt), color='orangered')#, label='TESS')
axs.plot(wavc, trc/np.max(trc), color='cornflowerblue')#, label='CHEOPS')

axs.set_xlim(np.min(wavc.value), np.max(wavt.value))
axs.set_ylim(0., 1.1)

axs.set_xlabel(r'Wavelength [$\mu$m]', fontsize=16)
axs.set_ylabel('Response function', fontsize=16)

plt.setp(axs.get_xticklabels(), fontsize=15)
plt.setp(axs.get_yticklabels(), fontsize=15)

axs.text(0.62, 0.2, 'TESS bandpass', color='orangered', fontsize=12)
axs.text(0.62, 0.1, 'CHEOPS bandpass', color='cornflowerblue', fontsize=12)

#plt.legend(loc='best', framealpha=0., fontsize=12)#, color='white')
plt.tight_layout()
#plt.show()
plt.savefig('Ch3/Fig32/res_func.pdf')
#plt.savefig('PhDTalk/res_func.png', dpi=500, transparent=True)