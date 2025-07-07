import numpy as np
import matplotlib.pyplot as plt
import plotstyles


# ------------------------------------------------------------------------------
#
#                    To plot three spectral timeseries
#
# ------------------------------------------------------------------------------


def gaus(x, amp, mu, std):
    return amp * np.exp(-0.5 * ( (x - mu) / (std))**2 )


l1, l2 = 588.995095*10, 589.592424*10
x = np.linspace(l1-10, l2+10, 10000)
vdop, c = 551 * 2e2, 299798452  # both is m/s

shift1, shift2 = vdop * l1 / c, -vdop * l1 / c

spec_base = -gaus(x, 1., l1, 1.) - gaus(x, 0.5, l2, 1.)
spec_blue = -gaus(x, 1., l1 + shift2, 1.) - gaus(x, 0.5, l2 + shift2, 1.)
spec_red = -gaus(x, 1., l1 + shift1, 1.) - gaus(x, 0.5, l2 + shift1, 1.)

fig, axs = plt.subplots(figsize=(16/2, 9/2))
axs.plot(x, spec_base, 'k-')
axs.plot(x, spec_blue+1.2, ls='-', color='cornflowerblue')
axs.plot(x, spec_red-1.2, ls='-', color='orangered')

axs.text(x[100], spec_base[100]-0.3, 'Baseline', fontsize=20)
axs.text(x[100], spec_blue[100]+1.2-0.3, 'Blueshift', fontsize=20, color='navy')
axs.text(x[100], spec_red[100]-1.2-0.3, 'Redshift', fontsize=20, color='maroon')

axs.set_yticks(ticks=np.array([]), labels=np.array([]))

axs.set_xlabel('Wavelength [Å]', fontsize=24)
axs.set_ylabel('Arbitrary units', fontsize=24)

plt.setp(axs.get_xticklabels(), fontsize=20)
plt.setp(axs.get_yticklabels(), fontsize=20)

plt.tight_layout()

#plt.show()
plt.savefig('Ch1/rv_spec_W43.pdf')#, dpi=250)