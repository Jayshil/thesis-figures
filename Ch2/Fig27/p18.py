import numpy as np
import matplotlib.pyplot as plt
import plotstyles

def gaussian(x, mu, sigma):
    amp = 1 / np.sqrt(2 * np.pi * sigma * sigma)
    exp = np.exp(-0.5 * ( (x-mu)**2 ) / sigma / sigma )
    return amp * exp

x1 = np.linspace(-0.5, 0.5, 1000)
y1 = gaussian(x=x1, mu=0., sigma=0.1)

plt.figure(figsize=(16/2, 9/2))
plt.plot(x1, y1, 'k-')
#plt.show()
plt.savefig('Ch2/Fig27/pos_gauss.pdf')

plt.figure(figsize=(16/2, 9/2))
plt.plot(x1, -y1, 'k-')
#plt.show()
plt.savefig('Ch2/Fig27/neg_gauss.pdf')