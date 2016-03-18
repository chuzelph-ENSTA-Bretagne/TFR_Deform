import scipy.interpolate as inter
import numpy as np
import pylab as plt
import tftb

fm, am, iflaw = tftb.generators.misc.doppler(512, 200.0, 65.0, 10.0, 50.0)
sig = np.real(am * fm)
plt.subplot(211), plt.plot(sig)
plt.subplot(212), plt.plot(iflaw)

test = inter.UnivariateSpline(range(512), iflaw)
test.set_smoothing_factor(0.0005)
xnew = np.linspace(0, 512, num=512*10, endpoint=True)
plt.figure()
plt.plot (range(512), iflaw)
plt.plot (xnew, test(xnew),'+')
plt.show()

