import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


fl = '3.4ms.cube'
d = np.genfromtxt(fl, dtype=[('f', 'i4'), ('t', 'i4'), ('I', 'f4')])
# print d.dtype
# print d.shape
# print d['I'][:5]
I = d['I'].reshape(1024, 500)
print I.shape

t = 4.0 * np.arange(500) # ms
f = np.genfromtxt('freq.txt', dtype='f4', delimiter=',')
print f.shape, f[0], f[-1]
f = f[::-1] # make ascending
f *= 1.0e-3 # GHz

# plot
plt.figure()
# plt.imshow(I.T[:,::-1], extent=(f[0], f[-1], t[-1], t[0]), aspect='auto', origin='lower')
plt.imshow(I, extent=(t[0], t[-1], f[0], f[-1]), aspect='auto', cmap='gray')
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig('I.png')
plt.close()
