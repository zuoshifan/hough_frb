import os
import numpy as np
import h5py
import brute


frb = 'FRB010125'
fl = '5.4ms.cube'

Dl = 600.0 # lower bound of DM
Dh = 1000.0 # uper bound of DM
ND = 2000 # number of DM

DM = 4.15
dl = Dl * DM
dh = Dh * DM
# ds = np.linspace(dl, dh, ND)

indir = '../frb_data/%s/' % frb
# outdir = '../accumulator/%s/' % frb

data = np.genfromtxt(indir+fl, dtype=[('f', 'i4'), ('t', 'i4'), ('I', 'f4')])
# print data.dtype
# print data.shape
I = data['I'].reshape(96, 500)
I = I[::-1, :] # first freq, then time, both ascending
print I.shape

dt = 4.0 # ms
t = dt * np.arange(I.shape[1]) # ms
f = np.genfromtxt(indir+'freq.txt', dtype='f4', delimiter=',')
print f.shape, f[0], f[-1]
f = f[::-1] # make ascending
f *= 1.0e-3 # GHz

B = brute.dedisp(I, t, f, dl, dh, ND)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# plot B
plt.figure()
plt.imshow(B.T, origin='lower', aspect='auto', extent=[Dl, Dh, 0, dt*I.shape[1]], cmap='gray')
plt.savefig('B.png')
plt.close()

Cl = 0
dC = 4.0
dD = (Dh - Dl) / ND
p = np.argmax(B.T)
Cp = p / ND
Dp = p % ND
print 'peak at: (DM, offset) = (%g, %g), with value %g' % (Dl+Dp*dD, Cl+Cp*dC, B.T.flatten()[p])