import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


outdir = '../sim/'

DM = 1000.0 # pc cm^-3
d = 4.15 * DM
t0 = 0.0 # time offset

fl = 0.4 # lower bound of freq, GHz
fh = 0.8 # uper bound of freq, GHz
Nf = 2048 # number of freq
Nt = 2048 # number of time

f = np.linspace(fl, fh, Nf)
t = d * f**-2
t -= t.min() # make time offset to 0
t += t0 # make time offset to t0
tl = t.min()
th = t.max()
dt = (th - tl) / Nt

sigma_n = 1.0 # sigma of noise
sigma_s = 3.0 # sigma of signal
mu_s = 3.0    # mean of signal
I = sigma_n * np.random.randn(Nt, Nf) # initialize I to noise
# I = np.zeros((Nt, Nf)) # initialize I to noise

for fi in range(Nf):
    ti = np.int(np.around((t[fi] - tl) / dt))
    ti = max(0, ti)
    ti = min(Nt-1, ti)
    I[ti, fi] += (mu_s + sigma_s * np.random.randn())

# save generated data
if not os.path.exists(outdir):
    os.mkdir(outdir)

# plot I
plt.figure()
plt.imshow(I.T, origin='lower', aspect='auto', extent=(t[0], t[-1], f[0], f[-1]), cmap='gray')
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig(outdir+'I.png')
plt.close()

# save I
with h5py.File(outdir+'I.hdf5', 'w') as f:
    f.create_dataset('I', data=I)
    f.attrs['axes'] = '(time, freq)'
    f.attrs['DM'] = DM
    f.attrs['t0'] = t0
    f.attrs['fl'] = fl
    f.attrs['fh'] = fh
    f.attrs['Nf'] = Nf
    f.attrs['tl'] = tl
    f.attrs['th'] = th
    f.attrs['Nt'] = Nt
    f.attrs['sigma_n'] = sigma_n
    f.attrs['sigma_s'] = sigma_s
    f.attrs['mu_n'] = mu_s
