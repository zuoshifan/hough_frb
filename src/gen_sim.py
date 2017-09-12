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
# Nf = 512 # number of freq
# Nt = 512 # number of time
# Nf = 200 # number of freq
# Nt = 200 # number of time

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
# mu_s = 2.5    # mean of signal
I = sigma_n * np.random.randn(Nt, Nf) # initialize I to noise
# I = np.zeros((Nt, Nf)) # initialize I to noise

for fi in range(Nf):
    ti = np.int(np.around((t[fi] - tl) / dt))
    ti = max(0, ti)
    ti = min(Nt-1, ti)
    val = (mu_s + sigma_s * np.random.randn())
    I[ti, fi] += val
    # for i in range(-1, 2):
    #     if 0 <= ti+i and ti+i <= Nt-1:
    #         I[ti+i, fi] += val

# save generated data
if not os.path.exists(outdir):
    os.mkdir(outdir)

# plot I
plt.figure()
plt.imshow(I.T, origin='lower', aspect='auto', extent=(tl, th, f[0], f[-1]), cmap='gray')
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
# plt.colorbar()
plt.savefig(outdir+'I.png')
plt.close()

# mask data based on a given threshold
threshold = 3.0
med = np.median(I)
abs_diff = np.abs(I - med)
mad = np.median(abs_diff) / 0.6745
Im = np.where(abs_diff>threshold*mad, I-med, np.nan) # subtract median

# plot Im
plt.figure()
plt.imshow(Im.T, origin='lower', aspect='auto', extent=(tl, th, f[0], f[-1]), cmap='gray')
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
# plt.colorbar()
plt.savefig(outdir+'Im.png')
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
