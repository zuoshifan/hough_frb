import os
import numpy as np
import h5py
from caput import mpiutil
import hough


outdir = '../pmr/threshold/'

# save generated data
if not os.path.exists(outdir):
    os.makedirs(outdir)


DM = 1000.0 # pc cm^-3
d = 4.15 * DM
t0 = 0.0 # time offset

Dl = 800.0 # lower bound of DM
Dh = 1200.0 # uper bound of DM
ND = 4000 # number of DM
NC = 2000 # number of time offeset
dl = 4.15 * Dl
dh = 4.15 * Dh

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

time = tl + (th - tl)/Nt * np.arange(Nt) # ms
freq = fl + (fh - fl)/Nf * np.arange(Nf) # GHz

sigma_n = 1.0 # std of noise
sigma_s = 3.0 # std of signal
mu_ss = [ 0.5, 1.0, 2.0, 3.0 ] # mean of signal
thresholds = [ 0.5 * i for i in range(24, -1, -1) ] # tau


pmrs = []
for i, mu_s in enumerate(mu_ss):
    pmrs.append([])

    I = np.random.normal(loc=0.0, scale=sigma_n, size=(Nt, Nf)) # initialize I to noise

    for fi in range(Nf):
        ti = np.int(np.around((t[fi] - tl) / dt))
        ti = max(0, ti)
        ti = min(Nt-1, ti)
        I[ti, fi] += np.random.normal(loc=mu_s, scale=sigma_s)

    for j, threshold in enumerate(thresholds):

        if mpiutil.rank0:
            print '%d of %d...' % (i*len(mu_ss)+j, len(mu_ss)*len(thresholds))

        A, dl, dh, Cl, Ch = hough.hough_transform(I, time, freq, dl, dh, nd=ND, nt0=NC, threshold=threshold, comm=mpiutil._comm)

        # compute peak-to-median ratio (PMR)
        if mpiutil.rank0:
            mx = np.max(A)
            med = np.median(A[A>0])
            pmr = mx / med
            pmrs[i].append(pmr)

if mpiutil.rank0:
    with h5py.File(outdir+'pmr_taus_mus.hdf5', 'w') as f:
        f.create_dataset('pmr', data=np.array(pmrs))
        f.create_dataset('mu_s', data=np.array(mu_ss))
        f.create_dataset('tau', data=thresholds)