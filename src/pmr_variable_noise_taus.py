import os
import numpy as np
import h5py
from caput import mpiutil
import hough


outdir = '../pmr/noise/'

# save generated data
if not os.path.exists(outdir):
    os.mkdir(outdir)


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

sigma_s = 3.0 # std of signal
mu_s = 3.0    # mean of signal
sigma_ns = [ 0.5 * i for i in range(1, 13) ] # list of std of noise
thresholds = [ 1.0, 2.0, 3.0, 4.0 ] # tau

pmrs = []
for i, sigma_n in enumerate(sigma_ns):
    pmrs.append([])

    I = np.random.normal(loc=0.0, scale=sigma_n, size=(Nt, Nf)) # initialize I to noise

    for fi in range(Nf):
        ti = np.int(np.around((t[fi] - tl) / dt))
        ti = max(0, ti)
        ti = min(Nt-1, ti)
        I[ti, fi] += np.random.normal(loc=mu_s, scale=sigma_s)

    for j, threshold in enumerate(thresholds):

        if mpiutil.rank0:
            print '%d of %d...' % (i*len(sigma_ns)+j, len(sigma_ns)*len(thresholds))

        A, dl, dh, Cl, Ch = hough.hough_transform(I, time, freq, dl, dh, nd=ND, nt0=NC, threshold=threshold, comm=mpiutil._comm)

        # compute peak-to-median ratio (PMR)
        if mpiutil.rank0:
            mx = np.max(A)
            med = np.median(A[A>0])
            pmr = mx / med
            pmrs[i].append(pmr)

if mpiutil.rank0:
    with h5py.File(outdir+'pmr_sigmans_taus.hdf5', 'w') as f:
        f.create_dataset('pmr', data=np.array(pmrs))
        f.create_dataset('sigma_n', data=np.array(sigma_ns))
        f.create_dataset('tau', data=thresholds)