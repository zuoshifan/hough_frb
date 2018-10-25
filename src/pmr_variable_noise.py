import os
import numpy as np
import h5py
from caput import mpiutil
import hough
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


outdir = '../pmr/noise/'

# save generated data
if not os.path.exists(outdir):
    os.mkdir(outdir)


DM = 1000.0 # pc cm^-3
d = 4.15 * DM
t0 = 0.0 # time offset

threshold = 3.0
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

# sigma_n = 1.0 # std of noise
sigma_s = 3.0 # std of signal
mu_s = 3.0    # mean of signal

# list of std of noise
# sigma_ns = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
# sigma_ns = [ 1.0 * i for i in range(1, 11) ]
sigma_ns = [ 0.5 * i for i in range(1, 13) ]
pmrs = []
for sigma_n in sigma_ns:
    I = np.random.normal(loc=0.0, scale=sigma_n, size=(Nt, Nf)) # initialize I to noise

    for fi in range(Nf):
        ti = np.int(np.around((t[fi] - tl) / dt))
        ti = max(0, ti)
        ti = min(Nt-1, ti)
        I[ti, fi] += np.random.normal(loc=mu_s, scale=sigma_s)

    # # plot I
    # plt.figure()
    # plt.imshow(I.T, origin='lower', aspect='auto', extent=(tl, th, f[0], f[-1]), cmap='gray')
    # plt.xlabel('Time [ms]')
    # plt.ylabel('Frequency [GHz]')
    # # plt.colorbar()
    # plt.savefig(outdir+'I_%.1f.png' % sigma_n)
    # plt.close()

    A, dl, dh, Cl, Ch = hough.hough_transform(I, time, freq, dl, dh, nd=ND, nt0=NC, threshold=threshold, comm=mpiutil._comm)

    dD = (Dh - Dl) / ND
    dC = (Ch - Cl) / NC

    # peak finding
    p = np.argmax(A)
    Cp = p / ND
    Dp = p % ND
    print 'peak at: (DM, offset) = (%g, %g), with value %g' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p])
    # compute peak-to-median ratio (PMR)
    mx = np.max(A)
    med = np.median(A[A>0])
    pmr = mx / med
    pmrs.append(pmr)
    print 'peak-to-median ratio: %g' % pmr
    f = open(outdir+'peak_%.1f.txt' % sigma_n, 'w')
    f.write('peak at: (DM, offset) = (%g, %g), with value %g\n' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p]))
    f.write('peak-to-median ratio: %g\n' % pmr)
    f.close()

    # plot A
    plt.figure()
    plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch, Cl], cmap='gray', vmax=min(mx, 2.0*med))
    plt.plot(Dl+Dp*dD, Cl+Cp*dC, 'r+', markersize=15.0)
    plt.xlim(Dl, Dh)
    plt.ylim(Ch, Cl)
    plt.xlabel('DM / pc cm${}^{-3}$')
    plt.ylabel('offset / ms')
    # plt.colorbar()
    plt.savefig(outdir+'A_%.1f.png' % sigma_n)
    plt.close()

f = open(outdir+'pmr.txt', 'w')
f.write('sigma_n: %s\n' % sigma_ns)
f.write('pmr: %s\n' % pmrs)
f.close()