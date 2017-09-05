import os
import time
import numpy as np
import h5py
import brute


def hough_transform(I, time, freq, dl, dh, nd=None, t0l=None, t0h=None, nt0=None, threshold=3.0):
    fl, fh = freq[0], freq[-1]
    Nf = freq.shape[0]
    tl, th = time[0], time[-1]
    Nt = time.shape[0]
    df = (fh - fl) / Nf
    dt = (th - tl) / Nt

    # mask data based on a given threshold
    med = np.median(I)
    if threshold > 0.0:
        abs_diff = np.abs(I - med)
        mad = np.median(abs_diff) / 0.6745
        Im = np.where(abs_diff>threshold*mad, I-med, np.nan) # subtract median
    else:
        Im = I - med

    inds = np.where(np.isfinite(Im)) # inds of non-masked vals

    if nd is None:
        dd = fl**2 * dt
        nd = np.int(np.ceil((dh - dl) / dd))
    d = np.linspace(dl, dh, nd)
    # compute the range of t0
    if t0l is None:
        t0l = -fl**-2 * dh + tl
    if t0h is None:
        t0h = -fh**-2 * dl + th
    # assert t0h > t0l, "Must have t0h (= %g) > t0l (= %g)" % (t0h, t0l)
    if nt0 is None:
        nt0 = nd
    dt0 = (t0h - t0l) / nt0

    # initialize the accumulator
    A = np.zeros((nt0, nd)) # the accumulator

    # accumulat the accumulator
    for (ti, fi) in zip(inds[0], inds[1]):
        fv = fl + fi * df
        tv = tl + ti * dt
        t0v = -fv**-2 * d + tv
        for di in xrange(nd):
            t0i = np.int(np.around((t0v[di] - t0l) / dt0))
            t0i = max(0, t0i)
            t0i = min(nt0-1, t0i)
            A[t0i, di] += Im[ti, fi]

    return A, dl, dh, t0l, t0h


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
Nt = t.shape[0]
f = np.genfromtxt(indir+'freq.txt', dtype='f4', delimiter=',')
print f.shape, f[0], f[-1]
f = f[::-1] # make ascending
f *= 1.0e-3 # GHz
Nf = f.shape[0]


N = 10

tms = []
# timing for brute force method
dts = []
for i in range(N):
    t1 = time.time()
    brute.dedisp(I, t, f, dl, dh, ND)
    t2 = time.time()
    dts.append(t2 - t1)
dt_mean = np.sum(dts) / N
print dt_mean
tms.append(dt_mean)


# timing for Hough transform method

I = I.T # first time, then freq, both ascending

thresholds = [ 0.0, 1.0, 2.0, 3.0 ]
# thresholds = [ 3.0 ]
for td in thresholds:
    dts = []
    for i in range(N):
        t1 = time.time()
        hough_transform(I, t, f, dl, dh, nd=ND, nt0=Nt, threshold=td)
        t2 = time.time()
        dts.append(t2 - t1)
    dt_mean = np.sum(dts) / N
    print dt_mean
    tms.append(dt_mean)

print tms

# A, dl, dh, t0l, t0h = hough_transform(I, t, f, dl, dh, nd=ND, threshold=threshold)

# # plot
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

# plt.figure()
# plt.imshow(A, aspect='auto', extent=[Dl, Dh, t0h, t0l], cmap='gray')
# plt.savefig('A_test.png')
# plt.close()