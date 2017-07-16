import os
import numpy as np
import h5py
from caput import mpiutil
import hough


frb = 'FRB140514'
fl = '1.raw.cube'

threshold = 3.0
Dl = 400.0 # lower bound of DM
Dh = 800.0 # uper bound of DM
ND = 4000 # number of DM
NC = 2000 # number of time offeset

DM = 4.15
dl = Dl * DM
dh = Dh * DM

indir = '../frb_data/%s/' % frb
outdir = '../accumulator/%s/' % frb

d = np.genfromtxt(indir+fl, dtype=[('f', 'i4'), ('t', 'i4'), ('I', 'f4')])
# print d.dtype
# print d.shape
I = d['I'].reshape(1024, 1024)
I = I.T[:, ::-1] # first time, then freq, both ascending
print I.shape

t = 2.219968e3 / 1024 * np.arange(1024) # ms
f = np.genfromtxt(indir+'freq.txt', dtype='f4', delimiter=',')
print f.shape, f[0], f[-1]
f = f[::-1] # make ascending
f *= 1.0e-3 # GHz

A, dl, dh, t0l, t0h = hough.hough_transform(I, t, f, dl, dh, nd=ND, nt0=NC, threshold=threshold, comm=mpiutil._comm)
# A, dl, dh, t0l, t0h = hough.hough_transform(I, t, f, dl, dh, threshold=3.0, comm=mpiutil._comm)

if mpiutil.rank0:
    print A.shape
    # save data
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with h5py.File(outdir+'A.hdf5', 'w') as f:
        f.create_dataset('A', data=A)
        f.attrs['axes'] = '(offset, DM)'
        f.attrs['Dl'] = dl / DM
        f.attrs['Dh'] = dh / DM
        f.attrs['ND'] = A.shape[1]
        f.attrs['Cl'] = t0l
        f.attrs['Ch'] = t0h
        f.attrs['NC'] = A.shape[0]
