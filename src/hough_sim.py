import os
import numpy as np
import h5py
from caput import mpiutil
import hough


threshold = 3.0
Dl = 800.0 # lower bound of DM
Dh = 1200.0 # uper bound of DM
ND = 4000 # number of DM
NC = 2000 # number of time offeset

DM = 4.15
dl = Dl * DM
dh = Dh * DM

indir = '../sim/'
outdir = '../sim/'

with h5py.File(indir+'I.hdf5', 'r') as f:
    I = f['I'][:] # first time, then freq, both ascending
    fl = f.attrs['fl']
    fh = f.attrs['fh']
    Nf = f.attrs['Nf']
    tl = f.attrs['tl']
    th = f.attrs['th']
    Nt = f.attrs['Nt']

t = tl + (th - tl)/Nt * np.arange(Nt) # ms
f = fl + (fh - fl)/Nf * np.arange(Nf) # GHz

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
