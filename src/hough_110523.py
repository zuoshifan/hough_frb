import os
import numpy as np
import h5py
from caput import mpiutil
import hough

frb = 'FRB110523'
fl = 'filtered_short.npy'

threshold = 3.0
Dl = 500.0 # lower bound of DM
Dh = 700.0 # uper bound of DM
ND = 4000 # number of DM
NC = 2000 # number of time offeset

DM = 4.15
dl = Dl * DM
dh = Dh * DM

indir = '../frb_data/%s/' % frb
outdir = '../accumulator/%s/' % frb

I = np.load(indir+fl)
I = I[:, 0].T[:, ::-1] # only I, first time, then freq, both ascending

# # plot
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# plt.figure()
# plt.imshow(I, aspect='auto', cmap='gray')
# plt.savefig('I523.png')
# plt.close()
# err

f = np.load(indir+'freq.npy')
f = f[::-1] # small to large
f *= 1.0e-3 # GHz
t = np.load(indir+'time_short.npy')
t = 1.0e3 * (t - t.min()) # ms

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
