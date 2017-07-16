import os
import numpy as np
import h5py
from caput import mpiutil
import hough


psr = 'B1929+10'
fl = 'I.hdf5'

threshold = 3.0
Dl = 0.0 # lower bound of DM
Dh = 100.0 # uper bound of DM
ND = 4000 # number of DM
NC = 2000 # number of time offeset

DM = 4.15
dl = Dl * DM
dh = Dh * DM

indir = '../pulsar_data/%s/' % psr
outdir = '../accumulator/%s/' % psr

with h5py.File(indir+fl, 'r') as f:
    I = f['I'][:]
    freq = f['freq'][:]
    time = f['time'][:]

# # mask
# threshold = 4.0
# med = np.median(I)
# abs_diff = np.abs(I - med)
# mad = np.median(abs_diff) / 0.6745
# Im = np.where(abs_diff>threshold*mad, I-med, np.nan) # subtract median
# # plot Im
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# plt.figure()
# plt.imshow(Im, aspect='auto', cmap='gray')
# plt.colorbar()
# plt.savefig('Im.png')
# plt.close()
# err


I = I.T[:, ::-1] # first time, then freq, both ascending
freq = freq[::-1] # make ascending
freq *= 1.0e-3 # GHz
print freq.shape, freq[0], freq[-1]

A, dl, dh, t0l, t0h = hough.hough_transform(I, time, freq, dl, dh, nd=ND, nt0=NC, threshold=threshold, comm=mpiutil._comm)

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
