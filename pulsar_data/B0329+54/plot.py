import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


fl = 'I.hdf5'
with h5py.File(fl, 'r') as f:
    I = f['I'][:]
    freq = f['freq'][:]
    time = f['time'][:]


freq *= 1.0e-3 # GHz

# plot
plt.figure()
plt.imshow(I, aspect='auto', cmap='gray', extent=(time[0], time[-1], freq[-1], freq[0]))
# plt.colorbar()
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig('I.png')
plt.close()