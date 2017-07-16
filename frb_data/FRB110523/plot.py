import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


indir = './'
outdir = './'

I = np.load(indir+'filtered_short.npy')
# I = I[:, 0][:, ::-1] # only I
I = I[:, 0] # only I
freq = np.load(indir+'freq.npy')
freq = freq[::-1] # small to large
freq *= 1.0e-3 # GHz
time = np.load(indir+'time_short.npy')
time = 1.0e3 * (time - time.min()) # ms

# plot
plt.figure()
plt.imshow(I, extent=(time[0], time[-1], freq[0], freq[-1]), aspect='auto', cmap='gray')
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig('I.png')
plt.close()
