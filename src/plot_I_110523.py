import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot image.')
parser.add_argument('-m', '--vmax', type=float, default=None, help='Max value of plot.')
parser.add_argument('-n', '--vmin', type=float, default=None, help='Min value of plot.')
parser.add_argument('-c', '--cmap', type=str, default='jet', help='Color map.')
parser.add_argument('-b', '--colorbar', action='store_false', help='Whether to show colorbar, default show.')
args = parser.parse_args()

frb = '110523'
indir = '../frb_data/FRB%s/' % frb
outdir = '../frb_data/FRB%s/' % frb

fl = 'filtered_short.npy'
I = np.load(indir+fl)
I = I[:, 0] # only I
freq = np.load(indir+'freq.npy')
freq = freq[::-1] # small to large
freq *= 1.0e-3 # GHz
time = np.load(indir+'time_short.npy')
time = 1.0e3 * (time - time.min()) # ms

# plot
plt.figure()
plt.imshow(I, aspect='auto', cmap=args.cmap, extent=(time[0], time[-1], freq[0], freq[-1]), vmin=args.vmin, vmax=args.vmax)
if args.colorbar:
    plt.colorbar()
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig(outdir+'Ic.png')
plt.close()