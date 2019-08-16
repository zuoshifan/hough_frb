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

psr = 'B2319+60'
indir = '../pulsar_data/%s/' % psr
outdir = '../pulsar_data/%s/' % psr

fl = 'I.hdf5'
with h5py.File(indir+fl, 'r') as f:
    I = f['I'][:]
    freq = f['freq'][:]
    time = f['time'][:]


freq *= 1.0e-3 # GHz

# plot
plt.figure()
plt.imshow(I, aspect='auto', cmap=args.cmap, extent=(time[0], time[-1], freq[-1], freq[0]), vmin=args.vmin, vmax=args.vmax)
if args.colorbar:
    plt.colorbar()
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig(outdir+'Ic.png')
plt.close()