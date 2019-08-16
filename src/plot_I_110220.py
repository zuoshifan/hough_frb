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

frb = '110220'
indir = '../frb_data/FRB%s/' % frb
outdir = '../frb_data/FRB%s/' % frb

fl = '3.4ms.cube'
d = np.genfromtxt(indir+fl, dtype=[('f', 'i4'), ('t', 'i4'), ('I', 'f4')])
# print d.dtype
# print d.shape
# print d['I'][:5]
I = d['I'].reshape(1024, 500)
print I.shape

t = 4.0 * np.arange(500) # ms
f = np.genfromtxt(indir+'freq.txt', dtype='f4', delimiter=',')
print f.shape, f[0], f[-1]
f = f[::-1] # make ascending
f *= 1.0e-3 # GHz

# plot
plt.figure()
plt.imshow(I, aspect='auto', cmap=args.cmap, extent=(t[0], t[-1], f[0], f[-1]), vmin=args.vmin, vmax=args.vmax)
if args.colorbar:
    plt.colorbar()
plt.xlabel('Time [ms]')
plt.ylabel('Frequency [GHz]')
plt.savefig(outdir+'Ic.png')
plt.close()