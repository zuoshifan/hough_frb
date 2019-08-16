import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot accumulator.')
parser.add_argument('-m', '--vmax', type=float, default=None, help='Max value of plot.')
parser.add_argument('-n', '--vmin', type=float, default=None, help='Min value of plot.')
parser.add_argument('-c', '--cmap', type=str, default='jet', help='Color map.')
parser.add_argument('-b', '--colorbar', action='store_false', help='Whether to show colorbar, default show.')
parser.add_argument('-s', '--save_peak', action='store_false', help='Whether to save peak value to file, default save.')
args = parser.parse_args()


indir = '../sim/'
outdir = '../sim/'

# read in data
with h5py.File(indir+'A.hdf5', 'r') as f:
    A = f['A'][:]
    Dl = f.attrs['Dl']
    Dh = f.attrs['Dh']
    ND = f.attrs['ND']
    Cl = f.attrs['Cl']
    Ch = f.attrs['Ch']
    NC = f.attrs['NC']

dD = (Dh - Dl) / ND
dC = (Ch - Cl) / NC

# create dir
if not os.path.exists(outdir):
    os.makedirs(outdir)

# peak finding
p = np.argmax(A)
Cp = p / ND
Dp = p % ND
# compute peak-to-median ratio (PMR)
mx = np.max(A)
med = np.median(A[A>0])
pmr = mx / med
msg = 'peak at: (DM, offset) = (%g, %g), with value %g, PMR %g' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p], pmr)
print msg
# compute peak-to-median ratio (PMR)
pmr = np.max(A) / np.median(A[A>0])
print 'peak-to-median ratio: %g' % pmr
if args.save_peak:
    f = open(outdir+'peak.txt', 'w')
    f.write('%s\n' % msg)
    f.close()

# plot A
plt.figure()
# plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch, Cl], cmap=args.cmap, vmax=args.vmax)
# crop A to display better
A = A[600:1600]
A = np.where(A>150, 150, A)
Ch = Cl + 1600 * dC
Cl = Cl + 600 * dC
plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch, Cl], cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
plt.xlabel('DM / pc cm${}^{-3}$')
plt.ylabel('offset / ms')
if args.colorbar:
    plt.colorbar()
plt.savefig(outdir+'Ac1_crop.png')
plt.close()