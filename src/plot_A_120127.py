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


frb = 'FRB120127'
indir = '../accumulator/%s/' % frb
outdir = '../plot/%s/' % frb

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
print 'peak at: (DM, offset) = (%g, %g), with value %g' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p])
if args.save_peak:
    f = open(outdir+'peak.txt', 'w')
    f.write('peak at: (DM, offset) = (%g, %g), with value %g\n' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p]))
    f.close()

# plot A
plt.figure()
plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch, Cl], cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
if args.colorbar:
    plt.colorbar()
# plt.plot(Dl+Dp*dD, Cl+Cp*dC, 'r+', markersize=15.0)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC, marker='+', c='r', s=60, linewidth=1.5)
plt.xlim(Dl, Dh)
plt.ylim(Ch, Cl)
plt.xlabel('DM / pc cm${}^{-3}$', fontsize=16)
plt.ylabel('offset / ms', fontsize=16)
plt.savefig(outdir+'Ac2.png')
plt.close()