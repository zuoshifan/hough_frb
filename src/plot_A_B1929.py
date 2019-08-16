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


psr = 'B1929+10'
indir = '../accumulator/%s/' % psr
outdir = '../plot/%s/' % psr

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
# p = np.argmax(np.abs(A))
Cp = p / ND
Dp = p % ND
print 'peak at: (DM, offset) = (%g, %g), with value %g' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p])
if args.save_peak:
    f = open(outdir+'peak.txt', 'w')
    f.write('peak at: (DM, offset) = (%g, %g), with value %g\n' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p]))
    f.close()

# plot A
plt.figure()
# A = np.abs(A) # negative peak in original data
# crop A to display better
A = A[60:1960]
Ch1 = Cl + 1960 * dC
Cl1 = Cl + 60 * dC
plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch1, Cl1], cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
if args.colorbar:
    plt.colorbar()
val = 228
plt.scatter(Dl+Dp*dD, Cl+Cp*dC, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val*6, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val*5, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val*4, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val*3, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val*2, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC-val, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC+val, marker='+', c='r', s=60, linewidth=1.5)
plt.scatter(Dl+Dp*dD, Cl+Cp*dC+val*2, marker='+', c='r', s=60, linewidth=1.5)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val*6, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val*5, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val*4, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val*3, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val*2, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC-val, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC+val, 'r+', markersize=15.0)
# plt.plot(Dl+Dp*dD, Cl+Cp*dC+val*2, 'r+', markersize=15.0)
# plt.scatter(Dl+Dp*dD, Cl+Cp*dC, s=50, c='r', marker='+')
plt.xlim(Dl, Dh)
# plt.ylim(Ch, Cl)
plt.ylim(Ch1, Cl1)
plt.xlabel('DM / pc cm${}^{-3}$')
plt.ylabel('offset / ms')
plt.savefig(outdir+'Ac1.png')
plt.close()