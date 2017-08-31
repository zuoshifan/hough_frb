import os
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot accumulator.')
parser.add_argument('-m', '--vmax', type=float, default=None, help='Max value of plot.')
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
print 'peak at: (DM, offset) = (%g, %g), with value %g' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p])
f = open(outdir+'peak.txt', 'w')
f.write('peak at: (DM, offset) = (%g, %g), with value %g\n' % (Dl+Dp*dD, Cl+Cp*dC, A.flatten()[p]))
f.close()

# plot A
plt.figure()
plt.imshow(A, aspect='auto', extent=[Dl, Dh, Ch, Cl], cmap='gray', vmax=args.vmax)
plt.xlabel('DM / pc cm${}^{-3}$')
plt.ylabel('offset / ms')
# plt.colorbar()
plt.savefig(outdir+'A_1.png')
plt.close()