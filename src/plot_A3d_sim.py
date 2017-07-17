import os
# import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt


# parser = argparse.ArgumentParser(description='Plot accumulator.')
# parser.add_argument('-m', '--vmax', type=float, default=None, help='Max value of plot.')
# args = parser.parse_args()


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

# create dir
if not os.path.exists(outdir):
    os.makedirs(outdir)

D = np.linspace(Dl, Dh, ND)
C = np.linspace(Cl, Ch, NC)
# C, D = np.meshgrid(C, D)
D, C = np.meshgrid(D, C)
print C.shape, D.shape, A.shape

# plot A in 3d
ax = plt.subplot(111, projection='3d')
# ax.plot_surface(D, C, A, rstride=20, cstride=20, cmap='gray')
# ax.plot_wireframe(D, C, A, rstride=20, cstride=20, color='black')
ax.plot_wireframe(D, C, A, rstride=80, cstride=80, color='black')
# ax.set_zlim(-100, 4000)
# ax.view_init(30, -45)
ax.view_init(30, 30)
ax.set_xlabel('DM / pc cm${}^{-3}$')
ax.set_ylabel('offset / ms')
# plt.colorbar(ax=ax)
plt.savefig(outdir+'A3d.png')
plt.close()
