import os
# import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm


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
fig = plt.figure()
ax = plt.subplot(111, projection='3d')
# ax.plot_surface(D, C, A, rstride=20, cstride=20, cmap='gray')
# ax.plot_surface(D, C, A, rstride=20, cstride=20, cmap='jet')
# ax.plot_wireframe(D, C, A, rstride=20, cstride=20, color='black')
# ax.plot_wireframe(D, C, A, rstride=80, cstride=80, color='black')
# ax.plot_wireframe(D, C, A, rstride=80, cstride=80, cmap='jet')
# ax.plot_surface(D[600:1600], C[600:1600], A[600:1600], rstride=80, cstride=80, cmap='jet', alpha=0.3)
D = D[600:1600]
C = C[600:1600]
A = A[600:1600]
ax.scatter(D, C, A, s=2, c=A, cmap='jet', edgecolors='none')
# surf = ax.plot_surface(D, C, A, rstride=80, cstride=80, cmap='jet', shade=False)
# ax.plot_surface(D[600:1600], C[600:1600], A[600:1600], rstride=80, cstride=80, linewidth=0.0, alpha=0.25)
# surf.set_facecolor((0,0,0,0))
# ax.plot_surface(D, C, A, rstride=80, cstride=80, cmap='gray')
# ax.plot_surface(D, C, A, cmap='binary')
# ax.contour(D, C, A, cmap='gray', stride=5)
# ax.set_xlim(D.min(), D.max())

ax.set_xlim(D.min(), D.max())
ax.set_ylim(C.min(), C.max())
ax.set_zlim(A.min(), A.max())
# ax.axis('off')
ax.grid(False)

# st = 1
# D = D[::st, ::st].reshape(-1)
# C = C[::st, ::st].reshape(-1)
# A = A[::st, ::st].reshape(-1)
# p = ax.scatter(D, C, A, s=2, c=A, cmap='hot')
# fig.colorbar(p)
# ax.set_zlim(-100, 4000)
# ax.view_init(-45, 45)
ax.view_init(30, 30)
ax.set_xlabel('DM / pc cm${}^{-3}$')
ax.set_ylabel('offset / ms')
# plt.colorbar(ax=ax)
plt.savefig(outdir+'A3dc1.png')
plt.close()
