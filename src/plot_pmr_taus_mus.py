import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


parser = argparse.ArgumentParser(description='Plot accumulator.')
parser.add_argument('-t', '--type', type=str, default='lin', help='Whether to save peak value to file, default save.')
args = parser.parse_args()


indir = '../pmr/threshold/'
outdir = '../pmr/threshold/'


with h5py.File(indir+'pmr_taus_mus.hdf5', 'r') as f:
    mu_s = f['mu_s'][:]
    pmr = f['pmr'][:]
    tau = f['tau'][:]

colors = [ 'r', 'g', 'b', 'k' ]

threshold = tau[12:]
pmr = pmr[:, 12:]

# plot
fig, ax = plt.subplots()
for i in range(4):
    if args.type == 'lin':
        plt.plot(threshold, pmr[i], '%s' % colors[i], label=r'$\mu \, = \, %.1f$' % mu_s[i])
        plt.plot(threshold, pmr[i], '%so' % colors[i])
    elif args.type == 'log':
        plt.semilogy(threshold, pmr[i], '%s' % colors[i], label=r'$\mu \, = \, %.1f$' % mu_s[i])
        plt.semilogy(threshold, pmr[i], '%so' % colors[i])
    elif args.type == 'dB':
        plt.plot(threshold, 10*np.log10(pmr[i]), '%s' % colors[i], label=r'$\mu \, = \, %.1f$' % mu_s[i])
        plt.plot(threshold, 10*np.log10(pmr[i]), '%so' % colors[i])
plt.legend()
plt.xlim(-0.2, 6.2)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plt.xlabel(r'$\tau$', fontsize=16)
if args.type == 'dB':
    plt.ylabel('PMR / dB', fontsize=16)
else:
    plt.ylabel('PMR', fontsize=16)
plt.savefig(outdir+'pmr_taus_mus_%s.png' % args.type)
plt.close()