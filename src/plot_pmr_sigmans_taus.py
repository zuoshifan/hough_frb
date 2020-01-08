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


indir = '../pmr/noise/'
outdir = '../pmr/noise/'

with h5py.File(indir+'pmr_sigmans_taus.hdf5', 'r') as f:
    pmr = f['pmr'][:]
    sigma_n = f['sigma_n'][:]
    tau = f['tau'][:]

colors = [ 'r', 'g', 'b', 'k' ]

# plot
fig, ax = plt.subplots()
for i in range(4):
    if args.type == 'lin':
        plt.plot(sigma_n, pmr[:, i], '%s' % colors[i], label=r'$\tau \, = \, %.1f$' % tau[i])
        plt.plot(sigma_n, pmr[:, i], '%so' % colors[i])
    elif args.type == 'log':
        plt.semilogy(sigma_n, pmr[:, i], '%s' % colors[i], label=r'$\tau \, = \, %.1f$' % tau[i])
        plt.semilogy(sigma_n, pmr[:, i], '%so' % colors[i])
    elif args.type == 'dB':
        plt.plot(sigma_n, 10*np.log10(pmr[:, i]), '%s' % colors[i], label=r'$\tau \, = \, %.1f$' % tau[i])
        plt.plot(sigma_n, 10*np.log10(pmr[:, i]), '%so' % colors[i])
plt.legend()
plt.xlim(0.3, 6.2)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plt.xlabel(r'$\sigma_n$', fontsize=16)
if args.type == 'dB':
    plt.ylabel('PMR / dB', fontsize=16)
else:
    plt.ylabel('PMR', fontsize=16)
plt.savefig(outdir+'pmr_sigmans_tau_%s.png' % args.type)
plt.close()