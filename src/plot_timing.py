import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# frb = 'FRB010125'
frb = 'FRB110220'
outdir = '../timing/%s/' % frb

# create dir
if not os.path.exists(outdir):
    os.makedirs(outdir)

# for FRB010125
if frb == 'FRB010125':
    tms = [ 802943.875000, 2276113.250000, 1466629.250000, 783429.875000, 320572.000000, 117574.398438, 33982.300781, 24825.900391, 15664.099609, 15677.599609, 15674.700195, 15669.599609 ] # us
# for FRB110220
elif frb == 'FRB110220':
    tms = [ 9173782.000000, 24444358.000000, 12563451.000000, 7944341.500000, 4444145.500000, 2148044.750000, 1101158.500000, 453303.406250, 304178.906250, 254014.500000, 254026.000000, 254057.093750 ] # us
else:
    raise ValueError('Unknown %s' % frb)

tms = np.array(tms)
tms *= 1.0e-3 # to ms
ths = [ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 ]

# plot
fig, ax = plt.subplots()
plt.axhline(tms[0], label='brute force')
plt.plot(ths, tms[1:], 'g-', label='Hought transform')
plt.plot(ths, tms[1:], 'ro')
plt.xlim(-0.2, 5.2)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plt.xlabel(r'$\tau$')
plt.ylabel('Time [ms]')
plt.legend()
plt.savefig(outdir+'timing.png')
plt.close()