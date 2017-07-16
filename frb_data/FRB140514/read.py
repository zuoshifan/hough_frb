import numpy as np
import psrchive
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


fl = '2014-05-14-17:05:42_0000032532221952.000000_01.ar'
arch = psrchive.Archive_load(fl)
print arch.get_source()
data = arch.get_data()
print data.shape
err
# subint = arch.get_Integration(1)
# print subint
# print subint.get_duration()

# plot
plt.figure()
# plt.imshow(I, extent=(t[0], t[-1], f[0], f[-1]), aspect='auto', cmap='gray')
# plt.xlabel('Time [ms]')
# plt.ylabel('Grequency [GHz]')
int_data = data[0, 0].reshape(1024, 50, 625).sum(axis=1)
plt.imshow(int_data, aspect='auto')
plt.colorbar()
plt.savefig('I.png')
plt.close()


