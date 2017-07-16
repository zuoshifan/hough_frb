import numpy as np
import psrchive
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# fl = '5.ar'
fl = '5.fil'
arch = psrchive.Archive_load(fl)
print arch.get_source()
data = arch.get_data()
print data.shape
subint = arch.get_Integration(0)
print subint
# print subint.get_duration()

# # plot
# plt.figure()
# # plt.imshow(I, extent=(t[0], t[-1], f[0], f[-1]), aspect='auto', cmap='gray')
# # plt.xlabel('Time [ms]')
# # plt.ylabel('Grequency [GHz]')
# int_data = data[0, 0].reshape(96, 32, 500).sum(axis=1)
# plt.imshow(int_data, aspect='auto')
# plt.colorbar()
# plt.savefig('I.png')
# plt.close()


