import os
import numpy as np
import pyfits
import h5py


psr = 'B1929+10'
fl = '/home/nch/Data/pulsar_data/%s-00db.fits' % psr
outdir = './%s/' % psr

hdulist = pyfits.open(fl)
# print hdulist.info()
# print repr(hdulist[0].header)
# print repr(hdulist[1].header)
data = hdulist[1].data
# print data.shape
# print data[0][-1].shape

I = data[0][-1][:, 0, :, 0].T
# I = data[1][-1][:, 0, :, 0].T
freq = data[0][12] # MHz, descending
time = 0.001024e3 * np.arange(I.shape[1]) # ms

# save data
if not os.path.exists(outdir):
    os.makedirs(outdir)

with h5py.File(outdir+'I.hdf5', 'w') as f:
    I = f.create_dataset('I', data=I)
    I.attrs['axes'] = '(freq, time)'
    freq = f.create_dataset('freq', data=freq)
    freq.attrs['unit'] = 'MHz'
    time = f.create_dataset('time', data=time)
    time.attrs['unit'] = 'ms'
