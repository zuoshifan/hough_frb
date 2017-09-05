import numpy as np


def shift(a, n, pad=0):
    """Shift an array `a`. Right shift if n>0, left shift if n<0."""
    a1 = np.empty_like(a)
    if n > 0:
        a1[:n] = pad
        a1[n:] = a[:-n]
    elif n < 0:
        a1[:n] = a[-n:]
        a1[n:] = pad
    else:
        a1[:] = a

    return a1


def dedisp(I, time, freq, dl, dh, nd=None):
    """Brute force dedispersion. `I` should be first freq, then time and both are ascending."""
    fl, fh = freq[0], freq[-1]
    nf = freq.shape[0]
    cf = freq[nf/2] # central frequency
    tl, th = time[0], time[-1]
    nt = time.shape[0]
    dt = (th - tl) / nt

    if nd is None:
        dd = fl**2 * dt
        nd = np.int(np.ceil((dh - dl) / dd))
    ds = np.linspace(dl, dh, nd)

    B = np.zeros((nd, nt)) # to save result of brute force transformation
    for di, d in enumerate(ds):
        I1 = np.zeros_like(I)
        for fi in range(nf):
            delta_t = d * (freq[fi]**-2 - cf**-2) # time difference relative to cf
            I1[fi] = shift(I[fi], -delta_t/dt)
        B[di] = I1.sum(axis=0)

    return B