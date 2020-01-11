import numpy as np


def hough_transform(I, time, freq, dl, dh, nd=None, t0l=None, t0h=None, nt0=None, threshold=3.0, comm=None):
    """Hough transform algorithm.

    Parameters
    ----------
    I : 2D np.ndarray
        Input data image, row time, column freq, both in ascending. order
    time : 1D np.ndarray
        The corresponding time of `I`, in ms.
    freq : 1D np.ndarray
        The corresponding frequency of `I`, in GHz.
    dl : float
        Lowest d = 4.15 DM.
    dh : float
        High d = 4.15 DM.
    nd : integer
        Number of d.
    t0l : float
        Lowest time offset t0, in ms.
    t0h : float
        Highest time offset t0, in ms.
    nt0 : integer
        Number of time offset t0.
    threshold : float
        How many sigmas to truncate the data.
    comm : MPI communicator
        MPI communicator. Required if executed parallelly by using MPI.
    """

    try:
        from caput import mpiutil
    except ImportError:
        # no mpiutil, can not use MPI to speed up
        comm = None

    fl, fh = freq[0], freq[-1]
    Nf = freq.shape[0]
    tl, th = time[0], time[-1]
    Nt = time.shape[0]
    df = (fh - fl) / Nf
    dt = (th - tl) / Nt

    # mask data based on a given threshold
    med = np.median(I)
    if threshold > 0.0:
        abs_diff = np.abs(I - med)
        mad = np.median(abs_diff) / 0.6745
        Im = np.where(abs_diff>threshold*mad, I-med, np.nan) # subtract median
    else:
        Im = I - med

    inds = np.where(np.isfinite(Im)) # inds of non-masked vals

    if nd is None:
        dd = fl**2 * dt
        nd = np.int(np.ceil((dh - dl) / dd))
    d = np.linspace(dl, dh, nd)
    # compute the range of t0
    if t0l is None:
        t0l = -fl**-2 * dh + tl
    if t0h is None:
        t0h = -fh**-2 * dl + th
    assert t0h > t0l, "Must have t0h (= %g) > t0l (= %g)" % (t0h, t0l)
    if nt0 is None:
        nt0 = nd
    dt0 = (t0h - t0l) / nt0

    # initialize the accumulator
    A = np.zeros((nt0, nd)) # the accumulator

    # accumulat the accumulator
    if comm is None or comm.size == 1:
        linds = zip(inds[0], inds[1])
    else:
        linds = mpiutil.mpilist(zip(inds[0], inds[1]), comm=comm)
    for (ti, fi) in linds:
        fv = fl + fi * df
        tv = tl + ti * dt
        t0v = -fv**-2 * d + tv
        for di in xrange(nd):
            t0i = np.int(np.around((t0v[di] - t0l) / dt0))
            t0i = max(0, t0i)
            t0i = min(nt0-1, t0i)
            A[t0i, di] += Im[ti, fi]

    # gather and accumulat A
    if not (comm is None or comm.size == 1):
        A = A.reshape((1, nt0, nd))
        Ac = mpiutil.gather_array(A, comm=comm)
        if mpiutil.rank0:
            A = Ac.sum(axis=0)
        else:
            A = None

    return A, dl, dh, t0l, t0h
