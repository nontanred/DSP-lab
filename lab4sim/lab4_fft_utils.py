import numpy as np


# generates  array of twiddle factor
def twiddle_factor_table(N):
    # the array of twiddle factors W is   exp(-j2*pi*n/N)    where  0<=n<N
    W = np.exp(-1j * 2 * np.pi * np.arange(____) / ____)  # FILL IN
    return W


# bit reverse order of integers from 0 to N-1
def bit_reverse_order(N):
    br_order = np.zeros(N).astype(int)
    K = int(np.log2(N))
    for n in range(N):
        b = '{:0{width}b}'.format(n, width=K)
        br_order[n] = int(b[::-1], 2)
    return br_order


# single basic butterfly. returns x1+x2*W and x1-x2*W
def butterfly(x1, x2, W, N, K, k, n):
    twiddle_factor = W[(n >> (K - k - 1)) << (K - k - 1)]  # select W for the butterfly
    temp = ____ * ____  # FILL IN
    y1 = ____ + ____  # FILL IN
    y2 = ____ - ____  # FILL IN
    return y1, y2


# decimation in time FFT
def fft(x, N, W=None, br_order=None):
    K = ____  # number of stages, calculated from FFT length (N) - FILL IN
    if N != 2 ** K:
        print("FFT length must be power of 2")
        return
    if W is None:
        W = twiddle_factor_table(N)  # generate array of twiddle factors
    if br_order is None:
        br_order = bit_reverse_order(N)  # generate bit reverse indexing
    if np.size(x) < N:
        x = np.pad(x, (0, N-np.size(x)))  # zero-pad to N
    x = [x[n] for n in br_order]  # bit reverse order of input vector (decimation in time)

    for k in range(K):  # loop over FFT stages
        groupsize = 1 << (k + 1)
        numgroups = 1 << (K - k - 1)
        stepsize = 1 << k
        for n in range(N >> 1):  # loop over butterflies in stage
            n1 = (n % numgroups) * groupsize + (n >> (K - k - 1))
            n2 = n1 + stepsize
            x[n1], x[n2] = butterfly(x[n1], x[n2], W, N, K, k, n)
    return x
