import numpy as np
import matplotlib.pyplot as plt
import lab4_fft_utils
plt.rcParams['axes.grid'] = True


### 1 - impulse signal
N = 8
x = [0, 1, 0, 0, 0, 0, 0, 0]

# ### 2 - pulse
# N = 16
# n0 = 4
# x = (np.arange(N) < n0) * 1

# ### 3 – another pulse
# N = 16
# n0 = 4
# x = (np.abs(np.mod(np.arange(N)+N/2, N)-N/2) < n0) * 1

# ### 4 - sine wave
# N = 64
# w0 = 2 * np.pi / N
# x = np.sin(w0 * np.arange(N))

y = lab4_fft_utils.fft(x, N)

print('x[n]: ', x)
print('X[k]: ', np.round(y, 5))

fig, ax1 = plt.subplots(2, 1, sharex=True, sharey=True)
ax1[0].stem(np.real(x), use_line_collection=True)
ax1[1].stem(np.imag(x), use_line_collection=True)
ax1[0].set(ylabel='real(x[n])', title='x[n]')
ax1[1].set(ylabel='imag(x[n])', xlabel='n')

fig, ax = plt.subplots(3, 1, sharex=True, sharey=True)
ax[0].stem(np.abs(y), use_line_collection=True)
ax[1].stem(np.real(y), use_line_collection=True)
ax[2].stem(np.imag(y), use_line_collection=True)

ax[0].set(ylabel='abs(X[k])', title='FFT of x[n]')
ax[1].set(ylabel='real(X[k])')
ax[2].set(ylabel='imag(X[k])', xlabel='k')

plt.show()

