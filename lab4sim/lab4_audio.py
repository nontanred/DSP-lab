import numpy as np
import matplotlib.pyplot as plt
import lab4_fft_utils
import scipy.io.wavfile as wavefile

N = 2**16  # fft length
sample_rate, x = wavefile.read('Twinkle.wav')  # read audio file
sig_length = np.size(x)
y = lab4_fft_utils.fft(x, N)  # calculate FFT of the signal

plt.plot(np.abs(y))
plt.show()
