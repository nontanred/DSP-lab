import numpy as np
import matplotlib.pyplot as plt
import lab4_fft_utils
import scipy.io.wavfile as wavefile

sample_rate, x = wavefile.read('Twinkle.wav')
sig_length = np.size(x)

N = 2**10  # fft length – STFT frame length
N_frames = int(sig_length/N)  # number of fft frames in STFT
y_frames = np.zeros([N, N_frames])  # output of STFT is a 2-d array

for l in range(N_frames):  # loop over frames
    x_frames = x[range(l * N, (l + 1) * N)]
    y_frames[:, l] = np.abs(lab4_fft_utils.fft(x_frames, N))

ax = plt.gca()
ax.imshow(20*np.log10(np.maximum(y_frames, np.max(y_frames)/1e4)), aspect='auto',
          extent=[0, sig_length/sample_rate, sample_rate, 0])
# set y-axis limit to 0 – 3000 Hz
ax.set(ylim=[0, 3000], xlabel='t [sec]', ylabel='f [Hz]', title='Spectrogram')

