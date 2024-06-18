import matplotlib.pyplot as plt
import numpy as np
import scipy.io.wavfile as wavefile
import pyaudio
import lab4_fft_utils


sample_rate, x = wavefile.read('Hello.wav')  # read audio signal from file
Nsig = np.size(x)

h = np.loadtxt('FIRcoeffs.txt')  # read FIR coefficients from file
Nfilter = np.size(h)

Nfft = ____  # FFT length
L = Nfft - Nfilter + 1  # Overlap-Add frame length
Nframes = int(Nsig / L)  # number of frames

H = lab4_fft_utils.fft(h, Nfft)  # frequency coefficients of the filter

y = np.zeros(Nframes*L+Nfilter-1)
for f in range(Nframes):  # loop over frames
    x_frame = x[range(____, ____)]  # frame of input signal
    X_frame = lab4_fft_utils.fft(x_frame, Nfft)  # FFT over input frame
    Y_frame = np.multiply(____, ____)  # multiply signal and filter in frequency domain
    y_frame = np.real(np.divide(lab4_fft_utils.fft([Y_frame[i] for i in range(0, -Nfft, -1)], Nfft), Nfft))  # IFFT
    y[range(f*L, f*L+Nfft)] = y[range(f*L, f*L+Nfft)] + y_frame  # add output frame to output vector


y = np.zeros(Nframes*L+Nfilter-1)
for f in range(Nframes):  # loop over frames
    x_frame = x[range(f*L, (f+1)*L)]  # frame of input signal
    X_frame = lab4_fft_utils.fft(x_frame, Nfft)  # FFT over input frame
    Y_frame = np.multiply(X_frame, H)  # multiply signal and filter in frequency domain
    y_frame = np.real(np.divide(lab4_fft_utils.fft([Y_frame[i] for i in range(0, -Nfft, -1)], Nfft), Nfft))  # IFFT
    y[range(f*L, f*L+Nfft)] = y[range(f*L, f*L+Nfft)] + y_frame  # add output frame to output vector

# play audio to speakers/headphone
pa = pyaudio.PyAudio()
stream = pa.open(format = pa.get_format_from_width(2), channels = 1, rate = sample_rate, output=True)
stream.write(np.int16(y))
