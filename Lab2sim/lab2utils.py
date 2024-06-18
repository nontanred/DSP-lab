import numpy as np
from scipy import signal as signal
import lab2nco

def siggen(root):
    dspsmplfreq = root.SMPLFREQ * (1 + root.SMPLFREQALPHA)
    oversmplratio = root.OVERSMPLRATIO
    sigsmplfreq = dspsmplfreq * oversmplratio
    N0 = root.N0
    CHUNK = root.CHUNK
    n = np.arange(N0, N0+CHUNK)
    t = n / sigsmplfreq
    try:
        freq = root.sigfreq.get()
    except:
        freq = 0
    try:
        sigamp = root.sigamp.get()/2
    except:
        sigamp = 0
    amp = sigamp / 1000 / (1 + np.exp((freq - sigsmplfreq) * oversmplratio / sigsmplfreq))
    if root.sigshape.get() == 1:  # sine wave
        sig = np.sin(2 * np.pi * freq * t)
    elif root.sigshape.get() == 2:  # rectangular
        sig = np.sign(np.sin(2 * np.pi * freq * t))
    elif root.sigshape.get() == 3:  # triangular
        sig = abs(np.mod(4 * freq * t - 1, 4) - 2) - 1
    elif root.sigshape.get() == 4:  # saw-tooth
        sig = np.mod(2 * freq * t + 1, 2) - 1
    sig = amp * sig + np.random.normal(0, 0.005, np.size(sig))

    if root.pulsemod.get():
        root.pulsepw.get()
        root.pulsepri.get()
        envelope = np.mod(t*1e6, root.pulsepri.get()) <= root.pulsepw.get()
        sig = sig * envelope

    root.N0 = N0 + CHUNK
    return sig


def measuresig(root, sig1, sig2):
    measstr = ""
    for m in range(root.meascount):
        if root.measch[m] == "ch1":
            sigmeas = sig1
            # sigmeas = signal.lfilter(np.ones(root.triggerwindow)/root.triggerwindow, [1], sig1)
        else:
            sigmeas = sig2
            # sigmeas = signal.lfilter(np.ones(root.triggerwindow)/root.triggerwindow, [1], sig2)
        if root.meastype[m] == "pk-pk amp":
            measvar = np.diff(np.percentile(sigmeas, [0.1, 99.9]))[0] * 1000
            measunit = "mV"
        elif root.meastype[m] == "RMS amp":
            measvar = np.sqrt(np.mean(np.power(sigmeas, 2))) * 1000
            measunit = "mV"
        elif root.meastype[m] == "max amp":
            measvar = np.percentile(sigmeas, 99) * 1000
            measunit = "mV"
        elif root.meastype[m] == "min amp":
            measvar = np.percentile(sigmeas, 1) * 1000
            measunit = "mV"
        elif root.meastype[m] == "freq":
            if np.size(sigmeas) > 100000:
                sigmeas = sigmeas[range(np.size(sigmeas)-100000, np.size(sigmeas))]
            measvar = np.mean(np.diff(np.unwrap(np.angle(signal.hilbert(sigmeas))[range(2, np.size(sigmeas)-2)]))) / (2*np.pi) * root.SMPLFREQ * root.OVERSMPLRATIO
            measunit = "Hz"

        if not(np.isnan(root.measdata[m])):
            measvar = measvar * root.measalpha + root.measdata[m] * (1-root.measalpha)
        root.measdata[m] = measvar
        measstr = measstr + "(" + root.measch[m] + ") " + root.meastype[m] + ":" +\
                  format(measvar, " 9.2f") + " " + measunit + "\n"
    root.meastext.set_text(measstr[:-1])

def updatecursors(root):
    cursorstr = ""
    if root.cursorx1.isvisible():
        x1 = root.cursorx1.getvalue()
        cursorstr = cursorstr + "(x1): " + format(x1, "0.2f") + " usec\n"
    if root.cursorx2.isvisible():
        x2 = root.cursorx2.getvalue()
        cursorstr = cursorstr + "(x2): " + format(x2, "0.2f") + " usec\n"
    if root.cursorx1.isvisible() & root.cursorx2.isvisible():
        dx = x2 - x1
        cursorstr = cursorstr + "(dx): " + format(dx, "0.2f") + " usec\n"
    if root.cursory1.isvisible():
        y1 = root.cursory1.getvalue()
        cursorstr = cursorstr + "(y1): " + format(y1*1000, "0.1f") + " mV\n"
    if root.cursory2.isvisible():
        y2 = root.cursory2.getvalue()
        cursorstr = cursorstr + "(y2): " + format(y2*1000, "0.1f") + " mV\n"
    if root.cursory1.isvisible() & root.cursory2.isvisible():
        dy = y2 - y1
        cursorstr = cursorstr + "(dy): " + format(dy*1000, "0.1f") + " mV\n"
    root.cursortext.set_text(cursorstr[:-1])

def sigplot(root, sig1, sig2):
    if root.triggerlocked == False:
        if root.triggerchslctd.get() == 'ch1':
            sigtrig = sig1
        elif root.triggerchslctd.get() == 'ch2':
            sigtrig = sig2
        if root.sigtrigresvalid:
            sigtrig = np.concatenate((root.sigtrigres, sigtrig))
        Ntrig = np.size(sigtrig)
        sigtrigsmoothed = signal.lfilter(np.ones(root.triggerwindow), [1], sigtrig)
        zcrossed = (sigtrigsmoothed[range(root.triggerwindow-1, Ntrig-1)] < 0) & (sigtrigsmoothed[range(root.triggerwindow, Ntrig)] > 0)
        if np.any(zcrossed):
            trigidx = np.argmax(zcrossed) + (root.triggerwindow>>1)
            root.sigtrigresvalid = False
            root.triggerlocked = True
            sig1 = sig1[range(trigidx, np.size(sig1))]
            sig2 = sig2[range(trigidx, np.size(sig2))]
        else:
            root.sigtrigres = sigtrig[range(Ntrig-root.triggerwindow, Ntrig)]
            root.sigtrigresvalid = True

    if root.triggerlocked == True:
        sigch1buffer = np.concatenate((root.sigch1buffer, sig1))
        sigch2buffer = np.concatenate((root.sigch2buffer, sig2))
        Nsig = np.size(sigch1buffer)
        tbase = root.tbase.get()
        dspsmplfreq = root.SMPLFREQ * (1 + root.SMPLFREQALPHA)
        oversmplratio = root.OVERSMPLRATIO
        sigsmplfreq = dspsmplfreq * oversmplratio
        Tsig = Nsig / sigsmplfreq * 1e6

        if Tsig < tbase:
            root.sigch1buffer = sigch1buffer
            root.sigch2buffer = sigch2buffer
            if tbase < 1000000:
                skipdraw = True
            else:
                skipdraw = False
        else:
            root.sigch1buffer = []
            root.sigch2buffer = []
            root.triggerlocked = False
            skipdraw = False

        if skipdraw == False:
            t = np.arange(Nsig) / sigsmplfreq
            Nmax = 1e4
            N = np.size(t)
            if (N > Nmax):
                decidx = (np.random.rand(N)) < (Nmax/N)
                t = t[decidx]
            else:
                decidx = range(N)
            sigch1plot = sigch1buffer[decidx]
            sigch2plot = sigch2buffer[decidx]
            ch1vdiv = root.ch1vdiv.get()
            ch2vdiv = root.ch2vdiv.get()
            ch1vpos = root.ch1vpos.get()
            ch2vpos = root.ch2vpos.get()
            if root.ch1on.get():
                root.line1.set_data(t*1e6, sigch1plot*ch1vdiv + ch1vpos)
            else:
                root.line1.set_data([], [])
            if root.ch2on.get():
                root.line2.set_data(t*1e6, sigch2plot*ch2vdiv + ch2vpos)
            else:
                root.line2.set_data([], [])
            if root.runloopsingle.get():  # Single sweep - run once and turn Run/Stop to false
                root.runloopsingle.set(False)
                root.runloopon.set(False)
            measuresig(root, sigch1buffer, sigch2buffer)

    updatecursors(root)
    try:
        root.canvas.draw()
        root.canvas.flush_events()
    except:
        print("drawing error")
        root.runloopon.set(False)

def adc(root, sig):
    rho = 0.99992
    # rho = 0.99995
    b = [1, -2, 1]
    a = [1, -2 * rho, rho * rho]


    sigdcb, root.zidcb = signal.lfilter(b, a, sig, zi=root.zidcb)
    sigaaf, root.ziaaf = signal.lfilter(root.aafiltcoefs, 1, sigdcb, zi=root.ziaaf)
    sigadc = np.minimum(np.maximum(np.round((sigaaf[range(0, root.CHUNK, root.OVERSMPLRATIO)]) * 2**17), -2**15), 2**15)
    return sigadc


def dac(root, sig):
    # sigdac = np.kron(sig, np.ones(root.OVERSMPLRATIO)) / 2**15
    # sigdac = sigdac[range(0, root.CHUNK)]
    sigdac = (np.kron(sig, np.ones(root.OVERSMPLRATIO)) / 2**17)[range(0, root.CHUNK)]
    sigdacinterp, root.ziinterp = signal.lfilter(root.aafiltcoefs, 1, sigdac, zi=root.ziinterp)
    return sigdacinterp


def dsp(root, sig):
    dspmode = root.dspmode.get()
    if dspmode == 1:  # bypass
        sigdsp = sig
    elif dspmode == 2:  # sign-flip modulation
        sigdsp = np.multiply(sig, np.power(-1, list(range(0, np.size(sig)))))
    elif dspmode == 3:  # NCO
        sigdsp = sig * 0
        for n in range(int(root.CHUNK / root.OVERSMPLRATIO)):
            sigdsp[n] = lab2nco.NCO(int((root.N0 - root.CHUNK) / root.OVERSMPLRATIO) + n, int(root.ncoK.get()), root.sin_array)

    return sigdsp

def loadfiltercoefs():
    coef = []
    infile = open("filtercoefs.txt", "r")
    for line in infile:
        coef.append(float(line)/2**18 * np.exp(1) / 3)
    infile.close()
    return coef