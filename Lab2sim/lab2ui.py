import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import lab2utils
import lab2nco


# main run loop
def runloop(root):
    while root.runloopon.get():
        # if root.runloopsingle.get():  # Single sweep - run once and turn Run/Stop to false
        #     root.runloopsingle.set(False)
        #     root.runloopon.set(False)
        root.sigin = lab2utils.siggen(root)  # Generate signal batch
        root.sigadc = lab2utils.adc(root, root.sigin)  # simulate sampling process
        root.sigdsp = lab2utils.dsp(root, root.sigadc)  # perform DSP block
        root.sigout = lab2utils.dac(root, root.sigdsp)  # simulate reconstruction process
        lab2utils.sigplot(root, root.sigin, root.sigout)  # plot output


class Lab2UI:
    SMPLFREQ = 48000
    CHUNK = 10000
    SMPLFREQALPHA = np.exp(1) * 3e-5
    OVERSMPLRATIO = 10
    N0 = 0
    n = []
    t = []
    sigin = []
    sigadc = []
    sigdsp = []
    sigout = []
    sigch1buffer = []
    sigch2buffer = []
    triggerlocked = False
    triggerwindow = 16
    sigtrigres = []
    sigtrigresvalid = False
    MAXMEAS = 4
    zidcb = [0, 0]

    def __init__(self):
        self.window = tk.Tk()
        self.window.title("Lab2")
        self.window.resizable(0, 0)
        self.runloopon = tk.BooleanVar(False)
        self.runloopsingle = tk.BooleanVar(False)
        self.t0 = tk.DoubleVar()

        # Input selector and config
        frame_sigsrc = tk.LabelFrame(text="Signal source")
        frame_sigsrc.grid(column=0, row=0, sticky=tk.W + tk.N, pady=10, padx=10)

        self.inputsource = tk.IntVar()
        self.inputsource.set(1)
        rad_srcsiggen = tk.Radiobutton(frame_sigsrc, text='Signal Generator', value=1, variable=self.inputsource)
        rad_srcsiggen.grid(column=0, row=1, sticky=tk.W)
        frame_sigsrc_siggen = tk.Frame(frame_sigsrc)
        frame_sigsrc_siggen.grid(column=0, row=2, sticky=tk.W, padx=20)
        lbl_freq = tk.Label(frame_sigsrc_siggen, text="Frequency (Hz):")
        lbl_freq.grid(column=0, row=0, sticky=tk.W)
        self.sigfreq = tk.DoubleVar()
        self.sigfreq.set(1000.)
        spn_freq = tk.Spinbox(frame_sigsrc_siggen, from_=0, to=50000, width=7, textvariable=self.sigfreq)
        spn_freq.grid(column=1, row=0, sticky=tk.W)
        lbl_amp = tk.Label(frame_sigsrc_siggen, text="Amplitude (mV):")
        lbl_amp.grid(column=0, row=1, sticky=tk.W)
        self.sigamp = tk.DoubleVar()
        self.sigamp.set(500.)
        spn_sigamp = tk.Spinbox(frame_sigsrc_siggen, from_=0, to=2000, width=5, textvariable=self.sigamp)
        spn_sigamp.grid(column=1, row=1, sticky=tk.W)
        self.sigshape = tk.IntVar()
        self.sigshape.set(1)
        frame_sigsrc_sigtype = tk.Frame(frame_sigsrc_siggen)
        frame_sigsrc_sigtype.grid(column=0, row=2, sticky=tk.W)
        rad_sigshapesin = tk.Radiobutton(frame_sigsrc_sigtype, text="sin", indicatoron=0,
                                         value=1, variable=self.sigshape)
        rad_sigshaperect = tk.Radiobutton(frame_sigsrc_sigtype, text="rect", indicatoron=0,
                                          value=2, variable=self.sigshape)
        rad_sigshapetri = tk.Radiobutton(frame_sigsrc_sigtype, text="tri", indicatoron=0,
                                         value=3, variable=self.sigshape)
        rad_sigshapesaw = tk.Radiobutton(frame_sigsrc_sigtype, text="saw", indicatoron=0,
                                         value=4, variable=self.sigshape)
        rad_sigshapesin.grid(column=0, row=0, sticky=tk.W)
        rad_sigshaperect.grid(column=1, row=0, sticky=tk.W)
        rad_sigshapetri.grid(column=2, row=0, sticky=tk.W)
        rad_sigshapesaw.grid(column=3, row=0, sticky=tk.W)
        self.pulsemod = tk.BooleanVar()
        self.pulsemod.set(False)
        chk_modpulse = tk.Checkbutton(frame_sigsrc_siggen, text="Pulse modulation", variable=self.pulsemod)
        chk_modpulse.grid(column=0, row=3, sticky=tk.W, columnspan=2)
        lbl_pulsepw = tk.Label(frame_sigsrc_siggen, text="Pulse width [usec]")
        lbl_pulsepw.grid(column=0 ,row=4, sticky=tk.W)
        self.pulsepw = tk.DoubleVar()
        self.pulsepw.set(1000.)
        spn_pulsepw = tk.Spinbox(frame_sigsrc_siggen, from_=0, to=50000, width=7, textvariable=self.pulsepw)
        spn_pulsepw.grid(column=1, row=4, sticky=tk.W)
        lbl_pulsepri = tk.Label(frame_sigsrc_siggen, text="Pulse period [usec]")
        lbl_pulsepri.grid(column=0 ,row=5, sticky=tk.W)
        self.pulsepri = tk.DoubleVar()
        self.pulsepri.set(10000.)
        spn_pulsepri = tk.Spinbox(frame_sigsrc_siggen, from_=0, to=100000, width=7, textvariable=self.pulsepri)
        spn_pulsepri.grid(column=1, row=5, sticky=tk.W)

        rad_srcmic = tk.Radiobutton(frame_sigsrc, text='Microphone', value=2, variable=self.inputsource, state="disabled")
        rad_srcmic.grid(column=0, row=3, sticky=tk.W)
        rad_srcfile = tk.Radiobutton(frame_sigsrc, text='Input from file', value=3, variable=self.inputsource, state="disabled")
        rad_srcfile.grid(column=0, row=4, sticky=tk.W)

        # Signal processing config
        frame_dsp = tk.LabelFrame(text="Signal processing")
        frame_dsp.grid(column=0, row=1, sticky=tk.W + tk.N, pady=10, padx=10)
        self.dspmode = tk.IntVar()
        self.dspmode.set(1)
        rad_bypass = tk.Radiobutton(frame_dsp, text='Bypass', value=1, variable=self.dspmode)
        rad_bypass.grid(column=0, row=0, columnspan=2, sticky=tk.W)
        rad_signflip = tk.Radiobutton(frame_dsp, text="sign-flip modulation", value=2, variable=self.dspmode)
        rad_signflip.grid(column=0, row=1, columnspan=2, sticky=tk.W)
        rad_nco = tk.Radiobutton(frame_dsp, text="NCO, K:", value=3, variable=self.dspmode)
        rad_nco.grid(column=0, row=3, sticky=tk.W)
        self.ncoK = tk.DoubleVar()
        self.ncoK.set(4)
        spn_ncoK = tk.Spinbox(frame_dsp, from_=0, to=100, width=2, textvariable=self.ncoK)
        spn_ncoK.grid(column=1, row=3, sticky=tk.W)
        self.sin_array = lab2nco.sin_array()

        # Output and display
        frame_out = tk.LabelFrame(text="Signal display")
        frame_out.grid(column=1, row=0, rowspan=2, sticky=tk.W+tk.N+tk.E+tk.S, pady=10, padx=10)
        frame_out_NWpnl = tk.Frame(frame_out)
        frame_out_NWpnl.grid(column=0, row=0)
        frame_out_Npnl = tk.Frame(frame_out)
        frame_out_Npnl.grid(column=1, row=0)
        frame_out_Wpnl = tk.Frame(frame_out)
        frame_out_Wpnl.grid(column=0, row=1, padx=10, pady=10)
        frame_out_tbase = tk.Frame(frame_out_Wpnl)
        frame_out_tbase.grid(column=0, row=0, pady=10)
        frame_out_trigger = tk.Frame(frame_out_Wpnl)
        frame_out_trigger.grid(column=0, row=1, pady = 10)
        frame_out_meascrsr = tk.Frame(frame_out_Wpnl)
        frame_out_meascrsr.grid(column=0, row=2, pady=10)
        frame_out_plot = tk.Frame(frame_out)
        frame_out_plot.grid(column=1, row=1, padx=10, pady=10, sticky=tk.N+tk.W+tk.E+tk.S)
        frame_out_SWpnl = tk.Frame(frame_out)
        frame_out_SWpnl.grid(column=0, row=2)
        frame_out_Spnl = tk.Frame(frame_out)
        frame_out_Spnl.grid(column=1, row=2, padx=10, pady=10)
        frame_out_ch1 = tk.Frame(frame_out_Spnl)
        frame_out_ch1.grid(column=0, row=0, padx=10)
        frame_out_ch2 = tk.Frame(frame_out_Spnl)
        frame_out_ch2.grid(column=1, row=0, padx=10)

        if False:
            self.fig, self.ax = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(5, 4))
        else:
            self.fig = plt.figure(figsize=(5, 4))
            self.ax = plt.subplot2grid((3,1), (0,0), rowspan=2)
            plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.05, hspace=0.1)
            self.axtxt = plt.subplot2grid((5, 1), (4, 0))
            self.axtxt._frameon = False
            self.axtxt.set_xticks([])
            self.axtxt.set_yticks([])
            self.meastext = self.axtxt.text(0, 0, "", fontsize=8)
            self.cursortext = self.axtxt.text(0.5, 0, "", fontsize=8)

        self.ax.set_ylim([-2, 2])
        self.ax.set_yticklabels([])
        self.ax.grid()
        self.line1, = self.ax.plot([], [])
        self.line2, = self.ax.plot([], [])
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame_out_plot)  # A tk.DrawingArea.
        # self.canvas.get_tk_widget().grid(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.draw()
        # ch1&ch2 on/off control
        self.ch1on = tk.BooleanVar()
        self.ch1on.set(True)
        self.ch2on = tk.BooleanVar()
        self.ch2on.set(True)

        def ch1show():
            self.ch1on.set(not (self.ch1on.get()))

        def ch2show():
            self.ch2on.set(not (self.ch2on.get()))

        btn_ch1show = tk.Button(frame_out_ch1, text="Ch1 on", command=ch1show)
        btn_ch1show.grid(column=0, row=0, rowspan=2)
        btn_ch2show = tk.Button(frame_out_ch2, text="Ch2 on", command=ch2show)
        btn_ch2show.grid(column=0, row=0, rowspan=2)
        # ch1&ch2 V/Div & v-position control
        self.chvdivVec = np.reshape((np.array([[1], [2], [5]]) * np.power(10, range(4))).transpose(), 3 * 4) / 100
        self.ch1vdivIdx = tk.IntVar()
        self.ch1vdivIdx.set(7)
        self.ch2vdivIdx = tk.IntVar()
        self.ch2vdivIdx.set(7)
        self.ch1vdiv = tk.DoubleVar()
        self.ch2vdiv = tk.DoubleVar()
        self.ch1vpos = tk.DoubleVar()
        self.ch1vpos.set(1)
        self.ch2vpos = tk.DoubleVar()
        self.ch2vpos.set(-1)

        def ch1vdivinc():
            self.ch1vdivIdx.set(min(self.chvdivVec.__len__() - 1, self.ch1vdivIdx.get() + 1))
            vdivupdate()

        def ch1vdivdec():
            self.ch1vdivIdx.set(max(0, self.ch1vdivIdx.get() - 1))
            vdivupdate()

        def ch2vdivinc():
            self.ch2vdivIdx.set(min(self.chvdivVec.__len__() - 1, self.ch2vdivIdx.get() + 1))
            vdivupdate()

        def ch2vdivdec():
            self.ch2vdivIdx.set(max(0, self.ch2vdivIdx.get() - 1))
            vdivupdate()

        def vdivupdate():
            self.ch1vdiv.set(self.chvdivVec[self.ch1vdivIdx.get()])
            self.ch2vdiv.set(self.chvdivVec[self.ch2vdivIdx.get()])

        def ch1vposup():
            self.ch1vpos.set(self.ch1vpos.get() + 0.1)

        def ch1vposdn():
            self.ch1vpos.set(self.ch1vpos.get() - 0.1)

        def ch2vposup():
            self.ch2vpos.set(self.ch2vpos.get() + 0.1)

        def ch2vposdn():
            self.ch2vpos.set(self.ch2vpos.get() - 0.1)

        btn_ch1vdivinc = tk.Button(frame_out_ch1, text="V/Div+", command=ch1vdivinc)
        btn_ch1vdivinc.grid(column=1, row=0)
        btn_ch1vdivdec = tk.Button(frame_out_ch1, text="V/Div-", command=ch1vdivdec)
        btn_ch1vdivdec.grid(column=1, row=1)
        btn_ch1vposup = tk.Button(frame_out_ch1, text="v-pos+", command=ch1vposup)
        btn_ch1vposup.grid(column=2, row=0)
        btn_ch1vposdn = tk.Button(frame_out_ch1, text="v-pos-", command=ch1vposdn)
        btn_ch1vposdn.grid(column=2, row=1)
        btn_ch2vdivinc = tk.Button(frame_out_ch2, text="V/Div+", command=ch2vdivinc)
        btn_ch2vdivinc.grid(column=1, row=0)
        btn_ch2vdivdec = tk.Button(frame_out_ch2, text="V/Div-", command=ch2vdivdec)
        btn_ch2vdivdec.grid(column=1, row=1)
        btn_ch2vposup = tk.Button(frame_out_ch2, text="v-pos+", command=ch2vposup)
        btn_ch2vposup.grid(column=2, row=0)
        btn_ch2vposdn = tk.Button(frame_out_ch2, text="v-pos-", command=ch2vposdn)
        btn_ch2vposdn.grid(column=2, row=1)
        vdivupdate()

        # T/Div control
        def tbaseinc():
            self.tbaseIdx.set(min(self.tbaseVec.__len__() - 1, self.tbaseIdx.get() + 1))
            tbaseupdate()

        def tbasedec():
            self.tbaseIdx.set(max(0, self.tbaseIdx.get() - 1))
            tbaseupdate()

        def tbaseupdate():
            self.tbase.set(self.tbaseVec[self.tbaseIdx.get()])
            self.ax.set_xlim(0, self.tbase.get())
            xticks = np.linspace(0, self.tbase.get(), 11)
            self.ax.set_xticks(xticks)
            tbasemag = np.floor(np.log10(self.tbase.get() / 10) / 3).astype(int)
            xticksnorm = (xticks / np.power(1000, tbasemag)).astype(int)
            self.ax.set_xticklabels(xticksnorm)
            tbasemagstr = ["usec", "msec", "sec"]
            self.ax.set_xlabel(tbasemagstr[tbasemag])
            try:
                for c in self.cursors:
                    c.updatetbase()
            except:
                print()

        self.tbaseIdx = tk.IntVar()
        self.tbaseIdx.set(9)
        self.tbaseVec = 10 * np.reshape((np.array([[1], [2], [5]]) * np.power(10, range(8))).transpose(), 3 * 8)
        self.tbase = tk.DoubleVar()
        tbaseupdate()
        btn_tbaseinc = tk.Button(frame_out_tbase, text="T/Div+", command=tbaseinc)
        btn_tbaseinc.grid(column=0, row=0)
        btn_tbasedec = tk.Button(frame_out_tbase, text="T/Div-", command=tbasedec)
        btn_tbasedec.grid(column=1, row=0)
        # trigger control
        chlist = ["ch1",
                  "ch2"]
        self.triggerchslctd = tk.StringVar()
        self.triggerchslctd.set((chlist[0]))
        lbl_trigger = tk.Label(frame_out_trigger, text="trigger:")
        lbl_trigger.grid(column=0, row=0, sticky=tk.W)
        lst_chtrigger = tk.OptionMenu(frame_out_trigger, self.triggerchslctd, *chlist)
        lst_chtrigger.grid(column=1, row=0, sticky=tk.W)
        # measure & cursurs control
        measchslctd = tk.StringVar()
        measchslctd.set((chlist[0]))
        measlist = ["pk-pk amp",
                    "RMS amp",
                    "max amp",
                    "min amp",
                    "freq"]
        measslctd = tk.StringVar()
        measslctd.set(measlist[0])
        cursorstatelist = ["off",
                           "cursor 1",
                           "cursor 2",
                           "both"]
        self.xcursormode = 0
        self.ycursormode = 0
        self.meascount = 0
        self.measch = []
        self.meastype = []
        self.measdata = []
        self.measalpha = 0.3
        def addmeas():
            self.meascount = self.meascount + 1
            self.measch = np.concatenate(([measchslctd.get()], self.measch))
            self.meastype = np.concatenate(([measslctd.get()], self.meastype))
            self.measdata = np.concatenate(([np.nan], self.measdata))
            if self.meascount > self.MAXMEAS:
                self.meascount = self.MAXMEAS
                self.measch = self.measch[:-1]
                self.meastype = self.meastype[:-1]
                self.measdata = self.measdata[:-1]
        def xcursor():
            self.xcursormode = (self.xcursormode + 1)%4
            lbl_xcursor.config(text=cursorstatelist[self.xcursormode])
            self.cursorx1.visible((self.xcursormode % 2) == 1)  # ch1 on or both on
            self.cursorx2.visible((self.xcursormode >> 1) == 1)  # ch2 on or both on

        def ycursor():
            self.ycursormode = (self.ycursormode + 1)%4
            lbl_ycursor.config(text=cursorstatelist[self.ycursormode])
            self.cursory1.visible((self.ycursormode % 2) == 1)  # ch1 on or both on
            self.cursory2.visible((self.ycursormode >> 1) == 1)  # ch2 on or both on

        class cursor:
            def __init__(self, cursortype, line, root):
                self.line = line[0]
                self.cursortype = cursortype
                self.press = None
                self.visible(False)
                self.cidpress = self.line.figure.canvas.mpl_connect('button_press_event', self.on_press)
                self.cidmotion = self.line.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
                self.cidrelease = self.line.figure.canvas.mpl_connect('button_release_event', self.on_release)
                self.root = root
            def visible(self, visibility):
                self.line.set_visible(visibility)
            def isvisible(self):
                return self.line.get_visible()
            def updatetbase(self):
                if self.cursortype == 'x':
                    if self.line.axes.get_xlim()[1] < self.line.get_xdata()[0]:
                        self.line.set_xdata(self.line.axes.get_xlim()[1] * np.ones([2]) * 0.7)

            def getvalue(self):
                if self.cursortype == 'x':
                    val = self.line.get_xdata()[0]
                elif self.cursortype == 'y':
                    val = self.line.get_ydata()[0]
                return val

            def on_press(self, event):
                'on button press we will see if the mouse is over us and store some data'
                if event.inaxes != self.line.axes: return
                contains, attrd = self.line.contains(event)
                if not contains: return
                if not self.isvisible(): return
                x0 = self.line.get_xdata()[0]
                y0 = self.line.get_ydata()[0]
                self.press = x0, y0, event.xdata, event.ydata

            def on_motion(self, event):
                'on motion we will move the rect if the mouse is over us'
                if self.press is None: return
                if event.inaxes != self.line.axes: return
                x0, y0, xpress, ypress = self.press
                dx = event.xdata - xpress
                dy = event.ydata - ypress
                #      (x0, xpress, event.xdata, dx, x0+dx))
                if self.cursortype == "x":
                    self.line.set_xdata([x0 + dx, x0 + dx])
                elif self.cursortype == "y":
                    self.line.set_ydata([y0 + dy, y0 + dy])
                self.line.figure.canvas.draw()

            def on_release(self, event):
                'on release we reset the press data'
                self.press = None
                self.line.figure.canvas.draw()
                lab2utils.sigplot(self.root, self.root.sigin, self.root.sigout)

        self.cursorx1 = cursor('x', self.ax.plot(np.diff(np.array(self.ax.get_xlim())) * [0.3, 0.3] + self.ax.get_xlim()[0], np.array(self.ax.get_ylim()), 'g--'), root=self)
        self.cursorx2 = cursor('x', self.ax.plot(np.diff(np.array(self.ax.get_xlim())) * [0.7, 0.7] + self.ax.get_xlim()[0], np.array(self.ax.get_ylim()), 'g--'), root=self)
        self.cursory1 = cursor('y', self.ax.plot(np.array(self.ax.get_xlim()), np.diff(np.array(self.ax.get_ylim())) * [0.25, 0.25] + self.ax.get_ylim()[0], 'g--'), root=self)
        self.cursory2 = cursor('y', self.ax.plot(np.array(self.ax.get_xlim()), np.diff(np.array(self.ax.get_ylim())) * [0.75, 0.75] + self.ax.get_ylim()[0], 'g--'), root=self)
        self.cursors = [self.cursorx1, self.cursorx2, self.cursory1, self.cursory2]

        lbl_meas = tk.Label(frame_out_meascrsr, text="meas & cursor:")
        lbl_meas.grid(column=0, row=0, sticky=tk.W)
        lst_chmeas = tk.OptionMenu(frame_out_meascrsr, measchslctd, *chlist)
        lst_chmeas.grid(column=1, row=0, sticky=tk.W, pady=5)
        lst_measlist = tk.OptionMenu(frame_out_meascrsr, measslctd, *measlist)
        lst_measlist.grid(column=0, row=1, sticky=tk.W)
        btn_measadd = tk.Button(frame_out_meascrsr, text="Add", command=addmeas)
        btn_measadd.grid(column=1, row=1, sticky=tk.W)
        btn_xcursor = tk.Button(frame_out_meascrsr, text="x cursor", command=xcursor)
        btn_xcursor.grid(column=0, row=2, sticky=tk.W)
        lbl_xcursor = tk.Label(frame_out_meascrsr, text=cursorstatelist[0], width=10, anchor=tk.W)
        lbl_xcursor.grid(column=1, row=2, sticky=tk.W, padx=5)
        btn_ycursor = tk.Button(frame_out_meascrsr, text="y cursor", command=ycursor)
        btn_ycursor.grid(column=0, row=3, sticky=tk.W)
        lbl_ycursor = tk.Label(frame_out_meascrsr, text=cursorstatelist[0], width=10, anchor=tk.W)
        lbl_ycursor.grid(column=1, row=3, sticky=tk.W, padx=5)

        self.aafiltcoefs = lab2utils.loadfiltercoefs()
        self.ziaaf = np.zeros(np.size(self.aafiltcoefs)-1)
        self.ziinterp = np.zeros(np.size(self.aafiltcoefs)-1)

        # Run/Stop & Single sweep control
        def runstop():
            self.runloopon.set(not (self.runloopon.get()))
            if self.runloopon.get():
                runloop(self)
        def single():
            self.runloopsingle.set(True)
            self.runloopon.set(True)
            runloop(self)
        btn_runstop = tk.Button(frame_out_Npnl, text="Run/Stop", command=runstop)
        btn_runstop.grid(column=0, row=0, sticky=tk.E)
        btn_single = tk.Button(frame_out_Npnl, text="Single", command=single)
        btn_single.grid(column=1, row=0, sticky=tk.E)

        # save figure
        def savefigure():
            figfile = tk.filedialog.asksaveasfile(title = "Select file",
                                                      filetypes = (("png files","*.png"), ("jpeg files","*.jpg"), ("all files","*.*")))
            if figfile != None:
                plt.savefig(figfile.name)
        btn_savefigure = tk.Button(frame_out_Npnl, text="Save figure", command=savefigure)
        btn_savefigure.grid(column=2, row=0)


        self.canvas.draw()
        self.canvas.flush_events()
