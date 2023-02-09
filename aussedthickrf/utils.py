import numpy as np
import matplotlib
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt


def signoise(self):  # , winsig, winnoise, relative='onset'):
    """
    Determine signal noise ratio by dividing the maximum in the two windows.
    """
    st = self.stats
    self_copy = self.copy()
    self_copy.detrend().taper(max_percentage=0.05)
    self_copy.filter("bandpass", freqmin=0.1, freqmax=1)  # ,corners=2, zerophase=True)
    winsig = [-5, 25]  # signal window
    winnoise = [-45, -15]  # noise window
    rel_time = getattr(st, "onset")
    winsig0 = rel_time - st.starttime + winsig[0]
    winsig1 = rel_time - st.starttime + winsig[1]
    winnoise0 = rel_time - st.starttime + winnoise[0]
    winnoise1 = rel_time - st.starttime + winnoise[1]
    #
    t = np.arange(self.stats.npts) * 1.0 / st.sampling_rate
    datasig = self_copy.data[(t >= winsig0) * (t <= winsig1)]
    datanoise = self_copy.data[(t >= winnoise0) * (t <= winnoise1)]

    try:
        st.signoise = max(abs(datasig)) / max(abs(datanoise))
        return st.signoise
    except:
        st.signoise = 0
        return st.signoise


def calc_cc_conv(data_rf, data):
    """
    Returns rf stream with correlation coeffecient in tr.stats.cc_conv
    """
    rf_all = data_rf
    sig_all = data.copy()
    rf_all.sort(keys=["onset"])
    sig_all.sort(keys=["onset"])
    sig_all.rotate("NE->RT")
    sig_all.trim2(-25, 75, "onset")
    ###
    cc = []
    fit = []
    if len(rf_all) == len(sig_all):
        # print('True')
        for i in range(0, int(len(rf_all) / 3)):
            # for i in range(0,2)):
            try:
                rf = rf_all[3 * i : 3 * i + 3]
                sig = sig_all[3 * i : 3 * i + 3]
                # rf,sig are stream containing same three traces for RF and ZNE
                obs_Z = sig.select(component="Z")[0].copy()
                obs_R = sig.select(component="R")[0].copy()
                obs_rfR = rf.select(component="R")[0].copy()
                sr = obs_Z.stats.sampling_rate
                # Filter using SNR bandpass
                obs_Z.detrend().taper(max_percentage=0.05, max_length=2.0)
                obs_R.detrend().taper(max_percentage=0.05, max_length=2.0)
                obs_Z.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)
                obs_R.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)
                obs_rfR.filter("bandpass", freqmin=0.1, freqmax=1.0, zerophase=True)

                pred_R = obs_R.copy()
                pred_R.stats.channel = "PRR"
                st = pred_R.stats.starttime
                time_shift = pred_R.stats.onset - pred_R.stats.starttime  # rel P onset
                ind1 = int(np.ceil(time_shift * sr))

                ind2 = ind1 + len(obs_Z.data)  # [leadin:leadin + n]
                pred_R.data = np.convolve(obs_Z.data, obs_rfR.data, mode="full")[
                    ind1:ind2
                ]

                obs_Z.trim(st + 20, st + 45)
                obs_R.trim(st + 20, st + 45)
                pred_R.trim(st + 20, st + 45)
                rf.select(component="R")[0].stats.cc_conv = np.corrcoef(
                    obs_R.data * obs_R.data, pred_R.data * pred_R.data
                )[0][1]
                # rf.select(component='R')[0].stats.cc_conv = np.corrcoef(abs(obs_R.data), abs(pred_R.data))[0][1]
                cc.append(np.corrcoef(obs_R.data, pred_R.data)[0][1])
            except:
                print("  ")  # do nothing
    return cc, rf_all


def max_p_onset(self):
    """
    Determine max amplitude of r-rf around P-onset. winP_st and winP_end should be adjusted in case of very thick sediments i.e. > 2km.
    """
    onset = self.stats.onset
    self_copy = self.copy()
    self_copy.trim(onset - 5, onset + 30)
    self_copy.taper(max_percentage=0.05)
    st = self_copy.stats
    rel_time = getattr(st, "onset")

    winP_st = rel_time - st.starttime - 0.1  # P onset -1
    winP_end = rel_time - st.starttime + 1.5  # P onset +1
    t = np.arange(st.npts) * 1.0 / st.sampling_rate
    max_abs = max(
        abs(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    )  # gets max amp around P (-.1,1.5 sec)
    max_ = max(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    self.stats.max_abs_P = max_abs
    self.stats.max_P = max_
    return max_abs, max_


def plot_stacks(stream, fig_width=7., trace_height=0.5,
                stack_height=0.5, dpi=None,
                scale=1, trim=None,
                show_vlines=False):
    """
    Plot receiver functions.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param dpi: dots per inch for the created figure
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    :param show_vlines: If True, show vertical alignment grid lines on plot
        at positions of the major x-tick marks.
    """

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    stats = stream[0].stats
    times = stream[0].times() - (stats.onset - stats.starttime)
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH), dpi=dpi)
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])

    def _plot(ax, t, d, i, color, label):
        c1, c2 = (color, 'k')
        ax.text(-3, i + 0.2, label)
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k')

    max_ = np.array([np.max(np.abs(tr.data)) for tr in stream])

    delays = np.array([tr.stats.delay for tr in stream])
    norm = matplotlib.colors.Normalize(vmin=np.min(delays), 
                                       vmax=np.max(delays))

    for i, tr in enumerate(stream):
        rgba_color = cm.turbo(tr.stats.delay)
        _plot(ax1, times, tr.data / max_[i] * scale, i + 1,
              rgba_color,
              f"{tr.stats.network}.{tr.stats.station}")

    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    sm = plt.cm.ScalarMappable(cmap='turbo', norm=norm)

    cbar = plt.colorbar(sm, orientation='horizontal')
    cbar.set_label('Delay [s]')
    return fig
