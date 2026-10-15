import numpy as np
import pandas as pd
import re
import TraceIO as tio
import matplotlib.pyplot as plt
from lmfit import Model
import glob

# Global constants
ticklength_s = 1e-7  # Timebase (s)
routeindex = 'GR'  # order of APDs


def T3R_to_h5(filename, report=True):
    trace = tio.Traces(filename)

    print(f'FCS>>> Processing file: {filename}')
    # Read T3R file
    headersize = 160
    with open(filename, 'rb') as file:
        NumberOfRecords = int(np.fromfile(file, np.int32, headersize)[0])
        records = np.fromfile(file, np.uint32, NumberOfRecords)
    TimeTag = records & 65535  # the lowest 16 bits
    Channel = (records >> 16) & 4095  # the next 12 bits
    Route = ((records >> 28) & 3)  # the next 2 bits
    Valid = (records >> 30) & 1  # the next bit

    # Put data in Pandas dataframe
    data = {"TimeTag": TimeTag, "Channel": Channel, "Route": Route, "Valid": Valid}
    timetags = pd.DataFrame(data)
    timetags['Route'] = timetags['Route'].astype(np.uint8)

    timetags['overflow'] = (timetags['Channel'] == 2048) & (
            timetags['Valid'] == 0)  # If Valid=0 and Channel=bit(2048), overflow=True
    timetags['n_flows'] = timetags['overflow'].cumsum().astype('int64')  # Cummulative sum of overflow = True
    timetags['ArrTime'] = timetags['TimeTag'] + (timetags['n_flows'] * (2 ** 16))  # arrival time in 100ns resolution
    timetags.set_index('ArrTime', inplace=True)  # Reset (indices) of dataframe
    timetags['Det'] = timetags['Route'].map({1: 'R', 0: 'G'})  # Wrong offset when this line is deleted
    timetags.dropna(inplace=True)

    # Clean up dataframe
    timetags = timetags[timetags['Valid'] == 1]  # Delete all non-arrival-time records
    timetags.drop(['overflow', 'n_flows', 'Valid', 'TimeTag', 'Channel'], axis=1, inplace=True)
    timetags.index += timetags.index.min()  # Set first arrival time to zero

    # Align with excitation scheme
    pulses = get_pulses(filename)
    pulseperiod_ticks = pulses['length (ticks)'].sum()
    timetags['Exc_index'] = np.mod(timetags.index, pulseperiod_ticks).astype('int64')
    timetags_PIE = timetags.groupby(['Exc_index', 'Route']).size().unstack()
    timetags.index -= get_offset(timetags_PIE, pulses)
    timetags.drop(['Exc_index', 'Det'], axis=1, inplace=True)

    # Calculate photon statistics
    for route in routeindex:
        for phase in pulses.index.tolist():
            intensity = select_timetags(timetags, pulses, route=route, phase=phase).count()[0]
            intensity /= timetags.index.max() * ticklength_s
            intensity /= 1000
            pulses.loc[phase, f'{route}: I (kHz)'] = intensity
            pulses.loc[phase, f'{route}: bg (kHz)'] = None
    pulses = get_background(filename, timetags, pulses, plot=False)

    timetags.to_hdf(trace.filename(ext='h5'), '_timetags', mode='a', format='table')
    pulses.to_hdf(trace.filename(ext='h5'), '_pulsescheme', mode='a', format='table')

    # Save signals characteristics
    for label in ['G514', 'R514', 'G632', 'R632']:
        I = len(select_timetags(timetags, pulses, route=label[0], wavelength=label[1:]))
        if I > 0:
            I /= timetags.index[-1] * ticklength_s * 1000.
            trace.write_par(f't (s)', float(timetags.index[-1] * ticklength_s), label=label)
            trace.write_par(f'Iexc (%)', pulses[pulses['wavelength (nm)'] == int(label[1:])].iat[0, 1], label=label)
            trace.write_par(f'I (kHz)', float(I), label=label)
            background = float(pulses[pulses['wavelength (nm)'] == int(label[1:])][f'{label[0]}: bg (kHz)'])
            trace.write_par(f'bg (kHz)', background, label=label)
    return pulses


def get_pulses(filename):
    try:
        trace = tio.Traces(filename)
        pulses = pd.read_hdf(trace.filename(ext='h5'), '_pulsescheme')
    except:
        # Extract from filename the pulse_sequence, pulse_period and resolution in 100ns
        parts = re.split("ns|us|ms|cont|_", filename)
        if filename.find('cont') == -1:
            PIE_time = ['ns', 'us', 'ms']
            for time_res in PIE_time:
                if time_res in filename:
                    pulse_time__s = int(filename[filename.find('ALEX') + 4:filename.find(time_res)]) \
                                    * 1e-9 * 10 ** (PIE_time.index(time_res) * 3)
                    pulse_scheme = list(parts[-2])
                    pulse_seq = []  # pulse_seq is pulse_scheme as a list of excitation
                    for phase in pulse_scheme:  # wavelengths instead of letters
                        if phase == 'R':
                            pulse_seq.append(632)
                        if phase == 'G':
                            pulse_seq.append(514)
                        #    if pulse_let == 'C':
                        #        pulse_seq.append(488)
                        #    if pulse_let == 'B':
                        #        pulse_seq.append(404)
                        if phase == 'd':
                            pulse_seq.append(0)

                    pulse_length__ticks = []
                    for seq in pulse_seq:
                        pulse_length__ticks.append(int(pulse_time__s / ticklength_s))
        else:
            if parts[-2][0] == 'G':
                pulse_seq = [514]
            elif parts[-2][0] == 'R':
                pulse_seq = [632]
            pulse_length__ticks = [1]

        pulse_int = np.array(re.split('R|G|d', parts[-2].replace('d', 'd0'))[1:], dtype=np.int)

        pulses = pd.DataFrame(pulse_seq, columns=['wavelength (nm)'])
        pulses.index.name = 'phase'
        pulses['intensity (%)'] = pulse_int
        pulses['length (ticks)'] = pulse_length__ticks
    return pulses


def get_offset(timetags, pulses):
    if (pulses['wavelength (nm)'] == 514).any():
        wavelength = 514
        align_route = 0
    else:
        wavelength = 632
        align_route = 1

    template = []
    for _, pulse in pulses.iterrows():
        if (pulse['wavelength (nm)'] == wavelength):
            template.extend(np.ones(int(pulse['length (ticks)'])))
        else:
            template.extend(np.zeros(int(pulse['length (ticks)'])))

    align_score = []
    for i, _ in enumerate(template):
        align_score.append(np.sum(np.roll(template, -i) * timetags[align_route]))
    offset_ticks = len(template) - np.argmax(align_score)

    return offset_ticks


def read_timetags(filename, range_s=None, **kwargs):
    filename = filename.rsplit('.')[0] + '.h5'
    print(f'FCS>>> Opening {filename}')
    if range_s is None:
        timetags = pd.read_hdf(filename, '_timetags', mode='a', **kwargs)
    else:
        selection = f'index > {int(range_s[0]/ticklength_s)} &index < {int(range_s[1]/ticklength_s)}'
        timetags = pd.read_hdf(filename, '_timetags/data', mode='a', where=selection, **kwargs)
    return timetags


def select_timetags(timetags, pulses, route=None, phase=None, wavelength=None, time_s=None, non_dark=None):
    unselected = len(timetags)
    if route is not None:
        timetags = timetags[(timetags['Route'] == routeindex.index(route))]
    if phase is not None:
        timetags = timetags[timetags.index % len(pulses) == phase]
    if wavelength is not None:
        selected = []
        for i, pulse in pulses.iterrows():
            if (pulse['wavelength (nm)'] == int(wavelength)):
                selected.extend(i * pulse['length (ticks)'] + np.arange(pulse['length (ticks)']))
        timetags = timetags[np.isin(timetags.index % pulses['length (ticks)'].sum(), selected)]
    if time_s is not None:
        timetags = timetags[
            timetags.index > np.min(time_s) * ticklength_s & timetags.index < np.max(time_s) * ticklength_s]
    if non_dark is True:
        Exc = pulses[pulses['wavelength (nm)'] > 0].index.tolist()
        timetags = timetags[np.isin((timetags.index % len(pulses)), Exc)]
    # print(f'route: {route}, phase: {phase}, wavelength: {wavelength}, range (s) = {time_s}, selected: {len(timetags)} ({len(timetags)/unselected:g})')
    return timetags.copy()


def correlate(filename, labels):
    from pycorrelate import pcorrelate, make_loglags

    trace = tio.Traces(filename=filename)
    pulses = get_pulses(trace.filename())
    timetags = read_timetags(trace.filename())

    tau_s = make_loglags(-6, 0, 8)

    for label in labels:
        print(f'FCS>>> Correlating {label}')
        signals = label.split('x')
        try:
            s1 = select_timetags(timetags, pulses, route=signals[0][0], wavelength=signals[0][1:]).index.to_numpy()
            s2 = select_timetags(timetags, pulses, route=signals[1][0], wavelength=signals[1][1:]).index.to_numpy()
            G = pcorrelate(s1 * ticklength_s, s2 * ticklength_s, tau_s, normalize=True)
            trace.write_channel(data=G - 1, channel='G-1 (-)', group='correlation', label=label)
        except Exception as e:
            print(f'FCS>>> !Warning, missing signals for {label} in {filename}')

    tau_s = tau_s[1:] + tau_s[:-1] / 2  # take centers of bin
    trace.write_channel(data=tau_s, channel='tau (s)', group='correlation', label='')


def FitFun_correlation_single(tau, N=1, td_s=1e-3, w_um=0.5, P3=0.1, t3_s=1e-5, Pa=0.1, ta_s=1e-6):
    result = (N * ((1 + (tau / td_s)) * np.sqrt(1 + (w_um ** -2) * (tau / td_s)))) ** -1 \
             * (1 + (P3 / (1 - P3)) * np.exp(-tau / t3_s)) \
             * (1 + (Pa / (1 - Pa)) * np.exp(-tau / ta_s))
    return result


def fit_correlation(filename, labels):
    traces = tio.Traces(filename=filename)
    for label in labels:
        try:
            tau_s = traces.read_channel(channel='tau (s)', group='correlation', label='')
            G_1 = traces.read_channel(channel='G-1 (-)', group='correlation', label=label)

            # fit correlation model
            model = Model(FitFun_correlation_single)
            params = model.make_params()
            for p in params:
                params[p].set(vary=False)

            N = 1 / np.mean(G_1[(tau_s > 1e-4) & (tau_s < 1e-3)])  # initial estimate
            params['N'].set(N, vary=True, min=0)
            params['w_um'].set(0.5, vary=True, min=0)
            params['td_s'].set(1e-3, vary=True, min=0)

            result = model.fit(G_1, tau=tau_s, params=params)
            fit = model.eval(result.params, tau=tau_s)

            traces.write_fit('correlation_single', 'G-1 (-)', fit, result.params, label=label)
            R2 = 1 - result.residual.var() / np.var(G_1)
            R2, _ = tio.round_significance(R2, 0.01)
            traces.write_par('R2_single', R2, label=label)
        except:
            print(f'FCS>>> !Warning fitting {label} in file {traces.filename()}')


def get_color(label):
    if label == 'R632':
        return 'red'
    elif label == 'R632xR632':
        return 'red'
    elif label == 'G514':
        return 'green'
    elif label == 'G514xG514':
        return 'green'
    elif label == 'R514':
        return 'blue'
    elif label == 'R514xR514':
        return 'blue'
    elif (label == 'G514xR632') | (label == 'R632xG514'):
        return 'orange'
    elif (label == 'R514xR632') | (label == 'R632xR514'):
        return 'darkviolet'
    elif (label == 'R514xG514') | (label == 'G514xR514'):
        return 'saddlebrown'
    else:
        return 'grey'


def plot_correlation(filename, labels, save=None, yrange=None):
    traces = tio.Traces(filename=filename)

    if (len(labels) == 1):
        pars = traces.parameters(label=labels[0])
    else:
        pars = None

    G0 = []
    for label in labels:
        tau_s = traces.read_channel(channel='tau (s)', group='correlation', label='')
        try:
            G_1 = traces.read_channel(channel='G-1 (-)', group='correlation', label=label)
            plt.scatter(tau_s, G_1, marker='o', color=get_color(label), facecolor='none')
            G0.append(np.median(G_1[tau_s < 5e-4]))
        except:
            pass
        try:
            fit = traces.read_channel(channel='G-1 (-)', group='fit_correlation_single', label=label)
            plt.plot(tau_s, fit, color=get_color(label))
        except:
            pass
        else:
            pars = None

    if yrange is None:
        yrange = [0, np.max(G0) * 1.5]

    plt.ylim(yrange)

    if save is None:
        filename = None
    else:
        filename = traces.filename(ext='emf', addition='correlation')
    tio.format_plot('tau (s)', 'G-1 (-)', pars=pars, xrange=[1e-5, 1], title=traces.plotname(),
                    scale_page=0.5, aspect=0.5, xlog=True, save=filename)


def plot_timetrace(filename, range_s=[0, 1], step_s=0.001, save=None, dark=False, Imax_kHz=None):
    import subprocess

    trace = tio.Traces(filename=filename)

    if range_s is None:
        timetags = read_timetags(trace.filename())
        range_s = [0, timetags.index.max() * ticklength_s]
        step_s = 0.01
    else:
        timetags = read_timetags(trace.filename(), range=range_s)

    time_s = np.linspace(range_s[0], range_s[1], int((range_s[1] - range_s[0]) / step_s))

    pulses = pd.read_hdf(trace.filename(), '_pulsescheme', mode='a')

    Imed = []
    for route in routeindex:
        for i, wavelength in pulses['wavelength (nm)'].iteritems():
            photons = select_timetags(timetags, pulses, route=route, wavelength=wavelength).index
            I, _ = np.histogram(photons, time_s / ticklength_s)
            if route == 'R': I *= -1
            if not dark and wavelength == 0:
                pass
            else:
                plt.plot(time_s[:-1], 0.001 * I / step_s, color=get_color(route + str(wavelength)))
                Imed.append(np.percentile(0.001 * np.abs(I) / step_s, 95))

    if Imax_kHz is None:
        Imax_kHz = np.max(Imed) * 2.0

    yrange = [-Imax_kHz, Imax_kHz]

    try:
        timetags = timetags[timetags['Burst'] > 0]
        I, _ = np.histogram(timetags.index, time_s / ticklength_s)
        I = I == 0
        I = np.where(I == 1, np.NaN, I)
        plt.plot(time_s[:-1], I, color='k', linewidth=3)
    except:
        pass

    if save is not None:
        plotfilename = trace.filename(ext=save, addition='timetrace')
    else:
        plotfilename = None

    tio.format_plot('time (s)', 'I (kHz)', title=trace.plotname(), scale_page=1.0, aspect=0.5, xrange=range_s,
                    yrange=yrange, save=plotfilename)


def exp(t, a0, f0):
    return a0 * np.exp(-t * f0)


def get_background(filename, timetags=None, pulses=None, treshold_ms=2.5, plot=False):
    from scipy.optimize import curve_fit
    trace = tio.Traces(filename)
    if timetags is None:
        timetags = read_timetags(trace.filename())
    if pulses is None:
        pulses = get_pulses(trace.filename())

    for j, phase in pulses.iterrows():
        # Make lagtime bins in multples of pulse period
        pulseperiod = pulses['length (ticks)'].sum()
        binsize = pulseperiod * (1000 // pulseperiod)
        bins = np.arange(0, 5e-3 / ticklength_s, binsize)
        for route in routeindex:
            label = route + str(int(phase['wavelength (nm)']))
            selectedTimetags = select_timetags(timetags, pulses, wavelength=label[1:], route=label[0])
            Intensity, bins = np.histogram(np.diff(selectedTimetags.index), bins)
            tau = bins[1:] * ticklength_s * 1000
            tau = tau[Intensity != 0]  # discard empty bins
            Intensity = Intensity[Intensity != 0]
            if plot: plt.scatter(tau, Intensity, facecolors='none', edgecolors=get_color(label))
            try:  # fit exponential decay of photon lags
                p0 = [np.max(Intensity) / 2, phase[f'{route}: I (kHz)'] / 100]
                p, _ = curve_fit(exp, tau[tau > treshold_ms], Intensity[tau > treshold_ms],
                                 bounds=[[0, 0], [np.inf, 1e5]], p0=p0)
                pulses.loc[j, f'{route}: bg (kHz)'] = p[1] * phase['length (ticks)'] / pulseperiod
                if plot: plt.plot(tau, exp(tau, p[0], p[1]), color=get_color(label))
            except:  # | no photons with lag > treshold
                pulses.loc[j, f'{route}: bg (kHz)'] = None

    if plot:
        tio.format_plot(xtitle='lagtime (ms)', ytitle='#photons', title=trace.filename(), ylog=True,
                        xrange=[0, np.max(tau)])
    return pulses


def get_bursts(filename, maxGap_s=200e-6, minCounts=10):
    trace = tio.Traces(filename)
    timetags = read_timetags(filename)
    pulses = get_pulses(filename)

    # correct for dutycycle Pulsed Interleaved Excitation
    dutycycle = len(pulses[pulses['wavelength (nm)'] != 0]) / len(pulses)
    maxGap_s /= dutycycle
    minCounts *= dutycycle
    selectedTimetags = select_timetags(timetags, pulses, non_dark=True)

    # identify and number bursts with gaps smaller than maxGap_s
    inBursts = np.diff(selectedTimetags.index, prepend=0) < maxGap_s / ticklength_s
    startIndex = (np.diff(inBursts, prepend=0) == 1)
    burstNr = np.cumsum(startIndex) * inBursts
    selectedTimetags['Burst'] = 0
    selectedTimetags.loc[selectedTimetags.index, 'Burst'] = burstNr

    # remove bursts smaller than minCount and renumber
    burstSize = selectedTimetags['Burst'].value_counts()
    burstsToRemove = burstSize[burstSize <= minCounts].index.values
    selectedTimetags.loc[selectedTimetags['Burst'].isin(burstsToRemove), 'Burst'] = 0
    inBursts = (selectedTimetags['Burst'] > 0)
    startIndex = (np.diff(inBursts, append=0) == 1)
    burstNr = np.cumsum(startIndex) * inBursts
    timetags['Burst'] = 0
    timetags.loc[selectedTimetags.index, 'Burst'] = burstNr

    # save to file
    timetags.to_hdf(trace.filename(ext='h5'), '_timetags', mode='a', format='table')
    return timetags['Burst'].value_counts()


def get_ES(filename, gamma=1.0, save=None):
    trace = tio.Traces(filename)
    timetags = read_timetags(filename)
    pulses = get_pulses(filename)

    # Get burst intensities
    bursts = pd.DataFrame(columns=['G514', 'R514', 'R632'])
    labels = list(bursts)
    for label in labels:
        bursts[label] = select_timetags(timetags, pulses, wavelength=label[1:], route=label[0])[
            'Burst'].value_counts()

    # clean up
    bursts = bursts.iloc[1:]
    bursts.fillna(0, inplace=True)
    bursts.sort_index(inplace=True)
    bursts.index.name = 'number'

    # Calculate burst duration
    start = timetags.drop_duplicates('Burst', keep='first')[1:]
    end = timetags.drop_duplicates('Burst', keep='last')[:-1]
    bursts.insert(0, 't (s)', start.index * ticklength_s)
    bursts.insert(1, 'dt (ms)', 1000.0 * (end.index - start.index) * ticklength_s)

    # Subtract background
    for label in labels:
        try:
            background = float(pulses[pulses['wavelength (nm)'] == int(label[1:])][f'{label[0]}: bg (kHz)'])
            bursts[label] -= background * bursts['dt (ms)']
        except:
            pass  # No background detected

    bursts['E'] = bursts['R514'] / (gamma * bursts['G514'] + bursts['R514'])
    bursts['S'] = (gamma * bursts['G514'] + bursts['R514']) / (gamma * bursts['G514'] + bursts['R514'] + bursts['R632'])

    bursts.to_hdf(trace.filename(ext='h5'), '_bursts', mode='a', format='table')

    plt.hexbin(bursts['E'], bursts['S'], extent=[0, 1, 0, 1], gridsize=25, cmap="Blues")

    if save is not None:
        plotfilename = trace.filename(ext=save, addition='ES')
    else:
        plotfilename = None
    tio.format_plot('E', 'S', trace.plotname(), xrange=[0, 1], yrange=[0, 1], aspect=1, scale_page=0.5,
                    save=plotfilename)

    return bursts


# import fretbursts as fb
#
# fb.dplot()

if __name__ == "__main__":
    #   Intended as examples. Make your own script and import this module to your script.
    #   You can copy the code below as a start of your own script
    pd.set_option('display.expand_frame_repr', False)

    files = [
        # r'C:\Users\noort\data\\FCS\data\DNAFRET_nr1_ALEX100nsR10dG100d_000.t3r',
        # r'C:\Users\noort\data\\FCS\data\DNAFRET_nr1_contR30_000.t3r',
        r'C:\Users\noort\data\\FCS\data\DNAFRET_nr1_contG50_000.t3r',
    ]
    files = [rf'C:/Users/noort/data/FCS/data/DNAFRET_nr4_contG50_000.t3r']
    files = [r'C:\Users\noort\data\FCS\bursts\DNAcy3b_7p6pM_ALEX25usR30G50_000.t3r']
    # files = [rf'C:/Users/noort/data/FCS/data/DNAFRET_nr4_ALEX100nsR10dG100d_000.t3r']

    # files = glob.glob('C:/Users/noort/data/FCS/data/*.t3r', recursive=True)

    labels = ['R632xR632', 'G514xG514', 'R632xG514', 'R514xR514']  # , 'R514xG514', 'R632xR514', 'R632xR514']

    for i, filename in enumerate(files):
        # print(T3R_to_h5(filename).T)
        # print(f'{get_pulses(filename)}\n\n')
        correlate(filename, labels)
        fit_correlation(filename, labels)
        # plot_correlation(filename, labels, save='png')
        # get_bursts(filename, maxGap_s=2000e-6, minCounts=100)
        # print(get_ES(filename, gamma=1.0))
        plot_timetrace(filename, step_s=0.025, range_s=[0, 30], Imax_kHz=10)

    # d = tio.aggregate_parameters(files, save='results2')