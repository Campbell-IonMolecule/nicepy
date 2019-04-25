# -*- coding: utf-8 -*-
"""
Created on Fri Jul 1 10:55:41 2016

Data Maninpulation Tools

@author: Gary Chen
"""

import scipy as _s
import scipy.ndimage as _ndimage
import scipy.ndimage.filters as _filters
import numpy as _n
import labrad
from labrad.units import THz
import matplotlib.pyplot as _p
from scipy import constants

_c = labrad.connect()
_d = _c.data_vault


def cd_data_vault(year, month, day):
    month = '%02d' % month
    day = '%02d' % day
    _d.cd([''])
    _d.cd(['%s' % (year), '%s' % (month), '%s_%s_%s' % (year, month, day)])


def collect_reaction_data(year, month, day, run):
    """
    Parses data from data vault and returns lists with relevant data.

    returns: times, data, sums, axial_power, radial_power, par, date, run

    Parameters
    ----------
    year: year data set is taken
    month: month data set is taken
    day: day data set is taken
    run: data run number
    """

    month = '%02d' % month
    day = '%02d' % day
    date = '%s_%s_%s' % (year, month, day)
    _d.cd([''])
    _d.cd(['%s' % (year), '%s' % (month), '%s_%s_%s' % (year, month, day), 'Be_reaction', 'images_%s' % (run)])

    data = []
    sums = []
    times = []
    axial_power = []
    radial_power = []
    images = len(_d.dir()[1])

    _d.open('00001 - images')
    par = _d.get_parameters()

    for image in range(images):
        _d.open('%05d - images' % (image + 1))

        dat = list(_d.get())
        if len(dat) is not 1:

            try:
                t = _d.get_parameter('step time')
            except:
                t = 0
            try:
                pa = _d.get_parameter('power_axial')
            except:
                pa = 0
            try:
                pr = _d.get_parameter('power_radial')
            except:
                pr = 0
            times.append(t)
            data.append(dat)
            axial_power.append(pa)
            radial_power.append(pr)
            sums.append(_n.sum(dat))
        else:
            None

    return _n.array(times), _n.array(data), sums, axial_power, radial_power, par, date, run


def plot_arb_images(label, data, label_string):
    """
    Neatly displays arbitrary numbers of images from the camera

    returns fig

    Parameters:
    -----------
    label: array of values that each image is labeled by, e.g. time
    data: array of arrays of image data
    label_string: string describing label, e.g. 's'
    """
    length = len(data)
    columns = 10
    if length % columns is not 0:
        rows = length / columns + 1
    else:
        rows = length / columns

    fig = _p.figure()
    fig.set_figheight(rows * 5)
    fig.set_figwidth(10)

    for i in range(length):
        ax = fig.add_subplot(rows, columns, i + 1)
        ax.matshow(data[i], vmin=_n.min(data), vmax=_n.max(data))
        ax.set_title('%s\n%.1f%s' % (i, label[i], label_string))
        if i % 10 is 0:
            ax.set_xticks([])
            ax.set_ylabel('pixels')
        else:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    return fig


def plot_raw_images(array):
    """
    Plots raw camera images

    returns fig

    Parameters:
    -----------
    array: raw array produced by collect_reaction_data
    """

    fig = plot_arb_images(array[0], array[1], 's')

    return fig


def plot_counted_image_set(array, ratio=0.5, neighborhood_size=20):
    """
    Plots raw images with identified ions circled

    Parameters:
    -----------
    array: raw array produced by collect_reaction_data
    ratio: threshold ratio between max and min values in data
    neighborhood_size: searched region around identified peak (ions default is 20)
    """

    times = array[0]
    data = array[1]
    length = len(data)

    x = []
    y = []
    num_objects = []

    for i in range(length):
        a = count_ions(data[i], ratio, neighborhood_size)
        num_objects.append(a[1])
        x.append(a[2])
        y.append(a[3])

    columns = 10
    if length % columns is not 0:
        rows = length / columns + 1
    else:
        rows = length / columns

    fig = _p.figure()
    fig.set_figheight(rows * 5)
    fig.set_figwidth(columns)
    for i in range(length):
        ax = fig.add_subplot(rows, columns, i + 1)
        ax.matshow(data[i], vmin=_n.min(data), vmax=_n.max(data))
        ax.scatter(x[i], y[i], facecolors='none', edgecolors='r')
        ax.set_title('%s\n%.1fs\n%s' % (i, times[i], num_objects[i]))
        ax.autoscale(False)
        if i % 10 is 0:
            ax.set_xticks([])
            ax.set_ylabel('pixels')
        else:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    return fig


def count_ions(data, ratio=0.5, neighborhood_size=20):
    """
    counts ions in a 2D image

    threshold is automatically determined by the ratio of filtered maximum and minimum

    Parameters:
    -----------
    data: array of values
    ratio: threshold ratio between max and min values in data
    neighborhood_size: searched region around identified peak (ions default is 20)
    """

    data_max = _filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = _filters.minimum_filter(data, neighborhood_size)

    maximum = _n.max(data_max)
    minimum = _n.min(data_min)

    threshold = (maximum - minimum) * ratio + minimum

    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = _ndimage.label(maxima)

    slices = _ndimage.find_objects(labeled)

    x, y = [], []
    for dy, dx in slices:
        x_center = (dx.start + dx.stop - 1) / 2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1) / 2
        y.append(y_center)

    return labeled, num_objects, x, y


def plot_counted_image(data, ratio=0.5, neighborhood_size=20):
    """
    Plots an image with circled ions identified

    Parameters:
    data: array of values
    neighborhood_size: searched region around identified peak (ions default is 20)
    """
    labeled, num_objects, x, y = count_ions(data, ratio, neighborhood_size)

    fig = _p.figure()
    ax = fig.add_subplot(111)
    ax.matshow(data)

    ax.autoscale(False)
    ax.scatter(x, y, facecolors='none', edgecolors='r')

    ax.set_title('%s' % (num_objects))
    ax.set_xticks([])
    fig.tight_layout()

    return fig


def ind_image_background(data, top_range, bot_range):
    """
    Calculates background of an image by extrapolating values from empty regions
    reaching from both the top and bottom of the image

    Parameters:
    data: array of values
    top_range: extent that blank area reaching in from the top
    bot_range: extent that blank area reaches in from the bottom
    """

    data = _n.array(data)

    length = _n.shape(data)[0]
    width = _n.shape(data)[1]

    top = data[:top_range]
    bot = data[-bot_range:]
    summ = _n.sum(top) + _n.sum(bot)
    ratio = float(len(data)) / (len(top) + len(bot))
    tot_summ = ratio * summ

    num = length * width

    ind_background = tot_summ / num

    norm_data = _n.array(data) - ind_background

    return norm_data


def sum_ion_fluorescence(data, ratio=0.5, neighborhood_size=20):
    """
    Sums fluorescence per ion of a data set by identifying the each ion
    and summing background subtracted counts in a region around each

    returns avg, std, ion_sums, ion_regions

    Parameters:
    -----------
    data: array of values
    neighborhood_size: searched region around identified peak (ions default is 20)
    """

    data = _n.array(data)

    length = _n.shape(data)[0]
    width = _n.shape(data)[1]

    labeled, num_objects, x, y = count_ions(data, ratio, neighborhood_size)

    top_range = _n.min(y) - neighborhood_size
    bot_range = _n.max(y) + neighborhood_size

    norm_data = _n.array(ind_image_background(data, top_range, length - bot_range))

    ylim = [(yy - neighborhood_size / 2, yy + neighborhood_size / 2) for yy in y]
    xlim = [(xx - neighborhood_size / 2, xx + neighborhood_size / 2) for xx in x]

    ions = map(list, zip(*[ylim, xlim]))

    ion_regions = []
    for ion in ions:
        a = norm_data[ion[0][0]:ion[0][1], ion[1][0]:ion[1][1]]
        ion_regions.append(a)

    ion_sums = [_n.sum(ion) for ion in ion_regions]

    avg = _n.mean(ion_sums)
    std = _n.std(ion_sums)

    return avg, std, ion_sums, ion_regions


def sum_ion_fluorescence_set(data, ratio=0.5, neighborhood_size=20):
    avg = []
    std = []
    sums = []
    regions = []

    for dat in data:
        a = sum_ion_fluorescence(dat, ratio, neighborhood_size)
        avg.append(a[0])
        std.append(a[1])
        sums.append(a[2])
        regions.append(a[3])

    return avg, std, sums, regions


def plot_ion_regions(avg, std, sums, regions):
    fig = _p.figure()
    region_num = len(regions)
    ion_num = [len(region) for region in regions]
    m_ions = _n.max(ion_num)
    fig.set_figwidth(region_num * 1.5)
    fig.set_figheight(m_ions * 1.8)

    ions = []
    for i, region in enumerate(regions):
        for ii, ion in enumerate(region):
            iii = 1 + i + ii * region_num
            ax = fig.add_subplot(m_ions, region_num, iii)
            ax.matshow(ion)
            if ii is 0:
                ax.set_title('%s\n[%.2E]\n%.2E' % (i, avg[i], sums[i][ii]))
            else:
                ax.set_title('%.2E' % (sums[i][ii]))
            ax.set_xticks([])
            ax.set_yticks([])
    #    fig.tight_layout()
    return fig


def plot_fluorescence_per_ion(power, fluor, init=[1e6, 0.1], fit=False):
    fig = _p.figure()
    ax = fig.add_subplot(111)

    if fit is True:
        fitt, cov = _s.optimize.curve_fit(r_pp, power, fluor[0], p0=init)
        x = _n.linspace(0, 50, 100)
        ax.scatter(power, fluor[0] / fitt[0])
        ax.plot(x, r_pp(x, *fitt) / fitt[0], color='r')

        ax.set_title('State Fraction')
        ax.set_xlabel('dumped power (mW)')
        ax.set_ylabel('P-state fraction')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.set_xlim(0, 50)
        ax.set_ylim(0)

        fig.tight_layout()

        # p state fraction
        states = r_pp(power, *fitt) / fitt[0]
        frac = [power, states]

        return fig, fitt, cov, frac

    else:
        ax.scatter(power, fluor[0])
        ax.set_title('Fluorescence per ion')
        ax.set_xlabel('dumped power (mW)')
        ax.set_ylabel('p-state fraction')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.set_xlim(0, max(power))
        ax.set_ylim(0)

        fig.tight_layout()

        return fig


def collect_iso_ions(array, images, ratio=0.5, neighborhood_size=20):
    time = array[0]
    data = array[1]
    ti = [time[image] for image in images]
    dat = [data[image] for image in images]
    fluor = sum_ion_fluorescence_set(dat, ratio, neighborhood_size)

    fig1 = plot_counted_image_set([ti, dat])
    fig2 = plot_ion_regions(*fluor)

    if 0 in fluor[1]:
        weighted = _n.mean(fluor[0])
    else:
        weighted = wmom(fluor[0], _n.divide(1, fluor[1]))

    return fig1, fig2, weighted, fluor


def wmom(arrin, weights_in, inputmean=None, calcerr=True, sdev=False):
    """
    NAME:
      wmom()

    PURPOSE:
      Calculate the weighted mean, error, and optionally standard deviation of
      an input array.  By default error is calculated assuming the weights are
      1/err^2, but if you send calcerr=True this assumption is dropped and the
      error is determined from the weighted scatter.

    CALLING SEQUENCE:
     wmean,werr = wmom(arr, weights, inputmean=None, calcerr=False, sdev=False)

    INPUTS:
      arr: A numpy array or a sequence that can rate_constants converted.
      weights: A set of weights for each elements in array.
    OPTIONAL INPUTS:
      inputmean:
          An input mean value, around which them mean is calculated.
      calcerr=False:
          Calculate the weighted error.  By default the error is calculated as
          1/sqrt( weights.sum() ).  If calcerr=True it is calculated as sqrt(
          (w**2 * (arr-mean)**2).sum() )/weights.sum()
      sdev=False:
          If True, also return the weighted standard deviation as a third
          element in the tuple.

    OUTPUTS:
      wmean, werr: A tuple of the weighted mean and error. If sdev=True the
         tuple will also contain sdev: wmean,werr,wsdev

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU

   """

    # no copy made if they are already arrays
    arr = _n.array(arrin, ndmin=1, copy=False)

    # Weights is forced to rate_constants type double. All resulting calculations
    # will also rate_constants double
    weights = _n.array(weights_in, ndmin=1, dtype='f8', copy=False)

    wtot = weights.sum()

    # user has input a mean value
    if inputmean is None:
        wmean = (weights * arr).sum() / wtot
    else:
        wmean = float(inputmean)

    # how should error rate_constants calculated?
    if calcerr:
        werr2 = (weights ** 2 * (arr - wmean) ** 2).sum()
        werr = _n.sqrt(werr2) / wtot
    else:
        werr = 1.0 / _n.sqrt(wtot)

    # should output include the weighted standard deviation?
    if sdev:
        wvar = (weights * (arr - wmean) ** 2).sum() / wtot
        wsdev = _n.sqrt(wvar)
        return wmean, werr, wsdev
    else:
        return wmean, werr


def exponential(t, a, b, c):
    return a * _n.exp(b * t) + c


def r_pp(power, a, radius):
    i = power / (_n.pi * radius ** 2)
    s = i / 0.85
    ro = a / 2 * s / (1 + s + 4 * (1) ** 2)

    return ro


def image_fit(array, start, stop=100, init=[1e8, 0, 1e7]):
    """
    Fits and plots image sums to exponential fit

    Parameters:
    -----------
    array: raw array produced by collect_reaction_data
    start: start time for fit
    stop: stop time for fit (default 100)
    init: initial fit parameters
    """
    times = array[0]
    sums = array[2]
    par = array[5]
    date = array[6]
    run = array[7]

    length = len(sums)

    fig = _p.figure()
    ax = fig.add_subplot(111)
    ax.plot(times, sums, 'o')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Total Counts (counts)')
    ax.set_title('Be$^+$ Fluorescence %s run %s' % (date, run), fontsize=16)
    ax.set_ylim(_n.min(sums), _n.max(sums))
    fig.set_figheight(8)
    fig.set_figwidth(12)

    param = [array[5][i + 1] for i in range(len(array[5]) - 1)]

    idx = _n.searchsorted(times, start)
    idx2 = _n.searchsorted(times, stop)

    fit, cov = _s.optimize.curve_fit(exponential, times[idx:idx2], sums[idx:idx2], p0=init)

    domain = _n.linspace(times[0], times[-2], 100)
    ax.plot(domain, exponential(domain, *fit))

    fig.text(0.4, 0.9, param)

    fig.text(0.4, 0.85, 'fit parameters: a = %e, b = %e, c = %e' % (fit[0], fit[1], fit[2]))

    fig.text(0.4, 0.8, 'fit region: %.2fs to %.2fs' % (times[idx], times[idx2]))
    fig.tight_layout()
    _p.show()

    return fit, cov, fig


def collect_beam_data(year, month, day, axis, scan):
    """
    Parses data from data vault and returns lists with relevant data.

    returns: times, voltages, detunings, temperatures, parameters (of first trace)

    Parameters
    ----------
    year: year of data set (xxxx)
    month: month of data set (xx)
    day: day of data set (xx)
    axis: axis relative to beam
          0 = longitudinal
          1 = transverse
    scan: scan number
    """

    axes = ['longitudinal', 'transverse']

    month = '%02d' % month
    day = '%02d' % day
    _d.cd([''])
    _d.cd(['%s' % (year), '%s' % (month), '%s_%s_%s' % (year, month, day), '%s_yb_beam_linescan' % (axes[axis]),
           'scan_%s' % (scan)])

    traces = range(len(_d.dir()[1]))
    times = []
    voltages = []
    temperatures = []
    detunings = []

    # Set up initial parameters
    _d.open('00001 - traces')
    first = list(_d.get())
    points = range(len(first))
    channels = range(len(first[0]))
    channels.remove(0)
    par = _d.get_parameters()

    try:
        center_freq = _d.get_parameter('center frequency') / THz
    except:
        center_freq = 752.4518

    try:
        channel_locs = [_d.get_parameter('channel %s location' % (channel)) for channel in channels]
    except:
        channel_locs = _n.zeros(len(first[0]))

    buffer_species = _d.get_parameter('buffer species')
    flow = _d.get_parameter('buffer flow rate')

    # Collect data into arrays
    for trace in traces:
        _d.open('%05d - traces' % (trace + 1))
        dat = list(_d.get())
        try:
            d = _d.get_parameter('frequency')
        except:
            d = _d.get_parameter('detuning')
        ct = _d.get_parameter('cell temperature')
        try:
            st = _d.get_parameter('shield temperature')
        except:
            st = 0
        t = [dat[point][0] for point in points]
        c = []
        for channel in channels:
            v = [dat[point][channel] for point in points]
            c.append(v)

        voltages.append(c)
        times.append(t)
        temperatures.append([ct, st])
        detunings.append(d)

    # Average and normalize data
    zeros = [next(index for index, time in enumerate(times[trace]) if time > 0.0) for trace in traces]

    voltages = [[_n.array(voltages[trace][channel]) - _n.mean(voltages[trace][channel][0:zeros[trace] - 100])
                 for trace in traces] for channel in range(len(channels))]

    detunings = _n.array(detunings)
    real_174 = 145.505 * 1e-6  # THz
    real_det = detunings + real_174

    return times, voltages, real_det, temperatures, par, traces, center_freq, channel_locs, buffer_species, flow


def plot_traces(times, voltages, detunings, num):
    """
    Plots the raw data of each scan

    Parameters
    ----------
    times: times array in ns
    voltages: voltages array in mV
    detunings: detunings from center in THz
    num: number of traces
    """

    channels = range(len(voltages[0]))
    traces = range(num)
    colors = ['b', 'r']

    for trace in traces:
        for channel in channels:
            fig = _p.figure()

            ax = fig.add_subplot(111)
            ax.plot(times[trace], voltages[trace][channel], color='%s' % (colors[channel]))
            ax.set_ylabel('Absorption Signal (mV)', color='%s' % (colors[channel]))
            ax.set_xlabel('time (s)')
            ax.set_title('Channel %s Absorption @ %sTHz Detuning' % (channel + 1, detunings[trace]))


def average_traces(array, step=10, *args):
    """
    Calculates sums over an interval

    Parameters
    ----------
    array: array of data
    step: number of data points to sum over
    """

    times = array[0]
    voltages = array[1]
    detunings = array[2]
    temperatures = array[3]
    par = array[4]
    traces = array[5]
    center_freq = array[6]
    channel_locs = array[7]
    buffer_species = array[8]
    flow = array[9]

    channels = range(len(voltages))
    zeros = [next(index for index, time in enumerate(times[trace]) if time > 0.0) for trace in traces]
    points = range(zeros[0], len(times[0]), step)

    step_time = [times[trace][zeros[0]::step] for trace in traces]
    avgs = [[[_n.mean(voltages[channel][trace][point:point + step]) for point in points]
             for trace in traces] for channel in channels]

    return step_time, avgs, detunings, temperatures, par, traces, center_freq, channel_locs, buffer_species, flow


def plot_spectrum_2D(array, step, channel):
    traces = range(len(sums[0]))
    vals = [sums[channel][trace][step] for trace in traces]
    maximum = _n.amax(_n.array(sums[channel]))
    minimum = _n.amin(_n.array(sums[channel]))

    fig = _p.figure()
    ax = fig.add_subplot(111)
    ax.plot(detunings, vals)
    ax.set_title('Spectrum @ %ss' % (step_time[0][step]))
    ax.set_xlabel('Detuning (THz)')
    ax.set_ylabel('Absorption (arb)')
    ax.set_ylim([minimum, maximum])


def plot_spectrum_3D(array, title_font=20, axis_size=16):
    time = array[0]
    voltages = array[1]
    detunings = array[2]
    temperatures = array[3]
    center_freq = array[6]
    channel_locs = array[7]
    buffer_species = array[8]
    flow = array[9]

    t = _n.transpose(temperatures)
    channels = range(len(voltages))

    figs = []

    for channel in channels:
        x, y = _n.meshgrid(detunings, time[0])

        fig = _p.figure(figsize=(15, 10))
        ax = fig.add_subplot(111)
        c = ax.contourf(x, y, _n.transpose(voltages[channel]))
        ax.set_title('%ssccm %s @ %s' % (flow, buffer_species, channel_locs[channel]), fontsize=title_font)
        ax.set_ylabel('Time (s)', fontsize=axis_size)
        ax.set_xlabel('Detuning (THz)', fontsize=axis_size)
        fig.colorbar(c, orientation='vertical')
        figs.append(fig)

    fig = _p.figure()
    ax = fig.add_subplot(111)
    ax.plot(detunings, t[0], label='cell')
    ax.plot(detunings, t[1], label='shield')
    ax.set_ylabel('Temperatures')
    ax.set_xlabel('Detuning (THz)')
    ax.legend()

    return figs


def gaussian(f, f_0, amp, FWHM):
    return amp * _n.exp(-4 * _n.log(2) * (f - f_0) ** 2 / FWHM ** 2)


def yb_spectrum_fixed(f, ff, ww, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    isotopes = [('176', -509.310),
                ('173 (5/2)', -253.418),
                ('174', 0),
                ('173 (3/2)', 515.975),
                ('172', 533.309),
                ('173 (7/2)', 587.986),
                ('171 (3/2)', 832.436),
                ('171 (1/2)', 1153.696),
                ('170', 1192.393)]

    spectrum = gaussian(f, ff + isotopes[0][1], a0, ww) + \
               gaussian(f, ff + isotopes[1][1], a1, ww) + \
               gaussian(f, ff + isotopes[2][1], a2, ww) + \
               gaussian(f, ff + isotopes[3][1], a3, ww) + \
               gaussian(f, ff + isotopes[4][1], a4, ww) + \
               gaussian(f, ff + isotopes[5][1], a5, ww) + \
               gaussian(f, ff + isotopes[6][1], a6, ww) + \
               gaussian(f, ff + isotopes[7][1], a7, ww) + \
               gaussian(f, ff + isotopes[8][1], a8, ww)

    return spectrum