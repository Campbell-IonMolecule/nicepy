from os import listdir
import numpy as _np
import scipy as _sp
import copy
import matplotlib as _plt
from matplotlib import rcParams
import itertools
from nicepy.data import DataObj as _DataObj


class TofFiles(object):
    """
    Class for cleaning up and parsing TOF data files from directory

    Parsing with self.exact [(start, stop, string)]
    """

    def __init__(self, excludes=[], includes=[], exact=[], endswith='.txt', startswith='', delimiter='_',
                 remove=['.txt'], params={}):
        self.excludes = excludes
        self.includes = includes
        self.exact = exact
        self.endswith = endswith
        self.startswith = startswith
        self.delimiter = delimiter
        self.remove = remove
        self.params = params
        self.raw_files = {}
        self.parsed_files = []

    def _raw_files(self, directory):
        """
        Gathers raw files that satisfy start and end conditions
        :param directory: directory of TOF files
        :return: None
        """
        all_files = listdir(directory)
        for fstring in all_files:
            if fstring.endswith(self.endswith):
                if fstring.startswith(self.startswith):
                    fil = tuple(self._split_file(fstring))
                    self.raw_files[fil] = directory + fstring

    def _split_file(self, fstring):
        """
        Removes unwanted substrings and splits files according to delimiter
        :param fstring: file string
        :return: array of file substrings
        """
        for i in self.remove:
            if i in fstring:
                fstring = fstring.replace(i, '')
        fil = fstring.split(self.delimiter)
        return fil

    def _parse_files(self):
        """
        Parses wanted files according to excludes and includes arrays
        :return:
        """
        if type(self.excludes) is not list:
            self.excludes = [self.excludes]
        if type(self.includes) is not list:
            self.includes = [self.includes]
        if type(self.exact) is not list:
            self.exact = [self.exact]

        for fil, fstring in self.raw_files.items():
            good = True
            if len(self.includes) is not 0:
                for inc in self.includes:
                    if inc not in fil:
                        good = False
                        break
            if len(self.excludes) is not 0:
                for exc in self.excludes:
                    if exc in fil:
                        good = False
                        break
            if len(self.exact) is not 0:
                for i, ii, exa in self.exact:
                    if exa not in fil[i: ii]:
                        good = False
                        break

            if good is True:
                p = self._get_params(fil)
                t = _DataObj(fstring, p)
                self.parsed_files.append(t)

    def _get_params(self, fil):
        """
        Takes self.params key as key and value as index
        :param fil: array of file substrings
        :return:
        """
        params = {}
        for key, idx in self.params.items():
            if '.' in fil[idx]:
                value = float(fil[idx])
            else:
                try:
                    value = int(fil[idx])
                except ValueError:
                    value = fil[idx]

            params[key] = value
        return params

    def get_files(self, directory):
        self._raw_files(directory)
        self._parse_files()


class TofData(_DataObj):
    """
    General class for TOF data
    """

    def __init__(self, fstring, params, fluor=True, factor=0.92152588, offset=-0.36290086):
        """
        Init function
        :param fstring: file path string
        :param params: dictionary of data parameters
        """
        _DataObj.__init__(self, {}, params)
        self.factor = factor
        self.offset = offset
        self.fluor = fluor
        self.collect_data(fstring)

    def _file_data(self, fil):
        dat = list(_np.loadtxt(fil))
        fluor = dat.pop()
        if self.fluor is False:
            fluor = 1
        center = int(len(dat) / 2)

        self.d['time'] = _np.array([s for s in dat[:center]])
        self._time_to_mass()

        raw = _np.array([-i / fluor for i in dat[center:]])
        m = _np.mean(raw[-300:])
        self.d['volts'] = raw - m
        self.sum = _np.sum(self.d['volts'])

    def _time_to_mass(self):
        self.d['mass'] = [(t - self.factor) ** 2 + self.offset for t in self.d['time']]

    def collect_data(self, fstring):
        self._file_data(fstring)


class TofSet(TofFiles):
    """
    Class for sets of TOF data objects
    """

    def __init__(self, excludes=[], includes=[], exact=[], endswith='.txt', startswith='', delimiter='_',
                 remove=['.txt'], params={}, masses=None, filter=False, filter_f = 1, factor=0.92152588,
                 offset=-0.36290086,
                 fluor=True, fitting='Full', upper=30, lower=-20):
        TofFiles.__init__(self, excludes=excludes, includes=includes, exact=exact, endswith=endswith,
                          startswith=startswith, delimiter=delimiter, remove=remove, params=params)
        self.raw = []
        self.avg = []
        self.peaks = []
        self.masses = masses
        self.filter = filter
        self.filter_f = filter_f
        self.factor = factor
        self.offset = offset
        self.fluor = fluor
        self.fitting = fitting
        self._upper = upper
        self._lower = lower

    def _get_data(self):
        if len(self.parsed_files) == 0:
            print('No files to get data from')
            return
        for toffile in self.parsed_files:
            t = TofData(toffile.d, toffile.p, fluor=self.fluor, factor=self.factor, offset=self.offset)
            self.raw.append(t)

    def _filter_data(self):
        sums = [raw.sum for raw in self.raw]
        mean = abs(_np.mean(sums))
        std = _np.std(sums)
        for raw in self.raw:
            if self.filter is 'lower':
                if abs(raw.sum) < mean - std * self.filter_f:
                    self.raw.remove(raw)
            elif self.filter is 'upper':
                if abs(raw.sum) > mean + std * self.filter_f:
                    self.raw.remove(raw)

    def get_raw(self):
        self._get_data()
        if self.filter is not False:
            self._filter_data()

    @staticmethod
    def _group(data, omit=None):
        """
        Groups data by parameters
        :param data: list of _DataObj's
        :param omit: parameter to omit from grouping
        :return: list of _DataObj's with an attribute list of _DataObj's
        """
        temp = {}
        for dat in data:
            ps = copy.deepcopy(dat.p)
            if omit in ps.keys():
                del ps[omit]
            params = tuple(ps.values())
            if params not in temp.keys():
                temp[params] = _DataObj([], ps)
            temp[params].d.append(dat)
        templist = [value for value in temp.values()]
        return templist

    def _get_background(self, time, bkg):

        if self.fitting is 'Sum':
            sample = _np.sum([abs(self._lower), abs(self._upper)])
        else:
            sample = 25

        chunks = []
        for i in range(0, len(bkg), sample):
            chunks.append(bkg[i:i + sample])

        if self.fitting is 'Sum':
            areas = []
            for chunk in chunks:
                areas.append(_np.sum(chunk))
        else:
            if self.fitting is 'Full':
                funct = _profile
            if self.fitting is 'Iso':
                funct = _gaussian
            times = time[:sample]
            areas = []
            for chunk in chunks:
                try:
                    fit, cov = _sp.optimize.curve_fit(funct, times, chunk, p0=[0.0, times[6], 0.01])
                    areas.append(fit[0] * fit[2])
                except RuntimeError:
                    pass

        hist = _np.histogram(areas)

        x = []
        for i in range(len(hist[1])):
            if i < len(hist[1]) - 1:
                x.append(_np.mean([hist[1][i], hist[1][i + 1]]))
            else:
                pass

        w = (max(x)-min(x))/2
        p = max(hist[0])
        try:
            fit, cov = _sp.optimize.curve_fit(_gaussian, x, hist[0], p0=[p, 0, w])
            # x0 = fit[1]
            sigma = fit[2]
        except RuntimeError:
            fig = _plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x, hist[0])

            sigma = 0

        return sigma, hist

    def _average(self, group):
        p = group.p
        temp = {}
        for obj in group.d:
            for key, value in obj.d.items():
                if key not in temp.keys():
                    temp[key] = []
                temp[key].append(value)

        tempdata = {}
        for key, data in temp.items():
            m = _np.mean(data, axis=0)
            tempdata[key] = m

        if len(temp['volts']) == 1:
            tempdata['sterr'] = temp['volts'][0] / 20
        else:
            tempdata['sterr'] = _sp.stats.sem(temp['volts'], axis=0)
        avg = _DataObj(tempdata, p)
        self.avg.append(avg)

    def _fit_noise(self, end=-300):

        for group in self.avg:
            time = group[0].d['time']
            bkg = []
            for avg in group:
                bkg += list(avg.d['volts'][end:])

            sigma, hist = self._get_background(time, bkg)
            group.hist = hist
            for avg in group:
                # avg.d['volts'] = avg.d['volts']
                avg.d['noise'] = sigma * 1.3

    def get_averages(self):
        self.get_raw()
        groups = self._group(self.raw)
        for group in groups:
            self._average(group)
        self.avg = sorted(self.avg, key=lambda obj: obj.p['delay'])
        self.avg = self._group(self.avg, 'delay')
        self._fit_noise()

    @staticmethod
    def find(array, value):
        """
        Finds the index of a value in a list
        :param array: array of values
        :param value: desired value
        :return: index of value in array
        """
        idx = _np.argmin(_np.abs(_np.array(array) - value))
        return idx

    @staticmethod
    def _get_p0(time, idxs):
        order = []
        p = []
        for key, idx in idxs.items():
            p += [5e-10, time[idx], 3e-3]
            order.append(key)
        return p, order

    @staticmethod
    def _get_bounds(time, idxs, high=20, low=-20):

        upper = []
        lower = []
        for key, idx in idxs.items():
            # [peak, x0, width]
            upper += [_np.inf, time[idx + high], 1e0]
            lower += [-_np.inf, time[idx + low], 0]

        return (tuple(lower), tuple(upper))

    def _fit_profile(self, avg_obj, dataobj, x, y, sy, p0, order, bounds=(-_np.inf, _np.inf)):

        if self.fitting is 'Full':
            funct = composite_profile
        if self.fitting is 'Iso':
            funct = _gaussian

        try:
            fit, cov = _sp.optimize.curve_fit(funct, x, y,
                                              sigma=sy, p0=p0,
                                              bounds=bounds, absolute_sigma=True)
            cov = _np.sqrt(_np.diag(cov))
        except RuntimeError:
            fit = _np.zeros(len(p0))
            cov = _np.zeros(len(p0))
            print('fit error')

        dataobj.d['fit'].append(fit)
        dataobj.d['cov'].append(_np.sqrt(_np.diag(cov)))

        f = [fit[i:i + 3] for i in range(0, len(fit), 3)]
        c = [cov[i:i + 3] for i in range(0, len(fit), 3)]

        if 'ind fit' not in avg_obj.d.keys():
            avg_obj.d['ind fit'] = {}
        if 'ind cov' not in avg_obj.d.keys():
            avg_obj.d['ind cov'] = {}
        for i, species in enumerate(order):
            avg_obj.d['ind fit'][species] = f[i]
            avg_obj.d['ind cov'][species] = c[i]

    def _get_fit_area(self, avg_obj):

        temp_a = {}
        temp_sa = {}
        noise = avg_obj.d['noise']

        for species in self.masses.keys():
            f = avg_obj.d['ind fit'][species]
            c = avg_obj.d['ind cov'][species]
            a = f[0] * f[2]
            sa = _np.sqrt(_np.sum(_np.array([c[0] * f[2], c[2] * f[0]]) ** 2))

            if a > noise:
                temp_a[species] = a
                temp_sa[species] = sa
            else:
                temp_a[species] = 0
                temp_sa[species] = noise
        return temp_a, temp_sa

    def get_peaks(self):
        """
        gets peak values from averaged TOF data

        different methods may be used: 'Full' fits the full ion profile
        'Iso' fits only the gaussian around the peak
        'Sum' does a rudimentary summation of the peaks in a set range.

        :return: self.peaks
        """
        if type(self.masses) is not dict:
            print('TofSet.masses dictionary not correctly defined')
            return
        idxs = {mass: self.find(self.raw[0].d['mass'], value) for mass, value in self.masses.items()}
        # print(idxs)
        for group in self.avg:
            time = group.d[0].d['time']
            p0, order = self._get_p0(time, idxs)
            bounds = self._get_bounds(time, idxs)

            t = _DataObj({}, group.p)

            for key in ['delay', 'fit', 'cov', 'norm tot', 'norm peaks']:
                t.d[key] = []

            for species in self.masses.keys():
                t.d[species] = {}
                for key in ['peaks', 'sterr']:
                    t.d[species][key] = []

            for data in group.d:
                params = copy.deepcopy(data.p)
                delay = params['delay']
                t.d['delay'].append(delay)

                y = data.d['volts']
                sy = data.d['sterr']

                t.d['norm tot'].append(_np.sum(y))

                if self.fitting is 'Full':

                    self._fit_profile(data, t, time, y, sy, p0, order, bounds)
                    temp_a, temp_sa = self._get_fit_area(data)

                    t.d['norm peaks'].append(_np.sum(list(temp_a.values())))

                    for species in self.masses.keys():
                        t.d[species]['peaks'].append(temp_a[species])
                        t.d[species]['sterr'].append(temp_sa[species])

                if self.fitting is 'Iso':
                    temp_a = []

                    for species in self.masses.keys():
                        i = idxs[species] - 0.008/(time[1]-time[0]) # 0.008s to the edge of gaussian
                        ii = idxs[species] + 0.008/(time[1]-time[0])

                        p0 = [1e-9, time[idxs[species]], 3e-3]

                        self._fit_profile(data, t, time[i:ii], y[i:ii], sy[i:ii], p0, [species])
                        noise = data.d['noise']
                        f = data.d['ind fit'][species]
                        c = data.d['ind cov'][species]

                        a = f[0] * f[2]
                        sa = _np.sqrt(_np.sum(_np.array([c[0] * f[2], c[2] * f[0]]) ** 2))

                        if a > noise:
                            t.d[species]['peaks'].append(a)
                            t.d[species]['sterr'].append(sa)
                        else:
                            t.d[species]['peaks'].append(0)
                            t.d[species]['sterr'].append(noise)

                        temp_a.append(a)

                    t.d['norm peaks'].append(_np.sum(temp_a))

                if self.fitting is 'Sum':

                    norm = 0
                    data.d['span'] = {}
                    noise = data.d['noise']

                    for species in self.masses.keys():
                        i = idxs[species] + self._lower
                        ii = idxs[species] + self._upper

                        data.d['span'][species] = (i, ii)

                        a = _np.sum(data.d['volts'][i:ii])
                        sa = _np.sqrt(_np.sum(data.d['sterr'][i:ii]) ** 2)
                        if a > noise:
                            t.d[species]['peaks'].append(a)
                            t.d[species]['sterr'].append(sa)
                            norm += a
                        else:
                            t.d[species]['peaks'].append(0)
                            t.d[species]['sterr'].append(noise)

                    t.d['norm peaks'].append(norm)

            for species in self.masses.keys():
                t.d[species]['peaks'] = _np.array(t.d[species]['peaks'])
                t.d[species]['sterr'] = _np.array(t.d[species]['sterr'])
            t.d['delay'] = _np.array(t.d['delay'])
            t.d['norm tot'] = _np.array(t.d['norm tot'])
            t.d['norm peaks'] = _np.array(t.d['norm peaks'])
            self.peaks.append(t)

    def show_raw(self, xlabel='mass', select=None):
        fig = _plt.figure()
        ax = fig.add_subplot(111)
        if xlabel is 'mass':
            ax.set_xlabel('Mass (amu)')
        elif xlabel is 'time':
            ax.set_xlabel('Time (us)')
        ax.set_ylabel('Ion Detector Voltage (arb)')
        for r in self.raw:
            if type(select) is dict:
                good = True
                for key, value in select.items():
                    if type(value) is not list:
                        value = [value]
                    if r.p[key] not in value:
                        good = False
                        break
                if good is True:
                    x = r[xlabel]
                    y = r['volts']
                    ax.plot(x, y, label=r.p)
        return fig, ax

    @staticmethod
    def _return_fig_ax(fig, ax, group):
        figs = {}
        t = [val for val in group.p.values()]
        if len(t) == 1:
            t = t[0]
        else:
            t = tuple(t)
        figs['fig'] = fig
        figs['ax'] = ax
        return t, figs

    def show_avg(self, separate=False, omit=None, fmt=False, limit=False, lines=False):
        """

        :param separate: Boolean to separate delay plots
        :param omit: Dict with values wanted to be omitted: {parameter: [values]}
        :param fmt: Boolean to apply automatic formatting
        :param limit: display horizontal line for noise limit
        :return: Dict of {figure title dict values: (fig, ax)}
        """
        figs = {}
        xlim = None
        if self.masses is not None:
            m = [val for val in self.masses.values()]
            xlim = (min(m) - 1, max(m) + 1)
        if separate is False:
            for group in self.avg:
                good = True
                if omit is not None:
                    for key, value in omit.items():
                        if key in group.p.keys():
                            if type(value) is not list:
                                value = [value]
                            if group[key] in value:
                                good = False
                                break
                if good is True:
                    fig = _plt.figure()
                    ax = fig.add_subplot(111)
                    if xlim is not None:
                        ax.set_xlim(xlim[0], xlim[1])
                    ax.set_xlabel('Mass (m/z)')
                    ax.set_ylabel('Ion Detector Voltage (arb)')
                    ax.set_title(group.p)
                    for dataobj in group:
                        good = True
                        if omit is not None:
                            for key, value in omit.items():
                                if type(value) is not list:
                                    value = [value]
                                if key in dataobj.p.keys():
                                    if dataobj.p[key] in value:
                                        good = False
                                        break
                        if good is True:
                            x = dataobj.d['mass']
                            y = dataobj.d['volts']
                            ax.plot(x, y, label=dataobj.p['delay'])
                    if fmt is True:
                        format_fig(fig)
                        format_ax(ax)
                        # fig.tight_layout()
                    if lines is True:
                        for val in self.masses.values():
                            ax.axvline(val)
                    if limit is True:
                        ax.axhline(group[0].d['noise'])
                    t, f = self._return_fig_ax(fig, ax, group)
                    figs[t] = f
        elif separate is True:
            for group in self.avg:
                for dataobj in group:
                    good = True
                    if omit is not None:
                        for key, value in omit.items():
                            if type(value) is not list:
                                value = [value]
                            if key in group.p.keys():
                                if group.p[key] in value:
                                    good = False
                                    break
                            elif key in dataobj.p.keys():
                                if dataobj.p[key] in value:
                                    good = False
                                    break
                    if good is True:
                        fig = _plt.figure()
                        ax = fig.add_subplot(111)
                        x = dataobj.d['mass']
                        y = dataobj.d['volts']
                        ax.plot(x, y, color='k')
                        ax.set_title(dataobj.p)
                        ax.set_xlim(xlim[0], xlim[1])
                        ax.set_xlabel('Mass (m/z)')
                        ax.set_ylabel('Ion Detector Voltage (arb)')
                        if 'ind fit' in dataobj.d.keys():
                            for species in dataobj.d['ind fit'].keys():

                                t = dataobj.d['time']
                                l, = ax.plot(x, _gaussian(t, *dataobj.d['ind fit'][species]))
                                ax.fill_between(x, _gaussian(t, *dataobj.d['ind fit'][species]), 0, alpha=0.5,
                                                label=species, facecolor=l.get_color())
                                ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.04, 1))
                        elif 'span' in dataobj.d.keys():
                            cgen = itertools.cycle(clist)
                            for species, idx in dataobj.d['span'].items():
                                i = idx[0]
                                ii = idx[1]
                                ax.fill_between(x[i: ii], y[i: ii], label=species, alpha=0.5,
                                                facecolor=next(cgen)['color'])
                                ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.04, 1))
                        if lines is True:
                            for val in self.masses.values():
                                ax.axvline(val)
                        if fmt is True:
                            format_fig(fig)
                            format_ax(ax)
                            # fig.tight_layout()
                        if limit is True:
                            ax.axhline(dataobj.d['noise'])
                        t, f = self._return_fig_ax(fig, ax, dataobj)
                        figs[t] = f
        return figs

    def show_peaks(self, omit=None, fmt=False, norm=False, total=False):
        figs = {}
        for dataobj in self.peaks:
            good = True
            if omit is not None:
                for key, value in omit.items():
                    if type(value) is not list:
                        value = [value]
                    if key in dataobj.p.keys():
                        if dataobj.p[key] in value:
                            good = False
                            break
            if good is True:
                fig = _plt.figure()
                ax = fig.add_subplot(111)
                x = dataobj.d['delay']
                ax.set_xlim(min(x) - 2, max(x) + 2)
                ax.set_xlabel('Delay (s)')
                ax.set_ylabel('Ion Signal (arb)')
                ax.set_title(dataobj.p)
                summed = []
                ssummed = []
                tot = 1
                if norm is not False:
                    tot = dataobj.d['%s' % norm]
                for species in self.masses.keys():
                    y = dataobj.d[species]['peaks']/tot
                    yerr = dataobj.d[species]['sterr']/tot
                    ax.errorbar(x, y, yerr=yerr, label=species)
                    summed.append(y)
                    ssummed.append(yerr)
                if total is True:
                    y = _np.sum(summed, axis=0)
                    yerr = _np.sqrt(_np.sum(_np.array(ssummed) ** 2, axis=0))
                    ax.errorbar(x, y, yerr=yerr, label='Total')
                ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.04, 1))
                if fmt is True:
                    format_fig(fig)
                    format_ax(ax)
                    # fig.tight_layout()
                t, f = self._return_fig_ax(fig, ax, dataobj)
                figs[t] = f
        return figs


def format_fig(fig, height=8, width=8 * _sp.constants.golden, dpi=300):
    fig.set_size_inches(width, height)
    fig.set_dpi(dpi)


def format_ax(ax, font=30, label=24, box_out=True):
    ax.tick_params(labelsize=label)
    ax.yaxis.offsetText.set_fontsize(label)
    ax.xaxis.offsetText.set_fontsize(label)
    if box_out is True:
        ax.legend(fontsize=label, frameon=False, loc='upper left', bbox_to_anchor=(1.04, 1))
    elif box_out is False:
        ax.legend(fontsize=label, frameon=False, loc=0)
    ax.xaxis.label.set_size(label)
    ax.yaxis.label.set_size(label)
    ax.title.set_fontsize(font)


def _gaussian(x, peak, x0, sigma):
    g = peak * _np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    return g


def _linear(x, m, b, x0):
    l = m * (x - x0) + b
    return l


def _cosine(x, peak, f, p, x0):
    c = peak * _np.sin(2 * _np.pi * f * (x - x0) + p)
    return c


def _exponential(x, a, b, x0):
    return a * _np.exp(b * (x - x0))


def _profile(x, peak, x0, width):

    r = []
    edge1 = 0.008
    edge2 = 0.0325
    dip = _linear(x0, 1.00273587, 0.03094452, 0)

    dip_peak = -0.09
    osc_peak1 = 264.269
    osc_peak2 = -174.057
    f1 = 1.38323782
    f2 = 17.1945496
    p1 = -20.5309872e+01
    p2 = 41.10
    decay = -2.5
    w = width*1.2

    for i in x:

        if i < x0 + edge1:
        # if i:
            f = _gaussian(i, peak, x0, width)
            r.append(f)

        elif x0+edge1 <= i < x0+edge2:
            f = _linear(i, -3*peak, 0.1*peak, x0)
            r.append(f)

        elif x0 + edge2 <= i <= dip + w:
        # elif x0+edge2 < i:
            f = _gaussian(i, dip_peak*peak, dip, w)
            r.append(f)

        elif i > dip + w:
            f = (_gaussian(i, dip_peak * peak, dip + w, w)
                 + (_cosine(i, osc_peak1 * peak, f1, p1, x0) + _cosine(i, osc_peak2 * peak, f2, p2, x0)) *
                 _exponential(i, 0.0003 * peak, decay, x0))  # decay maybe -8.54
            r.append(f)

    return r


def composite_profile(x, *args):
    """
    composite profile function for arbitrary numbers of target particles
    :param x: array of time values
    :param args: list of function parameters
    :return: array of TOF voltage values
    """

    v = _np.zeros(len(x))
    for i in range(0, len(args), 3):
        sample = args[i:i + 3]
        v += _np.array(_profile(x, *sample))
    return v


def ratio(data, idxs, mass1, mass2):
    top = _np.mean(data[mass1]['peaks'][idxs[0]:idxs[1]])
    bot = _np.mean(data[mass2]['peaks'][idxs[0]:idxs[1]])
    return top / bot


def ratio_error(data, idxs, mass1, mass2):
    top = _np.mean(data[mass1]['peaks'][idxs[0]:idxs[1]])
    top_err = _np.std(data[mass1]['sterr'][idxs[0]:idxs[1]])
    bot = _np.mean(data[mass2]['peaks'][idxs[0]:idxs[1]])
    bot_err = _np.std(data[mass2]['sterr'][idxs[0]:idxs[1]])

    a = top_err / bot
    b = top * bot_err / bot ** 2

    e = _np.sqrt(_np.sum(_np.array([a, b]) ** 2))
    return e

# New stuff


def _time_to_mass(t, t0, m0):
    output = (t - t0) ** 2 + m0
    return output
