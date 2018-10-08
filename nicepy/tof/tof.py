import numpy as _np
from scipy.stats import sem as _sem
import pandas as _pd
import matplotlib.pyplot as _plt
from nicepy import format_fig as _ff, format_ax as _fa


class TofData:
    """
    General class for TOF data
    """

    def __init__(self, filename, params, norm=True, noise_range=(3, 8), bkg_range=(3, 8), fluor=True, factor=0.92152588, offset=-0.36290086):
        """

        :param filename:
        :param params:
        :param norm:
        :param noise_range:
        :param bkg_range:
        :param fluor:
        :param factor:
        :param offset:
        """
        self.filename = filename
        self.idx = False
        self.norm = norm
        self.noise_range = noise_range
        self.bkg_range = bkg_range
        self.factor = factor
        self.offset = offset
        self.fluor = fluor
        self._get_data(filename)
        self._subtract_bkg()
        self._get_noise()
        self._get_params(filename, params)
        self.peaks = None

    def _get_data(self, filename):
        """

        :param filename:
        :return:
        """
        dat = list(_np.loadtxt(filename))
        fluor = dat.pop()
        if self.fluor is False:
            fluor = 1
        center = int(len(dat) / 2)

        time = _np.array([s for s in dat[:center]])
        mass = self._time_to_mass(time, self.factor, self.offset)

        raw = _np.array([-i / fluor for i in dat[center:]]) / fluor

        raw = _pd.DataFrame({'Time': time, 'Mass': mass, 'Volts': raw})

        if self.norm is True:
            tot = raw['Volts'].sum()
            raw['Volts'] = raw['Volts']/tot

        self.raw = raw

    def _subtract_bkg(self):
        """

        :return:
        """
        temp = self._select_range('Mass', self.bkg_range[0], self.bkg_range[1])['Volts']
        m = temp.mean()
        self.raw['Volts'] = self.raw['Volts'] - m

    def _get_noise(self):
        """

        :return:
        """
        temp = self._select_range('Mass', self.noise_range[0], self.noise_range[1])['Volts']
        n = temp.std()
        self.noise = n

    def _get_params(self, filename, params):
        """

        :param filename:
        :param params:
        :return:
        """
        listed = filename.replace('.txt', '').split('_')
        temp = {key: listed[val] for key, val in params.items()}
        self.params = {}
        for key, val in temp.items():
            if key.lower() == 'version':
                val = val.lower()
                val = val.replace('v', '')
            else:
                pass
            if '.' in val:
                try:
                    val = float(val)
                except ValueError:
                    pass
            else:
                try:
                    val = int(val)
                except ValueError:
                    val = val.lower()
            self.params[key] = val

        self.params = _pd.Series(self.params)

    def _select_range(self, column, lower, upper):
        """
        Selects part of data that is between values upper and lower in column
        :param column: column name to be used to bound
        :param lower: lower value in column
        :param upper: upper value in column
        :return: parsed data frame
        """
        temp = self.raw[(self.raw[column] <= upper) & (self.raw[column] >= lower)]
        return temp

    def _get_closest(self, column, value):
        """

        :param column:
        :param value:
        :return:
        """
        temp = self.raw.loc[(self.raw[column] - value).abs().idxmin()]
        return temp

    def _get_range(self, mass, pk_range=(-80, 80)):
        """

        :param mass:
        :param pk_range:
        :return:
        """
        idx = self._get_closest('Mass', mass).name
        lower = idx + pk_range[0]
        if lower < 0:
            lower = 0
        upper = idx + pk_range[1]
        if upper > self.raw.index.max():
            upper = self.raw.index.max()

        return lower, upper

    def _get_peak(self, lower, upper):
        """

        :param lower:
        :param upper:
        :return:
        """
        temp = self.raw.loc[range(lower, upper + 1)]
        p = temp['Volts'].sum()
        if p < self.noise:
            p = 0

        return p

    def get_peaks(self, masses, **kwargs):
        """

        :param masses:
        :param kwargs:
        :return:
        """
        self.peaks = {}
        self.idx = {}
        for key, val in masses.items():
            lower, upper = self._get_range(val, **kwargs)
            self.peaks[key] = self._get_peak(lower, upper)
            self.idx[key] = (lower, val, upper)
        self.peaks = _pd.Series(self.peaks)
        self.idx = _pd.Series(self.idx)

    @staticmethod
    def _time_to_mass(time, factor, offset):
        mass = [(t - factor) ** 2 + offset for t in time]

        return mass

    def show(self, x='Mass', shade=True, **kwargs):
        """

        :param x:
        :param shade:
        :param kwargs:
        :return:
        """
        fig, ax = _plt.subplots()
        title = {key: val for key, val in self.params.items()}
        self.raw.plot.line(x=x, y='Volts', title='%s' % title, color='black', ax=ax, **kwargs)
        if shade is True:
            if self.idx is not False:
                for key, val in self.idx.items():
                    idx_range = self.raw.loc[range(val[0], val[2] + 1)]
                    ax.fill_between(idx_range[x], idx_range['Volts'], label=key, alpha=0.5)
                ax.legend(loc=0)
            else:
                pass
        else:
            pass
        ax.legend(loc=0)

        return fig, ax


class TofSet:

    def __init__(self, filenames, params, **kwargs):
        """

        :param filenames:
        :param params:
        :param kwargs:
        """
        self.filenames = filenames
        self.params = params
        self._get_tofs(**kwargs)
        self._get_raw()
        self.idx = False
        self.peaks = None

    def _get_tofs(self, **kwargs):
        """

        :param kwargs:
        :return:
        """
        tof_list = []
        for filename in self.filenames:
            t = TofData(filename, self.params, **kwargs)
            tof_list.append(t)

        tof_objs = []
        for t in tof_list:
            temp = t.params.copy()
            temp['tof'] = t
            tof_objs.append(temp)
        temp = _pd.DataFrame(tof_objs)
        temp.set_index(list(self.params.keys()), inplace=True)
        self.tof_objs = temp.sort_index()

    def _get_raw(self):
        """

        :return:
        """
        # temp = self.tof_objs.copy()
        # temp['raw'] = [t.raw for t in temp['tof']]
        # self.raw = temp.drop('tof', axis=1)

        temp = self.tof_objs.copy()
        raw = []
        for tof in temp['tof']:
            t = tof.raw
            for key, val in tof.params.items():
                t[key] = val
            raw.append(t)

        self.raw = _pd.concat(raw)
        self.raw.set_index(list(self.params.keys()), inplace=True)
        self.raw.sort_index(inplace=True)

    def _get_tof_peaks(self, masses, **kwargs):
        """

        :param masses:
        :param kwargs:
        :return:
        """
        for t in self.tof_objs['tof']:
            t.get_peaks(masses, **kwargs)
        self.idx = t.idx

    def get_peaks(self, masses, **kwargs):
        """

        :param masses:
        :param kwargs:
        :return:
        """
        self._get_tof_peaks(masses, **kwargs)
        temp_list = []
        for t in self.tof_objs['tof']:
            temp = _pd.concat([t.peaks, t.params])
            temp_list.append(temp)
        temp = _pd.concat(temp_list, axis=1)
        self.peaks = temp.transpose()
        self.peaks.set_index(list(self.params.keys()), inplace=True)
        self.peaks.sort_index(inplace=True)
        # self.peaks['total'] = self.peaks.sum(axis=1)

    def get_raw_means(self, levels=None):
        if levels is None:
            levels = []
            for key in ['version', 'delay']:
                if key in self.params.keys():
                    levels.append(key)
        else:
            pass

        self.levels = levels

        grouped = self.tof_objs.groupby(levels)
        temp_mean = []
        temp_error = []
        for indices, group in grouped:
            times = _np.mean([tof.raw['Time'] for tof in group['tof']], axis=0)
            masses = _np.mean([tof.raw['Mass'] for tof in group['tof']], axis=0)
            volts = _np.mean([tof.raw['Volts'] for tof in group['tof']], axis=0)
            errors = _sem([tof.raw['Volts'] for tof in group['tof']], axis=0)
            df_mean = _pd.DataFrame({'Time': times, 'Mass': masses, 'Volts': volts})
            df_error = _pd.DataFrame({'Time': times, 'Mass': masses, 'Volts': errors})
            if type(indices) is not tuple:
                indices = [indices]
            for key, index in zip(levels, indices):
                df_mean[key] = index
                df_error[key] = index
            temp_mean.append(df_mean)
            temp_error.append(df_error)
        self.raw_means = _pd.concat(temp_mean)
        self.raw_errors = _pd.concat(temp_error)
        self.raw_means.set_index(levels, inplace=True)
        self.raw_errors.set_index(levels, inplace=True)

    def get_peak_means(self, levels=None):
        """

        :param levels: list of index keys
        :return:
        """
        if levels is None:
            levels = []
            for key in ['version', 'delay']:
                if key in self.params.keys():
                    levels.append(key)
        else:
            pass
        ignore = [key for key in self.params.keys() if key not in levels]
        temp = self.peaks.reset_index(ignore, drop=True)
        self.peak_means = temp.groupby(levels).mean()
        self.peak_errors = temp.groupby(levels).sem()
        self.peak_total = self.peak_means.sum(axis=1)

    def show_means(self, x='Mass', levels=None, shade=True, fmt=False, box_out=True, **kwargs):
        self.get_raw_means(levels)

        levels = []
        figs = []
        axs = []

        for indices, group in self.raw_means.groupby(level=self.levels):
            if type(indices) is not tuple:
                idxs = [indices]
            fig, ax = _plt.subplots()
            ax.set_title('%s' % {key: val for key, val in zip(self.levels, idxs)})
            group.loc[indices].plot.line(x=x, y='Volts', ax=ax, color='black', **kwargs)
            if shade is True:
                if self.idx is not False:
                    for key, val in self.idx.items():
                        temp = group.loc[indices].reset_index()
                        idx_range = temp.loc[range(val[0], val[2] + 1)]
                        ax.fill_between(idx_range[x], idx_range['Volts'], label=key, alpha=0.5)
                    ax.legend(loc=0)
                else:
                    pass
            else:
                pass
            if fmt is True:
                _ff(fig)
                _fa(ax, box_out=box_out)
            levels.append(indices)
            figs.append(fig)
            axs.append(ax)

        l = list(self.raw_means.index.names)

        if type(indices) is tuple:
            multi = _pd.MultiIndex.from_tuples(levels, names=l)
            self.means_fig_ax = _pd.DataFrame({'fig': figs, 'ax': axs}, index=multi)
        else:
            self.means_fig_ax = _pd.DataFrame({'fig': figs, 'ax': axs}, index=levels)
            self.means_fig_ax.index.name = l[0]

    def show_peaks(self, norm=False, fmt=False, box_out=True, **kwargs):
        if norm is False:
            total = 1
        else:
            total = self.peak_total

        means = self.peak_means.div(total, axis=0)
        errors = self.peak_errors.div(total, axis=0)

        indexes = means.index
        if len(indexes.names) > 1:
            levels = []
            figs = []
            axs = []
            for indices, group in means.groupby(level=0):
                fig, ax = _plt.subplots()
                ax.set_title('%s' % (indices))
                group.loc[indices].plot.line(yerr=errors.loc[indices], marker='o', ax=ax, **kwargs)
                if fmt is True:
                    _ff(fig)
                    _fa(ax, box_out=box_out)
                levels.append(indices)
                figs.append(fig)
                axs.append(ax)
                l = list(self.peak_means.index.names)
                l.remove('delay')

                if type(indices) is tuple:
                    multi = _pd.MultiIndex.from_tuples(levels, names=l)
                    self.peaks_fig_ax = _pd.DataFrame({'fig': figs, 'ax': axs}, index=multi)
                else:
                    self.peaks_fig_ax = _pd.DataFrame({'fig': figs, 'ax': axs}, index=levels)
                    self.peaks_fig_ax.index.name = l[0]

        else:
            fig, ax = _plt.subplots()
            means.plot.line(yerr=errors, marker='o', ax=ax, **kwargs)
            if fmt is True:
                _ff(fig)
                _fa(ax, box_out=box_out)
            self.peaks_fig_ax = _pd.Series({'fig': fig, 'ax': ax})
