import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt


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
        self.idx = None

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
            if '.' in val:
                try:
                    val = float(val)
                except ValueError:
                    pass
            else:
                try:
                    val = int(val)
                except ValueError:
                    pass
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
        self.raw.plot.line(x=x, y='Volts', title='%s' % title, xlim=(3, 40), color='black', ax=ax, **kwargs)
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
        self.peaks = None

    def _get_tofs(self, **kwargs):
        """

        :param kwargs:
        :return:
        """
        self._tof_list = []
        for filename in self.filenames:
            t = TofData(filename, self.params, **kwargs)
            self._tof_list.append(t)

        tof_objs = []
        for t in self._tof_list:
            temp = t.params.copy()
            temp['tof'] = t
            tof_objs.append(temp)
        self.tof_objs = _pd.DataFrame(tof_objs)

    def _get_raw(self):
        """

        :return:
        """
        temp_list = []
        for t in self._tof_list:
            temp = t.raw
            for key, val in t.params.items():
                temp[key] = val
            temp_list.append(temp)
        self.raw = _pd.concat(temp_list)

    def _get_tof_peaks(self, masses, **kwargs):
        """

        :param masses:
        :param kwargs:
        :return:
        """
        for t in self._tof_list:
            t.get_peaks(masses, **kwargs)

    def get_peaks(self, masses, **kwargs):
        """

        :param masses:
        :param kwargs:
        :return:
        """
        self._get_tof_peaks(masses, **kwargs)
        temp_list = []
        for t in self._tof_list:
            temp = _pd.concat([t.peaks, t.params])
            temp_list.append(temp)
        temp = _pd.concat(temp_list, axis=1)
        self.peaks = temp.transpose()
