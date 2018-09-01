from scipy.stats import sem as _sem
from nicepy.data import DataObj as _DataObj
import matplotlib.pyplot as _plt

from nicepy.labrad import dv as _dv
from nicepy.labrad import *


class IonImages(_DataObj):
    def __init__(self):
        data = {}
        params = {}
        _DataObj.__init__(self, data, params)
        self.d['images'] = []
        self.d['sums'] = []
        self.d['sterr'] = []

    def __setitem__(self, key, item):
        self.p[key] = item

    def _get_params(self):
        p = _dv.get_parameters()
        if type(p) is not None:
            for i, ii in p:
                if i not in self.p.keys():
                    self.p[i] = []
                self.p[i].append(ii)
        else:
            pass

    def get_data(self, year, month, day, run):
        """
        Gets data and parameters from datavault
        stores raw data, sums, and standard deviation

        :param year: year of data int
        :param month: month of data int
        :param day: day of data int
        :param run: run of data int
        :return: data object
        """
        cd(year, month, day, run, kind='images')
        images = _dv.dir()[1]
        for image in images:
            _dv.open(image)
            data = list(_dv.get())
            self.d['images'].append(data)
            self.d['sums'].append(sum(data))
            self.d['sterr'].append(_sem(data))
            self._get_params()

    def show_images(self, idx=None):
        """
        Plots a set of images defined by idx
        :param idx: list of indices
        :return: fig, ax
        """
        ax = None
        if idx is not None:
            idx = list(idx)
            times = [self.p['step time'][i] for i in idx]
            images = [self.d['images'][i] for i in idx]
        else:
            times = self.p['step time']
            images = self.d['images']
        length = len(images)
        columns = 10
        if length % columns is not 0:
            rows = length / columns + 1
        else:
            rows = length / columns

        fig = _plt.figure()
        fig.set_figheight(rows * 5)
        fig.set_figwidth(10)
        expmin = min(images)
        expmax = max(images)
        for i, image in enumerate(images):
            ax = fig.add_subplot(rows, columns, i + 1)
            ax.matshow(image, vmin=expmin, vmax=expmax)
            ax.set_title('%s\n%.1f(s)' % (i, times[i]))
            if i % 10 is 0:
                ax.set_xticks([])
                ax.set_ylabel('pixels')
            else:
                ax.set_xticks([])
                ax.set_yticks([])

        return fig, ax

    def show_image(self, idx):
        """
        Plots a single image from data set
        :param idx: image index
        :return: fig, ax
        """
        fig = _plt.figure()
        ax = fig.add_subplot(111)
        ax.matshow(self.d['images'][idx])
        ax.set_xticks([])

        return fig, ax
