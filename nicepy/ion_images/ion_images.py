from scipy.stats import sem as _sem
from nicepy.data import DataObj as _DataObj
import matplotlib.pyplot as _plt
import numpy as _np
import pandas as _pd

from nicepy.labrad import dv as _dv
from nicepy import labrad


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
        labrad.cd(year, month, day, run, kind='images')
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


def collect_ion_data(year, month, day, run):

    labrad.cd(year, month, day, run)

    data = {}
    for image in _dv.dir()[1]:
        _dv.open(image)
        if not data:
            data['images'] = []
            for key, val in _dv.get_parameters():
                data[key] = [val]
        else:
            for key, val in _dv.get_parameters():
                data[key].append(val)
        dat = _np.array(_dv.get())
        data['images'].append(dat)

    df = _pd.DataFrame(data)

    return df


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
