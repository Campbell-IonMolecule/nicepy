from nicepy.data import DataObj as _DataObj
import numpy as _np
from nicepy.maths import reduced_chi_squared as _chi2
from scipy.optimize import curve_fit, OptimizeWarning
import matplotlib.pyplot as _plt
import pandas as _pd


class SharedFit(_DataObj):
    def __init__(self, x, vals, errors, functions, params=False, norm=False):
        """
        Imports data and functions and creates stacks for fitting
        vals, and errors need to rate_constants dicts with keys corresponding to function names

        IMPORTANT:
        Data should rate_constants sorted per the x list, _shared_functions iterates through function list with each non-repeated starting value, so data lists must reflect this ordering.
        e.g. x = [1,1,2,3,4,5,1,...]
        function[0] is passed for x values [1,1,2,3,4,5], the next 1 in x is then passed function[1], etc.

        :param x: independent variable
        :param vals: dict of dependent variables
        :param errors: dict of errors in vals
        :param functions: list of functions
        :param params: parameter description of data
        :param norm: boolean to normalize data
        """

        _DataObj.__init__(self, {'x': x, 'vals': vals, 'errors': errors}, params)
        self._length = len(x)
        self._stacks = {}
        self._names = []
        self.norm = norm
        self._collect_functions(functions)
        self._stack_data()
        self.fit = None
        self.cov = None
        self.chi2 = None

    def _stack_data(self):
        """
        stacks data
        :return:
        """
        self._stacks['vals'] = []
        self._stacks['errors'] = []
        self._stacks['x'] = []
        tot = 1
        if self.norm is True:
            temp = [self.d['vals'][key] for key in self.f.keys()]
            tot = _np.sum(temp, axis=0)
        for key in self.f.keys():
            self._names.append(key)
            self._stacks['vals'] += list(_np.array(self.d['vals'][key]) / tot)
            self._stacks['errors'] += list(_np.array(self.d['errors'][key]) / tot)
            self._stacks['x'] += list(self.d['x'])
        self._stacks['vals'] = _np.array(self._stacks['vals'])
        self._stacks['errors'] = _np.array(self._stacks['errors'])
        self._stacks['x'] = _np.array(self._stacks['x'])

    def _collect_functions(self, functions):
        """
        Stores functions as attributes

        Function names must match self.masses keys
        :param functions: list of functions to rate_constants used in fitting
        :return:
        """
        self.f = {funct.__name__: funct for funct in functions}

    def _shared_functions(self, *args):
        """
        Function passes stacked values through iterated functions
        :param args: fit parameters
        :return: fitting value
        """
        start = args[0][0]
        repeat = True
        output = []
        idx = 0
        for i in args[0]:
            if i == start:
                if repeat is False:
                    idx += 1
                    repeat = True
                else:
                    pass
            else:
                repeat = False
            output.append(self.f[self._names[idx]](i, *args[1:]))
        return output

    def get_fit(self, guess=None, absolute=False, bounds=(0, _np.inf)):
        """
        Get shared fit
        :param guess: fit parameter guess
        :param absolute: boolean to use absolute error
        :param bounds: list like tuples setting fit parameter bounds
        :return:
        """
        self.guess = guess
        if self.guess is None:
            print('Define initial guess in self.guess')
            return
        try:
            self.fit, self.cov = curve_fit(self._shared_functions, self._stacks['x'], self._stacks['vals'], sigma=self._stacks['errors'], p0=self.guess, absolute_sigma=absolute, bounds=bounds)
            self.chi2 = _chi2(self._shared_functions, self._stacks['x'], self._stacks['vals'],
                              self._stacks['errors'], self.fit, len(self.guess))
        except OptimizeWarning or RuntimeError:
            self.fit = self.guess
            self.cov = _np.ones((len(self.guess), len(self.guess))) * _np.inf
            print('no fit')

        self._get_fits()

    def print_fit(self):
        print(self.fit_string())

    def _get_fits(self):
        args = [i for i in self.f.values()]
        args = list(args[0].__code__.co_varnames)
        if 'self' in args:
            args.remove('self')

        self.fits = _pd.Series(self.fit, index=args[1 : len(self.fit)+1])
        self.covs = _pd.Series(_np.sqrt(_np.diag(self.cov)), index=args[1 : len(self.fit)+1])


    def fit_string(self):
        args = [i for i in self.f.values()]
        args = list(args[0].__code__.co_varnames)
        if 'self' in args:
            args.remove('self')
        string = ''
        for arg, fit, cov in zip(args[1:], self.fit, _np.sqrt(_np.diag(self.cov))):
            string += '%s: %.02e, %.02e\n' % (arg, fit, cov)
        if self.chi2 is not None:
            string += 'reduced chi squared: %.02e' % self.chi2
        else:
            string = string[:-1]
        return string

    def show(self, total=False, width=2, text=False, d=False, **kwargs):
        fig = _plt.figure()
        ax = fig.add_subplot(111)

        x = self.d['x']
        xfit = _np.linspace(min(x), max(x), 100)
        total_data = []
        total_error = []
        total_fit = []
        tot = 1
        if self.norm is True:
            temp = [self.d['vals'][key] for key in self.f.keys()]
            tot = _np.sum(temp, axis=0)
        for key in self.f.keys():
            if d:
                s = d[key]
            else:
                s = key
            line = ax.plot(xfit, self.f[key](xfit, *self.fit), label=s+'$^+$', linewidth=width)
            ax.errorbar(x, self.d['vals'][key]/tot,
                        yerr=self.d['errors'][key]/tot,
                        color=line[0].get_color(),
                        linestyle='', marker='o',
                        linewidth=width, capthick=width,
                        label='', **kwargs)

            total_data.append(self.d['vals'][key]/tot)
            total_error.append(self.d['errors'][key]/tot)
            total_fit.append(self.f[key](xfit, *self.fit))
        if total is True:
            total_data = _np.sum(total_data, axis=0)
            total_fit = _np.sum(total_fit, axis=0)
            total_error = _np.sqrt(_np.sum(_np.array(total_error) ** 2, axis=0))
            line = ax.plot(xfit, total_fit, label='Total')
            ax.errorbar(x, total_data, yerr=total_error, color=line[0].get_color())
        if text is True:
            string = self.fit_string().replace(',', ' $\pm$')
            ax.text(1, 0.5, string, transform=fig.transFigure)
        ax.set_xlabel('Delay (s)')
        ax.set_ylabel('Ion Signal (arb)')
        ax.set_title(self.p)
        ax.legend(frameon=False, loc=0)
        ax.set_xlim(min(x) - 1, max(x) + 1)
        return fig, ax
