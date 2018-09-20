import pandas as _pd
import numpy as _np
import matplotlib.pyplot as _plt
from nicepy.beam import pressure_or_density as _p_or_d
from nicepy import u as _u


class RGA:
    
    def __init__(self, filenames, bkg_range=(22, 26)):
        """
        Collects and cleans up RGA data
        :param filenames: paths to RGA data files in .txt format
        :param bkg_range: amu range to set background to
        """
        if type(filenames) is not list:
            filenames = [filenames]
        self.filenames = filenames
        self.data = None
        self._get_data()
        self._subtract_bkg(bkg_range)

    def _get_data(self):
        """
        Gets RGA data from .txt file
        :return: DataFrame of averaged rga files
        """
        data = []
        for filename in self.filenames:
            d = _pd.read_table(filename, header=17, delimiter=',', names=['Mass', 'Pressure'], index_col=False)
            data.append(d)

        if len(data) > 1:
            temp = _pd.concat(data, axis=1)
        else:
            temp = _pd.DataFrame(data)
        temp['Pressure (Torr)'] = temp['Pressure'].mean(axis=1)
        temp['Mass (m/z)'] = data[0]['Mass']
        temp['Sigma (Torr)'] = temp['Pressure'].sem(axis=1)

        self.data = temp[['Mass (m/z)', 'Pressure (Torr)', 'Sigma (Torr)']].copy()

    def _subtract_bkg(self, mass_range):
        """
        Subtracts background from data set by the mass range
        :param mass_range: tuple of range of masses to select for background subtraction
        :return: DataFrame with background subtracted
        """
        upper = max(mass_range)
        lower = min(mass_range)
        a = self._select_range('Mass (m/z)', lower, upper)
        bkg = abs(a['Pressure (Torr)'].mean())
        self.data['Pressure (Torr)'] = self.data['Pressure (Torr)'] - bkg
        # self.data['Pressure (Torr)'].loc[self.data['Pressure (Torr)'] < 0] = 0

    def show(self, **kwargs):
        """
        Outputs plot of RGA trace
        :param kwargs: plot arguments
        :return: fig, ax
        """

        fig, ax = _plt.subplots()

        self.data.plot.line(x='Mass (m/z)', y='Pressure (Torr)', yerr='Sigma (Torr)', ax=ax, **kwargs)

        return fig, ax

    def _select_range(self, column, lower, upper):
        """
        Selects part of data that is between values upper and lower in column
        :param column: column name to be used to bound
        :param lower: lower value in column
        :param upper: upper value in column
        :return: parsed data frame
        """
        temp = self.data[(self.data[column] <= upper) & (self.data[column] >= lower)]
        return temp

    def get_avg_pressure(self, lower, upper):
        """
        Gets an average pressure and standard error from a mass range
        :param lower: lower mass value
        :param upper: upper mass value
        :return: average and standard error for masses in defined mass range
        """
        temp = self._select_range('Mass (m/z)', lower, upper)['Pressure (Torr)']
        p = temp.mean()
        if len(temp) == 0:
            s = abs(p * 0.1)
        else:
            s = temp.sem()
        return p, s

    def water_iso(self, alpha=0.768, beta=0.185, gamma=0.047, h2o_range=(18.1, 18.2), hod_range=(19.1, 19.2), d2o_range=(20.1, 20.2)):

        masses = {'H2O': h2o_range, 'HOD': hod_range, 'D2O': d2o_range}

        def h2o(p1, p2, p3):
            a = 1 / alpha
            b = p1 - p2 * beta / (2 * alpha) - p3 * beta / alpha
    
            output = a * b
    
            return output
    
        def hod(p2):
            output = p2 / alpha
    
            return output
    
        def d2o(p3):
            output = p3 / alpha
    
            return output
    
        def h2o_error(sp1, sp2, sp3):
            a = sp1 / alpha
            b = (beta / alpha ** 2) * sp3
            c = (beta / (2 * alpha ** 2)) * sp2
    
            output = _np.sqrt(_np.sum(_np.array([a, b, c]) ** 2))
    
            return output
    
        def hod_error(sp2):
            a = sp2 / alpha
    
            output = _np.sqrt(_np.sum(_np.array([a]) ** 2))
    
            return output
    
        def d2o_error(sp3):
            a = sp3 / alpha
    
            output = _np.sqrt(_np.sum(_np.array([a]) ** 2))
    
            return output

        vals = {}
        errors = {}

        for key, val in masses.items():
            p, s = self.get_avg_pressure(val[0], val[1])
            vals[key] = p
            errors[key] = s

        data = {mass: {} for mass in masses.keys()}

        p1 = vals['H2O']
        p2 = vals['HOD']
        p3 = vals['D2O']
        sp1 = errors['H2O']
        sp2 = errors['HOD']
        sp3 = errors['D2O']
        data['H2O']['val'] = h2o(p1, p2, p3)
        data['HOD']['val'] = hod(p2)
        data['D2O']['val'] = d2o(p3)
        data['H2O']['error'] = h2o_error(sp1, sp2, sp3)
        data['HOD']['error'] = hod_error(sp2)
        data['D2O']['error'] = d2o_error(sp3)

        self.water_iso_pressures = _pd.DataFrame(data)
        temp = {key: _p_or_d(self.water_iso_pressures.loc['val'][key] * _u.torr, 300 * _u.K).magnitude for key in self.water_iso_pressures.loc['val'].keys()}
        self.water_iso_densities = _pd.Series(temp)
