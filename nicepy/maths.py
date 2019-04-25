from scipy import special as _special
import numpy as _np


def reduced_chi_squared(funct, x, y, yerr, fit, params):
    """
    Calculates reduced chi squared value
    :param funct: function fitted to
    :param x: x values
    :param y: y values
    :param yerr: y error values
    :param fit: fit parameters
    :param params: number of parameters
    :return:
    """
    diff = _np.array(funct(x, *fit)) - y
    top = _np.sum((diff / yerr) ** 2)
    bot = len(y) - params

    output = top / bot

    return output


def chi2_pdf(x, k):
    """

    :param x: chi2 value
    :param k: degrees of freedom
    :return: probability
    """
    top = x ** (k/2-1) * _np.exp(-x / 2)
    bot = 2**(k/2)*_special.gamma(k/2)

    output = top / bot

    return output


def chi2_cdf(x, k):
    """

    :param x: chi2 value
    :param k: degrees of freedom
    :return: cumulative probability of getting value associated with chi2
    """
    top = _special.gammainc(k/2, x/2)
    bot = _special.gamma(k/2)

    output = top / bot

    return output


def linear(x, m, b):
    """
    Linear model
    :param x: independent variable
    :param m: slope
    :param b: intercept
    :return: y value
    """
    x = _np.array(x)

    output = m * x + b

    return output


def constant(x, n):
    """
    Constant value model
    :param x: independent variable
    :param n: value
    :return: value array
    """

    output = [n for _ in x]

    return output


def gaussian(x, a, mu, sigma):
    """
    Gaussian model
    :param x: independent variable
    :param a: height
    :param mu: center
    :param sigma: standard deviation
    :return: y value
    """
    x = _np.array(x)

    output = a * _np.exp(- (x - mu)**2 / (2 * sigma ** 2))

    return output


def lorentzian(x, mu, gamma):
    """
    Lorentzian model
    :param x: independent variable
    :param mu: center
    :param gamma: width
    :return: y value
    """
    x = _np.array(x)

    top = gamma
    bot = 2 * _np.pi * ((x - mu) ** 2 + (gamma / 2) ** 2)

    output = top / bot

    return output


def weighted_mean_and_error(values, errors):
    """

    :param values:
    :param errors:
    :return:
    """
    values = _np.array(values)
    errors = _np.array(errors)
    avg = _np.average(values, weights=1 / errors)
    err = _np.sqrt(1 / _np.sum(1 / errors))

    return avg, err


def tsallis(ion_temp, avg_temp, n):
    """
    Non-normalized probability of an ion at ion-temp using a Tsallis distribution
    :param ion_temp: temperature of ion (K)
    :param avg_temp: average temperature of ions (K)
    :param n: average harmonic oscillator level
    :return: value
    """
    kb = 1.38e-23
    energy = ion_temp * kb
    top = (n - 3) * (n - 2) * (n - 1) * energy ** 2
    bot = 2 * (n * kb * avg_temp) ** 3 * (1 + energy / (n * kb * avg_temp)) ** n

    output = top / bot

    return output


def thermal(ion_temp, avg_temp):
    """
    Non-normalized probability of an ion at ion-temp using a thermal distribution
    :param ion_temp: temperature of ion (K)
    :param avg_temp: average temperature of ions (K)
    :return: value
    """
    kb = 1.38e-23
    energy = ion_temp * kb
    a = 2 / (kb * avg_temp) ** (3 / 2)
    b = _np.sqrt(energy / _np.pi)
    c = _np.exp(-energy / (kb * avg_temp))

    output = a * b * c

    return output


def ratio(vals):
    """
    Returns ratio of first value to sum of list
    :param vals: list or array of values
    :return: ratio of first value to sum of list
    """
    output = vals[0] / _np.sum(vals)

    return output


def ratio_error(vals, errors):
    tot = _np.sum(vals)
    tot_error = _np.sqrt(_np.sum(errors))

    a = errors[0] / tot
    b = tot_error * vals[0] / tot ** 2

    output = _np.sqrt(_np.sum(_np.array([a, b]) ** 2))

    return output