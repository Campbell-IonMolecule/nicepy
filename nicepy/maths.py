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


def weighted_mean(vals, sigma):
    """

    :param vals: values
    :param sigma: errors
    :return: error
    """
    vals = _np.array(vals)
    sigma = _np.array(sigma)

    top = _np.sum(vals / sigma)
    bot = _np.sum(1 / sigma)

    output = top / bot

    return output


def weighted_mean_error(sigma):
    """

    :param sigma: errors
    :return:
    """
    sigma = _np.array(sigma)

    top = len(sigma)
    bot = _np.sum(1 / sigma)

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


def gaussian(x, a, mu, sigma):
    """
    Gaussian model
    :param x: independent variablbe
    :param a: height
    :param mu: mean
    :param sigma: standard deviation
    :return:
    """
    x = _np.array(x)

    output = a * _np.exp((x - mu)**2 / (2 * sigma ** 2))

    return output
