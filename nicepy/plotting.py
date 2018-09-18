from matplotlib.pyplot import setp as _setp
from scipy.constants import golden as _golden
import numpy as _np


def format_plot(fig, ax=False, **kwargs):
    format_fig(fig, **kwargs)
    if ax:
        format_ax(ax, **kwargs)


def format_fig(fig, height=8, width=8 * _golden, dpi=300):
    """
    Generic formatting of figure
    :param fig:
    :param height:
    :param width:
    :param dpi:
    :return:
    """
    fig.set_size_inches(width, height)
    fig.set_dpi(dpi)


def format_ax(ax, font=30, label=24, tick=24, box_out=False, width=2):
    """
    Generic formatting of plot axes
    :param ax:
    :param font:
    :param label:
    :param tick:
    :param box_out:
    :param width:
    :return:
    """
    ax.tick_params(labelsize=tick)
    ax.yaxis.offsetText.set_fontsize(label)
    _setp(ax.get_lines(), linewidth=width)
    ax.xaxis.offsetText.set_fontsize(label)
    if box_out is True:
        ax.legend(fontsize=label, frameon=False, loc='upper left', bbox_to_anchor=(1.04, 1))
    elif box_out is False:
        ax.legend(fontsize=label, frameon=False, loc=0)
    ax.xaxis.label.set_size(label)
    ax.yaxis.label.set_size(label)
    ax.title.set_fontsize(font)


def print_fit(funct, fit, cov=False):
    """
    Prints the fit of a function with standard errors
    :param funct: function used for fitting
    :param fit: list of fit values
    :param cov: covariance matrix
    :return:
    """
    params = funct.__code__.co_varnames
    idx = 1
    if params[0] == 'self':
        idx = 2
    diag_cov = _np.sqrt(_np.diag(abs(cov)))
    params = params[idx:len(fit)+1]
    print('%s fit parameters' % funct.__name__)
    if cov is not False:
        for i, param in enumerate(params):
            print('%s: %.02e, %.02e' %(param, fit[i], diag_cov[i]))
    else:
        for i, param in enumerate(params):
            print('%s: %.02e' %(param, fit[i]))
