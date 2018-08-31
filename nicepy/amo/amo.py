from numpy import array as _array


def rho_power(power, a, b, detuning=10, lw=10):
    """
    Transition rate of a line addressed by a detuned laser as a function of laser power
    :param power: laser power
    :param a: overall scaling parameter
    :param b: intensity scaling parameter
    :param detuning: detuning from resonance (MHz)
    :param lw: transition linewidth (MHz)
    :return: transition rate
    """
    power = _array(power)
    s = power * b
    top = a * s
    bot = 2 * (1 + s + 4 * (detuning / lw) ** 2)

    output = top/bot

    return output


def rho_detuning(detuning, power, a, b, lw=10):
    """
    Transition rate of a line addressed bby a detuned laser as a function of detuning
    :param detuning: detuning from resonance (MHz)
    :param power: laser power
    :param a: overall scaling parameter
    :param b: intensity scaline parameter
    :param lw: transition linewidth (MHz)
    :return: transition rate
    """
    detuning = _array(detuning)
    s = power * b
    top = a * s
    bot = 2 * (1 + s + 4 * (detuning / lw) ** 2)

    output = top / bot

    return output
