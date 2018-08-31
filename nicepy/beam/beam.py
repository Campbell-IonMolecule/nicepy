from numpy import pi as _pi
from nicepy import u as _u


def mb_mean_velocity(temperature, mass):
    """
    Calculates the mean velocity of a Maxwell Boltzmann Distribution
    :param temperature: temperature of distribution
    :param mass: mass of gas particle
    :return: velocity
    """
    top = 8 * _u.k * temperature
    bot = mass * _pi

    return top / bot
