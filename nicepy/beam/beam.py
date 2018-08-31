from numpy import pi as _pi
from numpy import sqrt as _sqrt
from nicepy import u as _u


def mb_mean_velocity(temperature, mass):
    """
    Calculates the mean velocity of a Maxwell Boltzmann distribution
    :param temperature: temperature of distribution
    :param mass: mass of gas particle
    :return: velocity
    """
    top = 8 * _u.k * temperature
    bot = mass * _pi

    output = _sqrt(top / bot)

    return output


def mb_prob_velocity(temperature, mass):
    """
    Calculates the most probable velocity of a Maxwell Boltzmann distribution
    :param temperature: temperature of distribution
    :param mass: mass of gas particle
    :return: velocity
    """
    top = 2 * _u.k * temperature
    bot = mass

    output = _sqrt(top / bot)

    return output


def max_supersonic_velocity(temperature, mu, gamma):
    """
    Maximum velocity of a supersonic beam
    :param temperature: temperature of initial gas sample
    :param mu: reduced mass of gas particles
    :param gamma: heat capacity ratio (1+1/f) where f is the number of degrees of freedom
    :return: velocity
    """
    top = 2 * _u.k * temperature * gamma
    bot = mu * (gamma - 1)

    output = top / bot

    return output


def reaction_temperature(mu, velocity):
    """
    Calculates reaction temperature
    :param mu: reduced mass of reactants
    :param velocity: relative velocity of reactants
    :return: reaction temperature
    """
    top = mu * velocity ** 2
    bot = 2 * _u.k

    output = top / bot

    return output
