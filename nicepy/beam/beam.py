from numpy import pi as _pi
from numpy import sqrt as _sqrt
from numpy import array as _array
from nicepy import u as _u

_e = 4.8e-10  # CGS elementary charge
_kb = 1.38e-16  # CGS boltzmann constant


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

    output = _sqrt(top / bot)

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


def pressure_or_density(value, temperature):
    """
    Calculates the reciprocal of input value based on the ideal gas law
    :param value: value to convert (use pint units)
    :param temperature: temperature of gas
    :return: reciprocal of input value
    """

    kind = None
    units = {'pressure': _u.torr, 'density': _u.cm**-3}
    for key, val in units.items():
        a = value/val
        if a.units.dimensionless:
            kind = key
            break
    if kind == 'pressure':
        top = value
        bot = _u.k * temperature
        base = _u.cm**-3
    elif kind == 'density':
        top = value * _u.k * temperature
        bot = 1
        base = _u.torr
    elif kind is None:
        return 'Check input units'

    output = top / bot

    return output.to(base)


def pressure_to_density(pressure, temperature):
    """
    Calculates molecular density from a pressure
    :param pressure:
    :param temperature:
    :return:
    """
    pressure = _array(pressure)
    top = pressure
    bot = _u.k * temperature

    output = top/bot

    return output


def density_to_pressure(density, temperature):
    """
    Calculates pressure from a molecular density
    :param density:
    :param temperature:
    :return:
    """
    density = _array(density)

    output = density * _u.k * temperature

    return output


# def quadrupole(self):
#     a = self.langevin
#     b = 0.0146*self.vals[self.species]['Q']**2
#     c = _np.sqrt(self.kb * self.temp) * self.e * self.vals[self.species]['alpha']**1.5
#     self.quadrupole = a*b/c


# def aqo(self):
#     k = self.langevin + self.quadrupole
#     self.AQO = k
