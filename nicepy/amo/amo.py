from numpy import array as _array
from numpy import sqrt as _sqrt
from numpy import pi as _pi

_e = 4.8e-10  # CGS elementary charge
_kb = 1.38e-16  # CGS boltzmann constant


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


def reduced_mass(mass1, mass2):
    """
    Calculates reduced mass
    :param mass1: mass 1
    :param mass2: mass 2
    :return: reduced mass
    """
    top = mass1 * mass2
    bot = mass1 + mass2

    output = top / bot

    return output


def langevin(alpha, mu):
    """
    Calculates langevin rate constant
    :param alpha: polarizability of neutral reactant (A^3)
    :param mu: reduced mass of reactants (g)
    :return: ion-neutral interaction rate constant (cm^3/s)
    """
    a = 2 * _pi * _e / _sqrt(mu)
    b = _sqrt(alpha)

    output = a * b

    return output


def dipole(mu, mu_d, c, temperature):
    """
    Calculates ion-dipole interaction component to rate constant
    :param mu: reduced mass of reactants (g)
    :param mu_d: dipole moment of neutral reactant (D)
    :param c: unitless factor parameterized by mu_d and alpha
    :param temperature: collision temperature (K)
    :return: ion-dipole interaction rate constant (cm^3/s)
    """
    top = 2 * _sqrt(2) * _pi * _e * c * mu_d
    bot = _sqrt(mu * _pi * _kb * temperature)

    output = top / bot

    return output


def ado(alpha, mu, mu_d, c, temperature):
    """
    Calculates full ADO rate constant
    :param alpha: polarizability of neutral reactant (A^3)
    :param mu: reduced mass of reactants (g)
    :param mu_d: dipole moment of neutral reactant (D)
    :param c: unitless factor parameterized by mu_d and alpha
    :param temperature: collision temperature (K)
    :return: ADO rate constant (cm^3/s)
    """
    output = langevin(mu, alpha) + dipole(mu, mu_d, c, temperature)

    return output


# def quadrupole(self):
#     a = self.langevin
#     b = 0.0146*self.vals[self.species]['Q']**2
#     c = np.sqrt(self.kb * self.temp) * self.e * self.vals[self.species]['alpha']**1.5
#     self.quadrupole = a*b/c


# def aqo(self):
#     k = self.langevin + self.quadrupole
#     self.AQO = k
