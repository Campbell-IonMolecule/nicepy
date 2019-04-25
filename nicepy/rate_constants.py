from numpy import sqrt as _sqrt, pi as _pi
_e = 4.8e-10  # CGS elementary charge
_kb = 1.38e-16  # CGS boltzmann constant


def be_h2(pstate):
    output = 1.2e-9 * pstate

    return output


def be_h2o(pstate):
    output = 2.5e-9 * pstate + 2.2e-9

    return output


def be_o2(pstate):
    output = 5.5e-11 + pstate

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


def langevin(alpha, mass1, mass2):
    """
    :param alpha: alpha: polarizability of neutral reactant (A^3)
    :param mass1: mass of first particle (amu)
    :param mass2: mass of second particle (amu)
    :return: ion-neutral interaction rate constant (cm^3/s)
    """
    mu = reduced_mass(mass1, mass2) * 1.660538782e-24
    alpha = alpha * 1e-24
    top = 2 * _pi * _e * _sqrt(alpha)
    bot = _sqrt(mu)

    output = top / bot

    return output


def dipole(temperature, mass1, mass2, mu_d, c):
    """
    Calculates ion-dipole interaction component to rate constant
    :param temperature: collision temperature (K)
    :param mass1: mass of first particle (amu)
    :param mass2: mass of second particle (amu)
    :param mu_d: dipole moment of neutral reactant (D)
    :param c: unitless factor parameterized by mu_d and alpha
    :return: ion-dipole interaction rate constant (cm^3/s)
    """
    mu = reduced_mass(mass1, mass2) * 1.660538782e-24
    mu_d = mu_d *1e-18
    top = 2 * _sqrt(2) * _pi * _e * c * mu_d
    bot = _sqrt(mu * _pi * _kb * temperature)

    output = top / bot

    return output


def ado(temperature, mass1, mass2, mu_d, c, alpha):
    """
    Calculates full ADO rate constant
    :param alpha: polarizability of neutral reactant (A^3)
    :param mass1: mass of first particle (amu)
    :param mass2: mass of second particle (amu)
    :param mu_d: dipole moment of neutral reactant (D)
    :param c: unitless factor parameterized by mu_d and alpha
    :param temperature: collision temperature (K)
    :return: ADO rate constant (cm^3/s)
    """
    output = langevin(alpha=alpha, mass1=mass1, mass2=mass2) + dipole(temperature=temperature, mass1=mass2, mass2=mass2, mu_d=mu_d, c=c)

    return output


def aqo(temperature, mass1, mass2, Q, alpha):
    """

    :param temperature:
    :param mass1:
    :param mass2:
    :param Q:
    :param alpha:
    :return:
    """
    _e = 4.8e-10  # CGS elementary charge
    _kb = 1.38e-16  # CGS boltzmann constant
    kL = langevin(alpha=alpha, mass1=mass1, mass2=mass2)

    Q = Q * 1e-18 * 1e-8
    alpha = alpha * 1e-24

    top = Q ** 2
    bot = (_kb * temperature) ** 0.5 * _e * alpha ** 1.5

    output = kL * (1 + 0.0146 * top / bot)

    return output