import numpy as _np


class Direct:
    """
    Equations for Be^+ + HOD reaction without BeOH and BeOD interchange
    """

    def __init__(self, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2O density
        :param rho2: HOD density
        :param rho3: D2O density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, k1, k2, k3, eta, Be0, BeOH0, BeOD0):

        output = Be0 * (_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t)))

        return output
    
    def BeOH(self, t, k1, k2, k3, eta, Be0, BeOH0, BeOD0):
        aux0 = (_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t))) * ((-1. + (_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t)))) * ((k1 * self.rho1 + k2 * self.rho2) - (eta * k2 * self.rho2)))

        output = ((Be0 * aux0) + (BeOH0 * (k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)))) / (k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3))

        return output
    
    def BeOD(self, t, k1, k2, k3, eta, Be0, BeOH0, BeOD0):
        aux0 = (BeOD0 * ((_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t))) * (k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)))) + (Be0 * ((-1. + (_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t)))) * ((eta * k2 * self.rho2) + k3 * self.rho3)))

        output = ((_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t))) * aux0) / (k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3))

        return output

class FixK:
    """
    Equations for Be^+ + HOD reaction without BeOH and BeOD interchange but fixing the rate constant for Be^+ + H2O (k1)
    """

    def __init__(self, rho1=1, rho2=1, rho3=1):
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, k1, eta, Be0, BeOH0, BeOD0):

        output = Be0 * (_np.exp(((((-k1 * k1 * self.rho3) - k1 * k1 * self.rho2) - k1 * k1 * self.rho1) * t)))

        return output

    def BeOH(self, t, k1, eta, Be0, BeOH0, BeOD0):
        aux0 = (_np.exp(((((-k1 * k1 * self.rho3) - k1 * k1 * self.rho2) - k1 * k1 * self.rho1) * t))) * ((-1. + (_np.exp(((k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3)) * t)))) * ((k1 * k1 * self.rho1 + k1 * k1 * self.rho2) - (eta * k1 * k1 * self.rho2)))

        output = ((Be0 * aux0) + (BeOH0 * (k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3)))) / (k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3))

        return output

    def BeOD(self, t, k1, eta, Be0, BeOH0, BeOD0):
        aux0 = (BeOD0 * ((_np.exp(((k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3)) * t))) * (k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3)))) + (Be0 * ((-1. + (_np.exp(((k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3)) * t)))) * ((eta * k1 * k1 * self.rho2) + k1 * k1 * self.rho3)))

        output = ((_np.exp(((((-k1 * k1 * self.rho3) - k1 * k1 * self.rho2) - k1 * k1 * self.rho1) * t))) * aux0) / (k1 * k1 * self.rho1 + (k1 * k1 * self.rho2 + k1 * k1 * self.rho3))

        return output

class FixK1:
    """
    Equations for Be^+ + HOD reaction without BeOH and BeOD interchange but fixing the rate constant for Be^+ + H2O (k1)
    """

    def __init__(self, k1=2.2e-9, rho1=1, rho2=1, rho3=1):
        self.k1 = k1
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, alpha, k2, k3, eta, Be0, BeOH0, BeOD0):

        output = Be0 * (_np.exp(((((-k3 * alpha * self.rho3) - k2 * alpha * self.rho2) - self.k1 * alpha * self.rho1) * t)))

        return output

    def BeOH(self, t, alpha, k2, k3, eta, Be0, BeOH0, BeOD0):
        aux0 = (_np.exp(((((-k3 * alpha * self.rho3) - k2 * alpha * self.rho2) - self.k1 * alpha * self.rho1) * t))) * ((-1. + (_np.exp(((self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3)) * t)))) * ((self.k1 * alpha * self.rho1 + k2 * alpha * self.rho2) - (eta * k2 * alpha * self.rho2)))

        output = ((Be0 * aux0) + (BeOH0 * (self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3)))) / (self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3))

        return output

    def BeOD(self, t, alpha, k2, k3, eta, Be0, BeOH0, BeOD0):
        aux0 = (BeOD0 * ((_np.exp(((self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3)) * t))) * (self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3)))) + (Be0 * ((-1. + (_np.exp(((self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3)) * t)))) * ((eta * k2 * alpha * self.rho2) + k3 * alpha * self.rho3)))

        output = ((_np.exp(((((-k3 * alpha * self.rho3) - k2 * alpha * self.rho2) - self.k1 * alpha * self.rho1) * t))) * aux0) / (self.k1 * alpha * self.rho1 + (k2 * alpha * self.rho2 + k3 * alpha * self.rho3))

        return output

class Separate:

    def __init__(self, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2O density
        :param rho2: HOD density
        :param rho3: D2O density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, k1, k2, k3, Be0):
        output = Be0 * (_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t)))

        return output

    def BeOH(self, t, k1, k2, k3, eta, Be0, BeOH0):
        aux0 = (_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t))) * (
                    (-1. + (_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t)))) * (
                        (k1 * self.rho1 + k2 * self.rho2) - (eta * k2 * self.rho2)))

        output = ((Be0 * aux0) + (BeOH0 * (k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)))) / (
                    k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3))

        return output

    def BeOD(self, t, k1, k2, k3, eta, Be0, BeOD0):
        aux0 = (BeOD0 * ((_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t))) * (
                    k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)))) + (Be0 * (
                    (-1. + (_np.exp(((k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3)) * t)))) * (
                        (eta * k2 * self.rho2) + k3 * self.rho3)))

        output = ((_np.exp(((((-k3 * self.rho3) - k2 * self.rho2) - k1 * self.rho1) * t))) * aux0) / (
                    k1 * self.rho1 + (k2 * self.rho2 + k3 * self.rho3))

        return output