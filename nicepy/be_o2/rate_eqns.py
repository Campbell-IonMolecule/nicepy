import numpy as _np
from nicepy import rate_constants


class Full:
    """
    Equations for Be^+ + O2 reaction
    """

    def __init__(self, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2 density
        :param rho2: H2O density
        :param rho3: O2 density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, k1, k2, k3, k4, k5, Be0, BeH0, BeOH0, O20):

        output = Be0 * (_np.exp(((((-k3 * self.rho3) - (k2 * self.rho2)) - (k1 * self.rho1)) * t)))

        return output

    def BeH(self, t, k1, k2, k3, k4, k5, Be0, BeH0, BeOH0, O20):

        aux0 = (_np.exp((((-k5 * self.rho3) - (k4 * self.rho2)) * t))) * ((Be0 * (k1 * self.rho1)) + (BeH0 * ((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (k4 * self.rho2))))
        output = (aux0 - (Be0 * ((_np.exp(((((-k3 * self.rho3) - (k2 * self.rho2)) - (k1 * self.rho1)) * t))) * (k1 * self.rho1)))) / ((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (k4 * self.rho2))

        return output

    def BeOH(self, t, k1, k2, k3, k4, k5, Be0, BeH0, BeOH0, O20):

        aux0 = ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) * ((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (k4 * self.rho2))
        aux1 = (_np.exp((((-k5 * self.rho3) - (k4 * self.rho2)) * t))) * ((((BeH0 + BeOH0) * (_np.exp((((k4 * self.rho2) + (k5 * self.rho3)) * t)))) - BeH0) * aux0)
        aux2 = (_np.exp((((k1 * self.rho1) + (((k2 + k4) * self.rho2) + ((k3 + k5) * self.rho3))) * t))) * (((k1 * self.rho1) + (k2 * self.rho2)) * ((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (k4 * self.rho2)))
        aux3 = (_np.exp((((k4 * self.rho2) + (k5 * self.rho3)) * t))) * (((k4 - k2) * (self.rho2 * ((k1 * self.rho1) + (k2 * self.rho2)))) + (((k1 * (k5 * self.rho1)) + (k2 * ((k5 - k3) * self.rho2))) * self.rho3))
        aux4 = (_np.exp((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) * t))) * (k1 * (self.rho1 * ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))))
        aux5 = (_np.exp((((((-(k3 + k5) * self.rho3)) - ((k2 + k4) * self.rho2)) - (k1 * self.rho1)) * t))) * ((aux2 + aux3) - aux4)
        output = ((aux1 + (Be0 * aux5)) / ((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (k4 * self.rho2))) / ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))

        return output

    def O2(self, t, k1, k2, k3, k4, k5, Be0, BeH0, BeOH0, O20):
        
        aux0 = (k2 * (O20 * self.rho2)) + (k3 * (((Be0 + O20) - (Be0 * (_np.exp(((((-k3 * self.rho3) - (k2 * self.rho2)) - (k1 * self.rho1)) * t))))) * self.rho3))

        output = ((k1 * (O20 * self.rho1)) + aux0) / ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))

        return output


class PState:
    """
    Equations for Be^+ + O2 reaction
    """

    def __init__(self, pstate=0, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2 density
        :param rho2: H2O density
        :param rho3: O2 density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3
        self.k1 = rate_constants.be_h2.rate_constant(pstate)
        self.k2 = rate_constants.be_h2o.rate_constant(pstate)
        self.k4 = 4.01e-9

    def Be(self, t, k3, k5, Be0, BeH0, BeOH0, O20):
        
        output = Be0 * (_np.exp(((((-k3 * self.rho3) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t)))

        return output

    def BeH(self, t, k3, k5, Be0, BeH0, BeOH0, O20):
        aux0 = (_np.exp((((-k5 * self.rho3) - (self.k4 * self.rho2)) * t))) * ((Be0 * (self.k1 * self.rho1)) + (
        BeH0 * ((((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (self.k4 * self.rho2))))
        output = (aux0 - (
        Be0 * ((_np.exp(((((-k3 * self.rho3) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t))) * (self.k1 * self.rho1)))) / (
                 (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (self.k4 * self.rho2))

        return output

    def BeOH(self, t, k3, k5, Be0, BeH0, BeOH0, O20):
        aux0 = ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) * (
        (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (self.k4 * self.rho2))
        aux1 = (_np.exp((((-k5 * self.rho3) - (self.k4 * self.rho2)) * t))) * (
        (((BeH0 + BeOH0) * (_np.exp((((self.k4 * self.rho2) + (k5 * self.rho3)) * t)))) - BeH0) * aux0)
        aux2 = (_np.exp((((self.k1 * self.rho1) + (((self.k2 + self.k4) * self.rho2) + ((k3 + k5) * self.rho3))) * t))) * (
        ((self.k1 * self.rho1) + (self.k2 * self.rho2)) * (
        (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (self.k4 * self.rho2)))
        aux3 = (_np.exp((((self.k4 * self.rho2) + (k5 * self.rho3)) * t))) * (
        ((self.k4 - self.k2) * (self.rho2 * ((self.k1 * self.rho1) + (self.k2 * self.rho2)))) + (
        ((self.k1 * (k5 * self.rho1)) + (self.k2 * ((k5 - k3) * self.rho2))) * self.rho3))
        aux4 = (_np.exp((((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) * t))) * (
        self.k1 * (self.rho1 * ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))))
        aux5 = (_np.exp((((((-(k3 + k5) * self.rho3)) - ((self.k2 + self.k4) * self.rho2)) - (self.k1 * self.rho1)) * t))) * (
        (aux2 + aux3) - aux4)
        output = ((aux1 + (Be0 * aux5)) / (
        (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (k5 * self.rho3)) - (self.k4 * self.rho2))) / (
                 (self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))

        return output

    def O2(self, t, k3, k5, Be0, BeH0, BeOH0, O20):
        aux0 = (self.k2 * (O20 * self.rho2)) + (k3 * (((Be0 + O20) - (
        Be0 * (_np.exp(((((-k3 * self.rho3) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t))))) * self.rho3))

        output = ((self.k1 * (O20 * self.rho1)) + aux0) / ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))

        return output


class FullNoInterchange:
    """
    Equations for Be^+ + O2 reaction without BeH^+ + O2
    """

    def __init__(self, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2 density
        :param rho2: H2O density
        :param rho3: O2 density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3

    def Be(self, t, k1, k2, k3, k4, Be0, BeH0, BeOH0, O20):
        output = Be0 * (_np.exp(((((-k3 * self.rho3) - (k2 * self.rho2)) - (k1 * self.rho1)) * t)))

        return output

    def BeH(self, t, k1, k2, k3, k4, Be0, BeH0, BeOH0, O20):
        aux0 = Be0 * ((-1. + (_np.exp((((((k4 * self.rho2) - (k3 * self.rho3)) - (k2 * self.rho2)) - (k1 * self.rho1)) * t)))) * (k1 * self.rho1))
        aux1 = (-(_np.exp(((-k4 * (self.rho2 * t))))) * (aux0 - (BeH0 * (((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2)))))
        output = aux1 / (((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2))

        return output

    def BeOH(self, t, k1, k2, k3, k4, Be0, BeH0, BeOH0, O20):
        aux0 = (((BeH0 + BeOH0) * (_np.exp((k4 * (self.rho2 * t))))) - BeH0) * (
                    ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) * (((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2)))
        aux1 = (_np.exp((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) * t))) * (
                    ((k1 * self.rho1) + (k2 * self.rho2)) * (((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2)))
        aux2 = (_np.exp(((((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2)) * t))) * (
                    k1 * (self.rho1 * ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))))
        aux3 = (_np.exp((((((k4 * self.rho2) - (k3 * self.rho3)) - (k2 * self.rho2)) - (k1 * self.rho1)) * t))) * (
                    ((((k4 - k2) * (self.rho2 * ((k1 * self.rho1) + (k2 * self.rho2)))) + aux1) - aux2) - (k2 * (k3 * (self.rho2 * self.rho3))))
        aux4 = ((_np.exp(((-k4 * (self.rho2 * t))))) * (aux0 + (Be0 * aux3))) / (
                    ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3))) - (k4 * self.rho2))
        output = aux4 / ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))

        return output

    def O2(self, t, k1, k2, k3, k4, Be0, BeH0, BeOH0, O20):
        aux0 = (k2 * (O20 * self.rho2)) + (
                    k3 * (((Be0 + O20) - (Be0 * (_np.exp(((((-k3 * self.rho3) - (k2 * self.rho2)) - (k1 * self.rho1)) * t))))) * self.rho3));
        output = ((k1 * (O20 * self.rho1)) + aux0) / ((k1 * self.rho1) + ((k2 * self.rho2) + (k3 * self.rho3)))

        return output


class PstateNoInterchange:
    """
    Equations for Be^+ + O2 reaction without BeH^+ + O2
    """

    def __init__(self, pstate=0, rho1=1, rho2=1, rho3=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2 density
        :param rho2: H2O density
        :param rho3: O2 density
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3
        self.k1 = rate_constants.be_h2.rate_constant(pstate)
        self.k2 = rate_constants.be_h2o.rate_constant(pstate)
        self.k4 = 4.01e-9

    def Be(self, t, k3, Be0, BeH0, BeOH0, O20):
        output = Be0 * (_np.exp(((((-k3 * self.rho3) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t)))

        return output

    def BeH(self, t, k3, Be0, BeH0, BeOH0, O20):
        aux0 = Be0 * ((-1. + (_np.exp((((((self.k4 * self.rho2) - (k3 * self.rho3)) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t)))) * (self.k1 * self.rho1))
        aux1 = (-(_np.exp(((-self.k4 * (self.rho2 * t))))) * (aux0 - (BeH0 * (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2)))))
        output = aux1 / (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2))

        return output

    def BeOH(self, t, k3, Be0, BeH0, BeOH0, O20):
        aux0 = (((BeH0 + BeOH0) * (_np.exp((self.k4 * (self.rho2 * t))))) - BeH0) * (
                    ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) * (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2)))
        aux1 = (_np.exp((((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) * t))) * (
                    ((self.k1 * self.rho1) + (self.k2 * self.rho2)) * (((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2)))
        aux2 = (_np.exp(((((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2)) * t))) * (
                    self.k1 * (self.rho1 * ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))))
        aux3 = (_np.exp((((((self.k4 * self.rho2) - (k3 * self.rho3)) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t))) * (
                    ((((self.k4 - self.k2) * (self.rho2 * ((self.k1 * self.rho1) + (self.k2 * self.rho2)))) + aux1) - aux2) - (self.k2 * (k3 * (self.rho2 * self.rho3))))
        aux4 = ((_np.exp(((-self.k4 * (self.rho2 * t))))) * (aux0 + (Be0 * aux3))) / (
                    ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3))) - (self.k4 * self.rho2))
        output = aux4 / ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))

        return output

    def O2(self, t, k3, Be0, BeH0, BeOH0, O20):
        aux0 = (self.k2 * (O20 * self.rho2)) + (
                    k3 * (((Be0 + O20) - (Be0 * (_np.exp(((((-k3 * self.rho3) - (self.k2 * self.rho2)) - (self.k1 * self.rho1)) * t))))) * self.rho3));
        output = ((self.k1 * (O20 * self.rho1)) + aux0) / ((self.k1 * self.rho1) + ((self.k2 * self.rho2) + (k3 * self.rho3)))

        return output
