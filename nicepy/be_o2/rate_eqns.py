import numpy as _np


class Direct:
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