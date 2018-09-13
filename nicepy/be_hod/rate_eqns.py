import numpy as _np

# Equations for Be^+ + HOD reaction without BeOH and BeOD interchange


class Direct:

    def __init__(self, rho1=1, rho2=1, rho3=1):
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
