import numpy as _np

class FitK:

    def __init__(self, rho1=1, rho2=1):
        """
        Sets densities of H2O, HOD, and D2O
        Default values set to 1
        :param rho1: H2 density
        :param rho2: H2O density
        """
        self.rho1 = rho1
        self.rho2 = rho2

    def Be(self, t, k1, k2, k3, Be0, BeH0, BeOH0):
        output = Be0 * _np.exp(t * (-k1 * self.rho1 - k2 * self.rho2))

        return output

    def BeH(self, t, k1, k2, k3, Be0, BeH0, BeOH0):
        output = (Be0*k1*self.rho1*_np.exp(k3*self.rho2*t + t*(-k1*self.rho1 - k2*self.rho2)) - Be0*k1*self.rho1 - BeH0*k1*self.rho1 - BeH0*k2*self.rho2 + BeH0*k3*self.rho2)*_np.exp(-k3*self.rho2*t)/(-k1*self.rho1 - k2*self.rho2 + k3*self.rho2)

        return output

    def BeOH(self, t, k1, k2, k3, Be0, BeH0, BeOH0):
        output = (-Be0*k1*self.rho1*_np.exp(k3*self.rho2*t) + Be0*k1*self.rho1 - Be0*k2*self.rho2*_np.exp(k3*self.rho2*t) + Be0*k2*self.rho2*_np.exp(k3*self.rho2*t + t*(-k1*self.rho1 - k2*self.rho2)) + Be0*k3*self.rho2*_np.exp(k3*self.rho2*t) - Be0*k3*self.rho2*_np.exp(k3*self.rho2*t + t*(-k1*self.rho1 - k2*self.rho2)) - BeH0*k1*self.rho1*_np.exp(k3*self.rho2*t) + BeH0*k1*self.rho1 - BeH0*k2*self.rho2*_np.exp(k3*self.rho2*t) + BeH0*k2*self.rho2 + BeH0*k3*self.rho2*_np.exp(k3*self.rho2*t) - BeH0*k3*self.rho2 - BeOH0*k1*self.rho1*_np.exp(k3*self.rho2*t) - BeOH0*k2*self.rho2*_np.exp(k3*self.rho2*t) + BeOH0*k3*self.rho2*_np.exp(k3*self.rho2*t))*_np.exp(-k3*self.rho2*t)/(-k1*self.rho1 - k2*self.rho2 + k3*self.rho2)
        
        return output

class FitP:

    def __init__(self, pstate=0):
        """
        Sets P state fraction of Be^+
        :param pstate: P state fraction of Be^+
        """
        self.k1 = 1.2e-9 * pstate
        self.k2 = 2.5e-9 * pstate + 2.2e-9
        self.k3 = 1e-9

    def Be(self, t, rho1, rho2, Be0, BeH0, BeOH0):
        output = Be0 * _np.exp(t * (-self.k1 * rho1 - self.k2 * rho2))

        return output

    def BeH(self, t, rho1, rho2, Be0, BeH0, BeOH0):
        output = (Be0 * self.k1 * rho1 * _np.exp(self.k3 * rho2 * t + t * (
                    -self.k1 * rho1 - self.k2 * rho2)) - Be0 * self.k1 * rho1 - BeH0 * self.k1 * rho1 - BeH0 * self.k2 * rho2 + BeH0 * self.k3 * rho2) * _np.exp(
            -self.k3 * rho2 * t) / (-self.k1 * rho1 - self.k2 * rho2 + self.k3 * rho2)

        return output

    def BeOH(self, t, rho1, rho2, Be0, BeH0, BeOH0):
        output = (-Be0 * self.k1 * rho1 * _np.exp(
            self.k3 * rho2 * t) + Be0 * self.k1 * rho1 - Be0 * self.k2 * rho2 * _np.exp(
            self.k3 * rho2 * t) + Be0 * self.k2 * rho2 * _np.exp(
            self.k3 * rho2 * t + t * (-self.k1 * rho1 - self.k2 * rho2)) + Be0 * self.k3 * rho2 * _np.exp(
            self.k3 * rho2 * t) - Be0 * self.k3 * rho2 * _np.exp(
            self.k3 * rho2 * t + t * (-self.k1 * rho1 - self.k2 * rho2)) - BeH0 * self.k1 * rho1 * _np.exp(
            self.k3 * rho2 * t) + BeH0 * self.k1 * rho1 - BeH0 * self.k2 * rho2 * _np.exp(
            self.k3 * rho2 * t) + BeH0 * self.k2 * rho2 + BeH0 * self.k3 * rho2 * _np.exp(
            self.k3 * rho2 * t) - BeH0 * self.k3 * rho2 - BeOH0 * self.k1 * rho1 * _np.exp(
            self.k3 * rho2 * t) - BeOH0 * self.k2 * rho2 * _np.exp(
            self.k3 * rho2 * t) + BeOH0 * self.k3 * rho2 * _np.exp(self.k3 * rho2 * t)) * _np.exp(
            -self.k3 * rho2 * t) / (-self.k1 * rho1 - self.k2 * rho2 + self.k3 * rho2)

        return output
    
class FitAll:

    def __init__(self):
        pass

    def Be(self, t, k1, k2, k3, rho1, rho2, Be0, BeH0, BeOH0):
        output = Be0 * _np.exp(t * (-k1 * rho1 - k2 * rho2))

        return output

    def BeH(self, t, k1, k2, k3, rho1, rho2, Be0, BeH0, BeOH0):
        output = (Be0 * k1 * rho1 * _np.exp(k3 * rho2 * t + t * (
                    -k1 * rho1 - k2 * rho2)) - Be0 * k1 * rho1 - BeH0 * k1 * rho1 - BeH0 * k2 * rho2 + BeH0 * k3 * rho2) * _np.exp(
            -k3 * rho2 * t) / (-k1 * rho1 - k2 * rho2 + k3 * rho2)

        return output

    def BeOH(self, t, k1, k2, k3, rho1, rho2, Be0, BeH0, BeOH0):
        output = (-Be0 * k1 * rho1 * _np.exp(
            k3 * rho2 * t) + Be0 * k1 * rho1 - Be0 * k2 * rho2 * _np.exp(
            k3 * rho2 * t) + Be0 * k2 * rho2 * _np.exp(
            k3 * rho2 * t + t * (-k1 * rho1 - k2 * rho2)) + Be0 * k3 * rho2 * _np.exp(
            k3 * rho2 * t) - Be0 * k3 * rho2 * _np.exp(
            k3 * rho2 * t + t * (-k1 * rho1 - k2 * rho2)) - BeH0 * k1 * rho1 * _np.exp(
            k3 * rho2 * t) + BeH0 * k1 * rho1 - BeH0 * k2 * rho2 * _np.exp(
            k3 * rho2 * t) + BeH0 * k2 * rho2 + BeH0 * k3 * rho2 * _np.exp(
            k3 * rho2 * t) - BeH0 * k3 * rho2 - BeOH0 * k1 * rho1 * _np.exp(
            k3 * rho2 * t) - BeOH0 * k2 * rho2 * _np.exp(
            k3 * rho2 * t) + BeOH0 * k3 * rho2 * _np.exp(k3 * rho2 * t)) * _np.exp(
            -k3 * rho2 * t) / (-k1 * rho1 - k2 * rho2 + k3 * rho2)

        return output
    