import numpy as _np
from nicepy import rate_constants


class Full:
    """
    Equations for C^+ + H2O reaction
    """

    def __init__(self, rho=1):
        """
        Sets density of H2O
        Default values set to 1
        :param rho1: H2O density
        """
        self.rho = rho

    def C(self, t, k1, k2, C0, mz290, H3O0):
        t = _np.array(t)
        output = -C0 * (-k1 + k2) * _np.exp(-k1 * self.rho * t)/(k1 - k2)
        return output

    def H3O(self, t, k1, k2, C0, mz290, H3O0):
        t = _np.array(t)
        output = C0 * k2 * _np.exp(-k1 * self.rho * t)/(k1 - k2) + C0 + H3O0 + mz290 - (C0 * k1 + mz290 * (k1 - k2)) * _np.exp(-k2 * self.rho * t)/(k1 - k2)
        return output

    def HCO(self, t, k1, k2, C0, mz290, H3O0):
        t = _np.array(t)
        output = -C0 * k1 * _np.exp(-k1 * self.rho * t)/(k1 - k2) + (C0 * k1 + mz290 * (k1 - k2)) * _np.exp(-k2 * self.rho * t)/(k1 - k2)
        return output


class Scaled:
    """
    Equations for C^+ + H2O reaction
    """

    def __init__(self, rho=1):
        """
        Sets density of H2O
        Default values set to 1
        :param rho1: H2O density
        """
        self.rho = rho

    def C(self, t, k1, C0, mz290, H3O0):
        t = _np.array(t)
        output = -C0 * (-k1 + k1*1.004) * _np.exp(-k1 * self.rho * t)/(k1 - k1*1.004)
        return output

    def H3O(self, t, k1, C0, mz290, H3O0):
        t = _np.array(t)
        output = C0 * k1*1.004 * _np.exp(-k1 * self.rho * t)/(k1 - k1*1.004) + C0 + H3O0 + mz290 - (C0 * k1 + mz290 * (k1 - k1*1.004)) * _np.exp(-k1*1.004 * self.rho * t)/(k1 - k1*1.004)
        return output

    def HCO(self, t, k1, C0, mz290, H3O0):
        t = _np.array(t)
        output = -C0 * k1 * _np.exp(-k1 * self.rho * t)/(k1 - k1*1.004) + (C0 * k1 + mz290 * (k1 - k1*1.004)) * _np.exp(-k1*1.004 * self.rho * t)/(k1 - k1*1.004)
        return output