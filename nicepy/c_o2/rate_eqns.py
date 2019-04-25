import numpy as _np


class Full:

    def __init__(self):
        """
        Full fitting, no omissions to rate equations
        """
        pass

    def C(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        output = C0 * (_np.exp((((-k2) - k1) * (rho2 * t))))

        return output

    def CO(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (rho2 * t)))) * (
                    ((C0 * (k1 * rho2)) + (CO0 * (((k1 + k2) - k4) * rho2))) - (CO0 * (k3 * rho1)))
        aux1 = (_np.exp(((((-(k1 + (k2 + k4)) * rho2)) - (k3 * rho1)) * t))) * (
                    (C0 * ((_np.exp(((k3 * (rho1 * t)) + (k4 * (rho2 * t))))) * (k1 * rho2))) - aux0)
        output = aux1 / ((k3 * rho1) - (((k1 + k2) - k4) * rho2))

        return output

    def HCO(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (((_np.exp(((k1 + k2) * (rho2 * t)))) * ((k1 + k2) - k4)) + k4) - (
                    (_np.exp((((((k1 + k2) - k4) * rho2) - (k3 * rho1)) * t))) * (k1 + k2))
        aux1 = C0 * (k1 * (k3 * (rho1 * (((1. - (_np.exp(((k1 + k2) * (rho2 * t))))) * (k3 * rho1)) + (aux0 * rho2)))))
        aux2 = ((_np.exp(((k1 + k2) * (rho2 * t)))) * (((CO0 + HCO0) * (k3 * rho1)) + (HCO0 * (k4 * rho2)))) - (
                    CO0 * ((_np.exp((((((k1 + k2) - k4) * rho2) - (k3 * rho1)) * t))) * (k3 * rho1)))
        aux3 = (_np.exp((((-k2) - k1) * (rho2 * t)))) * (
                    aux1 + ((k1 + k2) * (((((k1 + k2) - k4) * rho2) - (k3 * rho1)) * aux2)))
        output = ((aux3 / ((k3 * rho1) + (k4 * rho2))) / ((((k1 + k2) - k4) * rho2) - (k3 * rho1))) / (k1 + k2)

        return output

    def O(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (rho2 * t)))) * ((k5 * (O0 * rho1)) - (((C0 * k2) + (((k1 + k2) - k6) * O0)) * rho2))
        aux1 = (_np.exp(((((-(k1 + (k2 + k6)) * rho2)) - (k5 * rho1)) * t))) * (
                    (C0 * ((_np.exp(((k5 * (rho1 * t)) + (k6 * (rho2 * t))))) * (k2 * rho2))) + aux0)
        output = aux1 / ((k5 * rho1) - (((k1 + k2) - k6) * rho2))

        return output

    def O2(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (k5 * (rho1 - ((_np.exp(((-k7 * (rho2 * t))))) * rho1))) + (
                    (1. - (_np.exp((((-k6 * rho2) - (k5 * rho1)) * t)))) * ((k6 - k7) * rho2))
        aux1 = (((_np.exp((((-k4 * rho2) - (k3 * rho1)) * t))) * (k1 * (k4 * (rho2 ** 2)))) / ((k3 * rho1) + (k4 * rho2))) / (
                    (k3 * rho1) - (((k1 + k2) - k4) * rho2))
        aux2 = ((_np.exp((((-k6 * rho2) - (k5 * rho1)) * t))) * (k2 * ((k6 - k7) * (rho2 ** 2)))) / (
                    (k5 * rho1) + ((k6 - k7) * rho2))
        aux3 = (k2 * (k3 * (k5 * (k7 * (rho1 ** 2))))) + (
                    (k1 + k2) * (((k1 * k4) + ((k2 - k4) * k6)) * (((k1 + k2) - k7) * (rho2 ** 2))))
        aux4 = ((k1 + k2) * ((k1 * (k4 * k5)) + (k2 * (k3 * k6)))) + (
                    (((k1 + k2) * ((k2 - k4) * k5)) - (k2 * (k3 * k6))) * k7)
        aux5 = ((_np.exp((((-k2) - k1) * (rho2 * t)))) * (aux3 - (aux4 * (rho1 * rho2)))) / (
                    (((k1 + k2) - k6) * rho2) - (k5 * rho1))
        aux6 = (aux2 / ((k5 * rho1) - (((k1 + k2) - k6) * rho2))) + (
                    ((aux5 / ((((k1 + k2) - k4) * rho2) - (k3 * rho1))) / ((k1 + k2) - k7)) / (k1 + k2))
        aux7 = (((_np.exp(((-k7 * (rho2 * t))))) * (k2 * (k5 * rho1))) / ((k5 * rho1) + ((k6 - k7) * rho2))) / ((k1 + k2) - k7)
        aux8 = C0 * (((1. + (aux1 + aux6)) - aux7) - (((k1 * (k3 * rho1)) / ((k3 * rho1) + (k4 * rho2))) / (k1 + k2)))
        aux9 = (CO0 * ((-1. + (_np.exp((((-k4 * rho2) - (k3 * rho1)) * t)))) * (k4 * rho2))) / ((k3 * rho1) + (k4 * rho2))
        output = ((O20 + (OH0 + (((O0 * aux0) / ((k5 * rho1) + ((k6 - k7) * rho2))) + aux8))) - aux9) - (
                    (_np.exp(((-k7 * (rho2 * t))))) * OH0)

        return output

    def OH(self, t, k1, k2, k3, k4, k5, k6, k7, rho1, rho2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = ((_np.exp(((k1 + k2) * (rho2 * t)))) * ((k5 * rho1) - (((k1 + k2) - k6) * rho2))) + (
                    (_np.exp((k7 * (rho2 * t)))) * (((k7 * rho2) - (k6 * rho2)) - (k5 * rho1)))
        aux1 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * rho2) - (k5 * rho1)) * t))) * (((k1 + k2) - k7) * rho2)) + aux0
        aux2 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * rho2) - (k5 * rho1)) * t))) * (k5 * (O0 * rho1))) - (
                    (_np.exp(((k1 + k2) * (rho2 * t)))) * ((k5 * ((O0 + OH0) * rho1)) + ((k6 - k7) * (OH0 * rho2))))
        aux3 = (C0 * (k2 * (k5 * (rho1 * aux1)))) + (((k1 + k2) - k7) * (((((k1 + k2) - k6) * rho2) - (k5 * rho1)) * aux2))
        aux4 = (((_np.exp(((((-k7) - k2) - k1) * (rho2 * t)))) * aux3) / (((k7 - k6) * rho2) - (k5 * rho1))) / (
                    (((k1 + k2) - k6) * rho2) - (k5 * rho1))
        output = aux4 / ((k1 + k2) - k7)
        
        return output


class FixK:

    def __init__(self, k3=7.5e-10, k4=1.2e-10, k5=1.7e-9, k6=1.9e-10, k7=5.9e-10, rho1=1, rho2=1):
        """
        Fixes k3-7 and pressures
        :param k3: 
        :param k4: 
        :param k5: 
        :param k6: 
        :param k7: 
        :param rho1: 
        :param rho2: 
        """
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.k7 = k7
        self.rho1 = rho1
        self.rho2 = rho2

    def C(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        output = C0 * (_np.exp((((-k2) - k1) * (self.rho2 * t))))

        return output

    def CO(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                ((C0 * (k1 * self.rho2)) + (CO0 * (((k1 + k2) - self.k4) * self.rho2))) - (CO0 * (self.k3 * self.rho1)))
        aux1 = (_np.exp(((((-(k1 + (k2 + self.k4)) * self.rho2)) - (self.k3 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((self.k3 * (self.rho1 * t)) + (self.k4 * (self.rho2 * t))))) * (k1 * self.rho2))) - aux0)
        output = aux1 / ((self.k3 * self.rho1) - (((k1 + k2) - self.k4) * self.rho2))

        return output

    def HCO(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k1 + k2) - self.k4)) + self.k4) - (
                (_np.exp((((((k1 + k2) - self.k4) * self.rho2) - (self.k3 * self.rho1)) * t))) * (k1 + k2))
        aux1 = C0 * (k1 * (self.k3 * (self.rho1 * (((1. - (_np.exp(((k1 + k2) * (self.rho2 * t))))) * (self.k3 * self.rho1)) + (aux0 * self.rho2)))))
        aux2 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * (((CO0 + HCO0) * (self.k3 * self.rho1)) + (HCO0 * (self.k4 * self.rho2)))) - (
                CO0 * ((_np.exp((((((k1 + k2) - self.k4) * self.rho2) - (self.k3 * self.rho1)) * t))) * (self.k3 * self.rho1)))
        aux3 = (_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (
                aux1 + ((k1 + k2) * (((((k1 + k2) - self.k4) * self.rho2) - (self.k3 * self.rho1)) * aux2)))
        output = ((aux3 / ((self.k3 * self.rho1) + (self.k4 * self.rho2))) / ((((k1 + k2) - self.k4) * self.rho2) - (self.k3 * self.rho1))) / (k1 + k2)

        return output

    def O(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                    (self.k5 * (O0 * self.rho1)) - (((C0 * k2) + (((k1 + k2) - self.k6) * O0)) * self.rho2))
        aux1 = (_np.exp(((((-(k1 + (k2 + self.k6)) * self.rho2)) - (self.k5 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((self.k5 * (self.rho1 * t)) + (self.k6 * (self.rho2 * t))))) * (k2 * self.rho2))) + aux0)
        output = aux1 / ((self.k5 * self.rho1) - (((k1 + k2) - self.k6) * self.rho2))

        return output

    def O2(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (self.k5 * (self.rho1 - ((_np.exp(((-self.k7 * (self.rho2 * t))))) * self.rho1))) + (
                (1. - (_np.exp((((-self.k6 * self.rho2) - (self.k5 * self.rho1)) * t)))) * ((self.k6 - self.k7) * self.rho2))
        aux1 = (((_np.exp((((-self.k4 * self.rho2) - (self.k3 * self.rho1)) * t))) * (k1 * (self.k4 * (self.rho2 ** 2)))) / (
                    (self.k3 * self.rho1) + (self.k4 * self.rho2))) / (
                       (self.k3 * self.rho1) - (((k1 + k2) - self.k4) * self.rho2))
        aux2 = ((_np.exp((((-self.k6 * self.rho2) - (self.k5 * self.rho1)) * t))) * (k2 * ((self.k6 - self.k7) * (self.rho2 ** 2)))) / (
                (self.k5 * self.rho1) + ((self.k6 - self.k7) * self.rho2))
        aux3 = (k2 * (self.k3 * (self.k5 * (self.k7 * (self.rho1 ** 2))))) + (
                (k1 + k2) * (((k1 * self.k4) + ((k2 - self.k4) * self.k6)) * (((k1 + k2) - self.k7) * (self.rho2 ** 2))))
        aux4 = ((k1 + k2) * ((k1 * (self.k4 * self.k5)) + (k2 * (self.k3 * self.k6)))) + (
                (((k1 + k2) * ((k2 - self.k4) * self.k5)) - (k2 * (self.k3 * self.k6))) * self.k7)
        aux5 = ((_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (aux3 - (aux4 * (self.rho1 * self.rho2)))) / (
                (((k1 + k2) - self.k6) * self.rho2) - (self.k5 * self.rho1))
        aux6 = (aux2 / ((self.k5 * self.rho1) - (((k1 + k2) - self.k6) * self.rho2))) + (
                ((aux5 / ((((k1 + k2) - self.k4) * self.rho2) - (self.k3 * self.rho1))) / ((k1 + k2) - self.k7)) / (k1 + k2))
        aux7 = (((_np.exp(((-self.k7 * (self.rho2 * t))))) * (k2 * (self.k5 * self.rho1))) / ((self.k5 * self.rho1) + ((self.k6 - self.k7) * self.rho2))) / (
                    (k1 + k2) - self.k7)
        aux8 = C0 * (((1. + (aux1 + aux6)) - aux7) - (((k1 * (self.k3 * self.rho1)) / ((self.k3 * self.rho1) + (self.k4 * self.rho2))) / (k1 + k2)))
        aux9 = (CO0 * ((-1. + (_np.exp((((-self.k4 * self.rho2) - (self.k3 * self.rho1)) * t)))) * (self.k4 * self.rho2))) / (
                    (self.k3 * self.rho1) + (self.k4 * self.rho2))
        output = ((O20 + (OH0 + (((O0 * aux0) / ((self.k5 * self.rho1) + ((self.k6 - self.k7) * self.rho2))) + aux8))) - aux9) - (
                (_np.exp(((-self.k7 * (self.rho2 * t))))) * OH0)

        return output

    def OH(self, t, k1, k2, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((self.k5 * self.rho1) - (((k1 + k2) - self.k6) * self.rho2))) + (
                (_np.exp((self.k7 * (self.rho2 * t)))) * (((self.k7 * self.rho2) - (self.k6 * self.rho2)) - (self.k5 * self.rho1)))
        aux1 = ((_np.exp((((((k1 + (k2 + self.k7)) - self.k6) * self.rho2) - (self.k5 * self.rho1)) * t))) * (((k1 + k2) - self.k7) * self.rho2)) + aux0
        aux2 = ((_np.exp((((((k1 + (k2 + self.k7)) - self.k6) * self.rho2) - (self.k5 * self.rho1)) * t))) * (self.k5 * (O0 * self.rho1))) - (
                (_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((self.k5 * ((O0 + OH0) * self.rho1)) + ((self.k6 - self.k7) * (OH0 * self.rho2))))
        aux3 = (C0 * (k2 * (self.k5 * (self.rho1 * aux1)))) + (
                    ((k1 + k2) - self.k7) * (((((k1 + k2) - self.k6) * self.rho2) - (self.k5 * self.rho1)) * aux2))
        aux4 = (((_np.exp(((((-self.k7) - k2) - k1) * (self.rho2 * t)))) * aux3) / (((self.k7 - self.k6) * self.rho2) - (self.k5 * self.rho1))) / (
                (((k1 + k2) - self.k6) * self.rho2) - (self.k5 * self.rho1))
        output = aux4 / ((k1 + k2) - self.k7)

        return output


class FixP:

    def __init__(self, rho1=1, rho2=1):
        """
        Fixes densities, set to 1 to fit for reaction rate Gamma[1/s]
        :param rho1: 
        :param rho2: 
        """
        self.rho1 = rho1
        self.rho2 = rho2

    def C(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        output = C0 * (_np.exp((((-k2) - k1) * (self.rho2 * t))))

        return output

    def CO(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                ((C0 * (k1 * self.rho2)) + (CO0 * (((k1 + k2) - k4) * self.rho2))) - (CO0 * (k3 * self.rho1)))
        aux1 = (_np.exp(((((-(k1 + (k2 + k4)) * self.rho2)) - (k3 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k3 * (self.rho1 * t)) + (k4 * (self.rho2 * t))))) * (k1 * self.rho2))) - aux0)
        output = aux1 / ((k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))

        return output

    def HCO(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k1 + k2) - k4)) + k4) - (
                (_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k1 + k2))
        aux1 = C0 * (k1 * (k3 * (self.rho1 * (((1. - (_np.exp(((k1 + k2) * (self.rho2 * t))))) * (k3 * self.rho1)) + (aux0 * self.rho2)))))
        aux2 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * (((CO0 + HCO0) * (k3 * self.rho1)) + (HCO0 * (k4 * self.rho2)))) - (
                CO0 * ((_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k3 * self.rho1)))
        aux3 = (_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (
                aux1 + ((k1 + k2) * (((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * aux2)))
        output = ((aux3 / ((k3 * self.rho1) + (k4 * self.rho2))) / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / (k1 + k2)

        return output

    def O(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                    (k5 * (O0 * self.rho1)) - (((C0 * k2) + (((k1 + k2) - k6) * O0)) * self.rho2))
        aux1 = (_np.exp(((((-(k1 + (k2 + k6)) * self.rho2)) - (k5 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k5 * (self.rho1 * t)) + (k6 * (self.rho2 * t))))) * (k2 * self.rho2))) + aux0)
        output = aux1 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))

        return output

    def O2(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = (k5 * (self.rho1 - ((_np.exp(((-k7 * (self.rho2 * t))))) * self.rho1))) + (
                (1. - (_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t)))) * ((k6 - k7) * self.rho2))
        aux1 = (((_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t))) * (k1 * (k4 * (self.rho2 ** 2)))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))) / (
                       (k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))
        aux2 = ((_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t))) * (k2 * ((k6 - k7) * (self.rho2 ** 2)))) / (
                (k5 * self.rho1) + ((k6 - k7) * self.rho2))
        aux3 = (k2 * (k3 * (k5 * (k7 * (self.rho1 ** 2))))) + (
                (k1 + k2) * (((k1 * k4) + ((k2 - k4) * k6)) * (((k1 + k2) - k7) * (self.rho2 ** 2))))
        aux4 = ((k1 + k2) * ((k1 * (k4 * k5)) + (k2 * (k3 * k6)))) + (
                (((k1 + k2) * ((k2 - k4) * k5)) - (k2 * (k3 * k6))) * k7)
        aux5 = ((_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (aux3 - (aux4 * (self.rho1 * self.rho2)))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        aux6 = (aux2 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                ((aux5 / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / ((k1 + k2) - k7)) / (k1 + k2))
        aux7 = (((_np.exp(((-k7 * (self.rho2 * t))))) * (k2 * (k5 * self.rho1))) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) / (
                    (k1 + k2) - k7)
        aux8 = C0 * (((1. + (aux1 + aux6)) - aux7) - (((k1 * (k3 * self.rho1)) / ((k3 * self.rho1) + (k4 * self.rho2))) / (k1 + k2)))
        aux9 = (CO0 * ((-1. + (_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t)))) * (k4 * self.rho2))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))
        output = ((O20 + (OH0 + (((O0 * aux0) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) + aux8))) - aux9) - (
                (_np.exp(((-k7 * (self.rho2 * t))))) * OH0)

        return output

    def OH(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0, OH0, O20, HCO0):
        aux0 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                (_np.exp((k7 * (self.rho2 * t)))) * (((k7 * self.rho2) - (k6 * self.rho2)) - (k5 * self.rho1)))
        aux1 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (((k1 + k2) - k7) * self.rho2)) + aux0
        aux2 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (k5 * (O0 * self.rho1))) - (
                (_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * ((O0 + OH0) * self.rho1)) + ((k6 - k7) * (OH0 * self.rho2))))
        aux3 = (C0 * (k2 * (k5 * (self.rho1 * aux1)))) + (
                    ((k1 + k2) - k7) * (((((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1)) * aux2))
        aux4 = (((_np.exp(((((-k7) - k2) - k1) * (self.rho2 * t)))) * aux3) / (((k7 - k6) * self.rho2) - (k5 * self.rho1))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        output = aux4 / ((k1 + k2) - k7)

        return output
    
    
class FixPOffset:

    def __init__(self, rho1=1, rho2=1, CO0=0, O0=0, OH0=0, O20=0, HCO0=0):
        """
        Fixes densities, set to 1 to fit for reaction rate Gamma[1/s]
        :param rho1: 
        :param rho2: 
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.CO0 = CO0
        self.O0 = O0
        self.OH0 = OH0
        self.O20 = O20
        self.HCO0 = HCO0

    def C(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        output = C0 * (_np.exp((((-k2) - k1) * (self.rho2 * t))))

        return output

    def CO(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                ((C0 * (k1 * self.rho2)) + (self.CO0 * (((k1 + k2) - k4) * self.rho2))) - (self.CO0 * (k3 * self.rho1)))
        aux1 = (_np.exp(((((-(k1 + (k2 + k4)) * self.rho2)) - (k3 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k3 * (self.rho1 * t)) + (k4 * (self.rho2 * t))))) * (k1 * self.rho2))) - aux0)
        output = aux1 / ((k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))

        return output

    def HCO(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        aux0 = (((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k1 + k2) - k4)) + k4) - (
                (_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k1 + k2))
        aux1 = C0 * (k1 * (k3 * (self.rho1 * (((1. - (_np.exp(((k1 + k2) * (self.rho2 * t))))) * (k3 * self.rho1)) + (aux0 * self.rho2)))))
        aux2 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * (((self.CO0 + self.HCO0) * (k3 * self.rho1)) + (self.HCO0 * (k4 * self.rho2)))) - (
                self.CO0 * ((_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k3 * self.rho1)))
        aux3 = (_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (
                aux1 + ((k1 + k2) * (((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * aux2)))
        output = ((aux3 / ((k3 * self.rho1) + (k4 * self.rho2))) / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / (k1 + k2)

        return output

    def O(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                    (k5 * (self.O0 * self.rho1)) - (((C0 * k2) + (((k1 + k2) - k6) * self.O0)) * self.rho2))
        aux1 = (_np.exp(((((-(k1 + (k2 + k6)) * self.rho2)) - (k5 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k5 * (self.rho1 * t)) + (k6 * (self.rho2 * t))))) * (k2 * self.rho2))) + aux0)
        output = aux1 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))

        return output

    def O2(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        aux0 = (k5 * (self.rho1 - ((_np.exp(((-k7 * (self.rho2 * t))))) * self.rho1))) + (
                (1. - (_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t)))) * ((k6 - k7) * self.rho2))
        aux1 = (((_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t))) * (k1 * (k4 * (self.rho2 ** 2)))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))) / (
                       (k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))
        aux2 = ((_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t))) * (k2 * ((k6 - k7) * (self.rho2 ** 2)))) / (
                (k5 * self.rho1) + ((k6 - k7) * self.rho2))
        aux3 = (k2 * (k3 * (k5 * (k7 * (self.rho1 ** 2))))) + (
                (k1 + k2) * (((k1 * k4) + ((k2 - k4) * k6)) * (((k1 + k2) - k7) * (self.rho2 ** 2))))
        aux4 = ((k1 + k2) * ((k1 * (k4 * k5)) + (k2 * (k3 * k6)))) + (
                (((k1 + k2) * ((k2 - k4) * k5)) - (k2 * (k3 * k6))) * k7)
        aux5 = ((_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (aux3 - (aux4 * (self.rho1 * self.rho2)))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        aux6 = (aux2 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                ((aux5 / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / ((k1 + k2) - k7)) / (k1 + k2))
        aux7 = (((_np.exp(((-k7 * (self.rho2 * t))))) * (k2 * (k5 * self.rho1))) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) / (
                    (k1 + k2) - k7)
        aux8 = C0 * (((1. + (aux1 + aux6)) - aux7) - (((k1 * (k3 * self.rho1)) / ((k3 * self.rho1) + (k4 * self.rho2))) / (k1 + k2)))
        aux9 = (self.CO0 * ((-1. + (_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t)))) * (k4 * self.rho2))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))
        output = ((self.O20 + (self.OH0 + (((self.O0 * aux0) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) + aux8))) - aux9) - (
                (_np.exp(((-k7 * (self.rho2 * t))))) * self.OH0)

        return output

    def OH(self, t, k1, k2, k3, k4, k5, k6, k7, C0):
        aux0 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                (_np.exp((k7 * (self.rho2 * t)))) * (((k7 * self.rho2) - (k6 * self.rho2)) - (k5 * self.rho1)))
        aux1 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (((k1 + k2) - k7) * self.rho2)) + aux0
        aux2 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (k5 * (self.O0 * self.rho1))) - (
                (_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * ((self.O0 + self.OH0) * self.rho1)) + ((k6 - k7) * (self.OH0 * self.rho2))))
        aux3 = (C0 * (k2 * (k5 * (self.rho1 * aux1)))) + (
                    ((k1 + k2) - k7) * (((((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1)) * aux2))
        aux4 = (((_np.exp(((((-k7) - k2) - k1) * (self.rho2 * t)))) * aux3) / (((k7 - k6) * self.rho2) - (k5 * self.rho1))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        output = aux4 / ((k1 + k2) - k7)

        return output
    

class FixPOffset2:

    def __init__(self, rho1=1, rho2=1, OH0=0, O20=0, HCO0=0):
        """
        Fixes densities, set to 1 to fit for reaction rate Gamma[1/s]
        :param rho1: 
        :param rho2: 
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.OH0 = OH0
        self.O20 = O20
        self.HCO0 = HCO0

    def C(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        output = C0 * (_np.exp((((-k2) - k1) * (self.rho2 * t))))

        return output

    def CO(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                ((C0 * (k1 * self.rho2)) + (CO0 * (((k1 + k2) - k4) * self.rho2))) - (CO0 * (k3 * self.rho1)))
        aux1 = (_np.exp(((((-(k1 + (k2 + k4)) * self.rho2)) - (k3 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k3 * (self.rho1 * t)) + (k4 * (self.rho2 * t))))) * (k1 * self.rho2))) - aux0)
        output = aux1 / ((k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))

        return output

    def HCO(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        aux0 = (((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k1 + k2) - k4)) + k4) - (
                (_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k1 + k2))
        aux1 = C0 * (k1 * (k3 * (self.rho1 * (((1. - (_np.exp(((k1 + k2) * (self.rho2 * t))))) * (k3 * self.rho1)) + (aux0 * self.rho2)))))
        aux2 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * (((CO0 + self.HCO0) * (k3 * self.rho1)) + (self.HCO0 * (k4 * self.rho2)))) - (
                CO0 * ((_np.exp((((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * t))) * (k3 * self.rho1)))
        aux3 = (_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (
                aux1 + ((k1 + k2) * (((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1)) * aux2)))
        output = ((aux3 / ((k3 * self.rho1) + (k4 * self.rho2))) / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / (k1 + k2)

        return output

    def O(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        aux0 = (_np.exp(((k1 + k2) * (self.rho2 * t)))) * (
                    (k5 * (O0 * self.rho1)) - (((C0 * k2) + (((k1 + k2) - k6) * O0)) * self.rho2))
        aux1 = (_np.exp(((((-(k1 + (k2 + k6)) * self.rho2)) - (k5 * self.rho1)) * t))) * (
                (C0 * ((_np.exp(((k5 * (self.rho1 * t)) + (k6 * (self.rho2 * t))))) * (k2 * self.rho2))) + aux0)
        output = aux1 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))

        return output

    def O2(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        aux0 = (k5 * (self.rho1 - ((_np.exp(((-k7 * (self.rho2 * t))))) * self.rho1))) + (
                (1. - (_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t)))) * ((k6 - k7) * self.rho2))
        aux1 = (((_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t))) * (k1 * (k4 * (self.rho2 ** 2)))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))) / (
                       (k3 * self.rho1) - (((k1 + k2) - k4) * self.rho2))
        aux2 = ((_np.exp((((-k6 * self.rho2) - (k5 * self.rho1)) * t))) * (k2 * ((k6 - k7) * (self.rho2 ** 2)))) / (
                (k5 * self.rho1) + ((k6 - k7) * self.rho2))
        aux3 = (k2 * (k3 * (k5 * (k7 * (self.rho1 ** 2))))) + (
                (k1 + k2) * (((k1 * k4) + ((k2 - k4) * k6)) * (((k1 + k2) - k7) * (self.rho2 ** 2))))
        aux4 = ((k1 + k2) * ((k1 * (k4 * k5)) + (k2 * (k3 * k6)))) + (
                (((k1 + k2) * ((k2 - k4) * k5)) - (k2 * (k3 * k6))) * k7)
        aux5 = ((_np.exp((((-k2) - k1) * (self.rho2 * t)))) * (aux3 - (aux4 * (self.rho1 * self.rho2)))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        aux6 = (aux2 / ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                ((aux5 / ((((k1 + k2) - k4) * self.rho2) - (k3 * self.rho1))) / ((k1 + k2) - k7)) / (k1 + k2))
        aux7 = (((_np.exp(((-k7 * (self.rho2 * t))))) * (k2 * (k5 * self.rho1))) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) / (
                    (k1 + k2) - k7)
        aux8 = C0 * (((1. + (aux1 + aux6)) - aux7) - (((k1 * (k3 * self.rho1)) / ((k3 * self.rho1) + (k4 * self.rho2))) / (k1 + k2)))
        aux9 = (CO0 * ((-1. + (_np.exp((((-k4 * self.rho2) - (k3 * self.rho1)) * t)))) * (k4 * self.rho2))) / (
                    (k3 * self.rho1) + (k4 * self.rho2))
        output = ((self.O20 + (self.OH0 + (((O0 * aux0) / ((k5 * self.rho1) + ((k6 - k7) * self.rho2))) + aux8))) - aux9) - (
                (_np.exp(((-k7 * (self.rho2 * t))))) * self.OH0)

        return output

    def OH(self, t, k1, k2, k3, k4, k5, k6, k7, C0, CO0, O0):
        aux0 = ((_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * self.rho1) - (((k1 + k2) - k6) * self.rho2))) + (
                (_np.exp((k7 * (self.rho2 * t)))) * (((k7 * self.rho2) - (k6 * self.rho2)) - (k5 * self.rho1)))
        aux1 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (((k1 + k2) - k7) * self.rho2)) + aux0
        aux2 = ((_np.exp((((((k1 + (k2 + k7)) - k6) * self.rho2) - (k5 * self.rho1)) * t))) * (k5 * (O0 * self.rho1))) - (
                (_np.exp(((k1 + k2) * (self.rho2 * t)))) * ((k5 * ((O0 + self.OH0) * self.rho1)) + ((k6 - k7) * (self.OH0 * self.rho2))))
        aux3 = (C0 * (k2 * (k5 * (self.rho1 * aux1)))) + (
                    ((k1 + k2) - k7) * (((((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1)) * aux2))
        aux4 = (((_np.exp(((((-k7) - k2) - k1) * (self.rho2 * t)))) * aux3) / (((k7 - k6) * self.rho2) - (k5 * self.rho1))) / (
                (((k1 + k2) - k6) * self.rho2) - (k5 * self.rho1))
        output = aux4 / ((k1 + k2) - k7)

        return output
