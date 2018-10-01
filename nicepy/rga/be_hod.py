import numpy as _np


class HOD:
    
    def __init__(self, k1=2.03e-9, k2=2.18e-9, k3=2.29e-9, sk1=0.04e-9, sk2=0.07e-9, sk3=0.05e-9):
        self.k1 = k1  # Be + H2O
        self.k2 = k2  # Be + HOD
        self.k3 = k3  # Be + D2O
        self.sk1 = sk1
        self.sk2 = sk2
        self.sk3 = sk3

    def eta(self, ratio, p1, p2, p3):
    
        a = ratio * (self.k1 * p1 + self.k2 * p2)
        b = self.k3 * p3
        c = 1 / (self.k2 * p2 * (1 + ratio))
    
        output = (a - b) * c
    
        return output
    
    def eta_error(self, ratio, sratio, p1, p2, p3, sp1, sp2, sp3):
    
        aux0 = (self.k2 ** -2.) * ((p1 ** 2) * ((p2 ** -2.) * ((ratio ** 2) * (((1. + ratio) ** -2.) * (self.sk1 ** 2)))))
        aux1 = (((self.k2 ** -2.) * ((((self.k1 * p1) + (self.k2 * p2)) * ratio) - (self.k3 * p3))) / (1. + ratio)) / p2
        aux2 = (self.k1 ** 2) * ((self.k2 ** -2.) * ((p2 ** -2.) * ((ratio ** 2) * (((1. + ratio) ** -2.) * (sp1 ** 2)))))
        aux3 = (((p2 ** -2.) * ((((self.k1 * p1) + (self.k2 * p2)) * ratio) - (self.k3 * p3))) / (1. + ratio)) / self.k2
        aux4 = ((((1. + ratio) ** -2.) * ((((self.k1 * p1) + (self.k2 * p2)) * ratio) - (self.k3 * p3))) / p2) / self.k2
        aux5 = (((((((self.k1 * p1) + (self.k2 * p2)) / (1. + ratio)) / p2) / self.k2) - aux4) ** 2) * (sratio ** 2)
        aux6 = ((self.k2 ** -2.) * ((self.k3 ** 2) * ((p2 ** -2.) * (((1. + ratio) ** -2.) * (sp3 ** 2))))) + aux5
        aux7 = ((self.k2 ** -2.) * ((p2 ** -2.) * ((p3 ** 2) * (((1. + ratio) ** -2.) * (self.sk3 ** 2))))) + (
                    aux2 + ((((((ratio / (1. + ratio)) / p2) - aux3) ** 2) * (sp2 ** 2)) + aux6))
    
        output = _np.sqrt((aux0 + ((((((ratio / (1. + ratio)) / self.k2) - aux1) ** 2) * (self.sk2 ** 2)) + aux7)))
    
        return output
