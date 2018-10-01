import numpy as _np

k1 = 2.03e-9
k2 = 2.18e-9
k3 = 2.29e-9
sk1 = 0.04e-9
sk2 = 0.07e-9
sk3 = 0.05e-9


def eta(t, ratio, k1, k2, k3, p1, p2, p3, Be_0, BeOH_0, BeOD_0):
    """

    :param t:
    :param ratio:
    :param k1:
    :param k2:
    :param k3:
    :param p1:
    :param p2:
    :param p3:
    :param Be_0:
    :param BeOH_0:
    :param BeOD_0:
    :return:
    """
    aux0 = ((((Be_0 + BeOH_0) * ((k1 * p1) + (k2 * p2))) + (BeOH_0 * (k3 * p3))) * ratio)-(BeOD_0 * ((k1 * p1) + ((k2
                                                                                                                   * p2) + (k3 * p3))))
    aux1 = (Be_0 * (k3 * p3)) + ((_np.exp((((k1 * p1) + ((k2 * p2) + (k3 * p3))) * t))) * (aux0 - (Be_0 * (k3 * p3))))
    aux2 = ((((aux1 - (Be_0 * (((k1 * p1) + (k2 * p2)) * ratio))) / (1. + ratio)) / p2) / k2) / (-1. + (_np.exp((((k1 *
                                                                                                                 p1) + ((k2 * p2) + (k3 * p3))) * t))))
    output = aux2 / Be_0

    return output


# def eta_error(ratio, sratio, p1, p2, p3, sp1, sp2, sp3):
#     """
#
#     :param ratio: ratio of BeOD/BeOH TOF signals
#     :param sratio: error in ratio of BeOD/BeOH TOF signals
#     :param p1: H2O pressure
#     :param p2: HOD pressure
#     :param p3: D2O pressure
#     :param sp1: error in H2O pressure
#     :param sp2: error in HOD pressure
#     :param sp3: error in D2O pressure
#     :return: error in BeOD branching ratio
#     """
#     aux0 = (k2 ** -2.) * ((p1 ** 2) * ((p2 ** -2.) * ((ratio ** 2) * (((1. + ratio) ** -2.) * (self.sk1 ** 2)))))
#     aux1 = (((k2 ** -2.) * ((((k1 * p1) + (k2 * p2)) * ratio) - (k3 * p3))) / (1. + ratio)) / p2
#     aux2 = (k1 ** 2) * ((k2 ** -2.) * ((p2 ** -2.) * ((ratio ** 2) * (((1. + ratio) ** -2.) * (sp1 ** 2)))))
#     aux3 = (((p2 ** -2.) * ((((k1 * p1) + (k2 * p2)) * ratio) - (k3 * p3))) / (1. + ratio)) / k2
#     aux4 = ((((1. + ratio) ** -2.) * ((((k1 * p1) + (k2 * p2)) * ratio) - (k3 * p3))) / p2) / k2
#     aux5 = (((((((k1 * p1) + (k2 * p2)) / (1. + ratio)) / p2) / k2) - aux4) ** 2) * (sratio ** 2)
#     aux6 = ((k2 ** -2.) * ((k3 ** 2) * ((p2 ** -2.) * (((1. + ratio) ** -2.) * (sp3 ** 2))))) + aux5
#     aux7 = ((k2 ** -2.) * ((p2 ** -2.) * ((p3 ** 2) * (((1. + ratio) ** -2.) * (self.sk3 ** 2))))) + (aux2 + ((((((ratio / (1. + ratio)) / p2) - aux3) ** 2) * (sp2 ** 2)) + aux6))
#
#     output = _np.sqrt((aux0 + ((((((ratio / (1. + ratio)) / k2) - aux1) ** 2) * (self.sk2 ** 2)) + aux7)))
#
#     return output
