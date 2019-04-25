import numpy as _np

# Theory Values
k1 = 2.03e-9
k2 = 2.18e-9
k3 = 2.29e-9
sk1 = 0.04e-9
sk2 = 0.07e-9
sk3 = 0.05e-9


def eta_no_initial(ratio, k1, k2, k3, p1, p2, p3):
    """

    :param ratio:
    :param k1:
    :param k2:
    :param k3:
    :param p1:
    :param p2:
    :param p3:
    :return:
    """
    top = k1 * p1 + k2 * p2 - k3 * p3 * ratio
    bot = k2 * p2 * (1 + ratio)

    output = top / bot

    return output


def eta_error_no_initial(ratio, k1, k2, k3, p1, p2, p3, sratio, sp1, sp2, sp3):
    """

    :param ratio:
    :param k1:
    :param k2:
    :param k3:
    :param p1:
    :param p2:
    :param p3:
    :param sratio:
    :param sp1:
    :param sp2:
    :param sp3:
    :return:
    """
    aux0 = ((((-((1. +ratio) ** -2.) * (((k1 * p1) + (k2 * p2)) - (k3 * (p3 *ratio))))) / p2) / k2) - ((((k3 * p3) / (1. +ratio)) / p2) / k2)
    aux1 = (((1. +ratio) ** -1.) / p2) - ((((p2 ** -2.) * (((k1 * p1) + (k2 * p2)) - (k3 * (p3 * ratio)))) / (1. + ratio)) / k2)
    aux2 = ((k2 ** -2.) * ((k3 ** 2) * ((p2 ** -2.) * ((sp3 ** 2) * ((ratio ** 2) * ((1. + ratio) ** -2.)))))) + (((sratio ** 2) * (aux0 ** 2)) + ((sp2 ** 2) * (aux1 ** 2)))
    
    output = _np.sqrt((((k1 ** 2) * ((k2 ** -2.) * ((p2 ** -2.) * ((sp1 ** 2) * ((1. +ratio) ** - 2.))))) + aux2))
    
    return output


def eta(t, ratio, k1, k2, k3, p1, p2, p3, Be0, BeOH0, BeOD0):
    """

    :param t:
    :param ratio:
    :param k1:
    :param k2:
    :param k3:
    :param p1:
    :param p2:
    :param p3:
    :param Be0:
    :param BeOH0:
    :param BeOD0:
    :return:
    """
    aux0 = (BeOH0 * (k1 * (p1 *ratio))) + ((BeOH0 * (k2 * (p2 *ratio))) + ((Be0 * (k3 * (p3 *ratio))) + (BeOH0 * (k3 * (p3 *ratio)))))
    aux1 = (Be0 * ((_np.exp(((((-k3 * p3) - (k2 * p2)) - (k1 * p1)) * t))) * (k1 * p1))) + ((Be0 * ((_np.exp(((((-k3 * p3) - (k2 * p2)) - (k1 * p1)) * t))) * (k2 * p2))) + aux0)
    aux2 = (aux1 - (Be0 * ((_np.exp(((((-k3 * p3) - (k2 * p2)) - (k1 * p1)) * t))) * (k3 * (p3 *ratio))))) - (BeOD0 * (k3 * p3))
    aux3 = ((((aux2 - (BeOD0 * (k2 * p2))) - (Be0 * (k2 * p2))) - (BeOD0 * (k1 * p1))) - (Be0 * (k1 * p1))) / (1. +ratio)
    
    output = (((aux3 / p2) / k2) / (-1. + (_np.exp(((((-k3 * p3) - (k2 * p2)) - (k1 * p1)) * t))))) / Be0
    
    return output


def eta_error(t, ratio, k1, k2, k3, p1, p2, p3, Be0, BeOH0, BeOD0, sratio, sp1, sp2, sp3, sBe0, sBeOH0, sBeOD0):
    """
    
    :param ratio: 
    :param sratio: 
    :param p1: 
    :param p2: 
    :param p3: 
    :param Be0: 
    :param BeOH0: 
    :param BeOD0: 
    :param sp1: 
    :param sp2: 
    :param sp3: 
    :param sBe0: 
    :param sBeOH0: 
    :param sBeOD0: 
    :return: 
    """
    
    aux0=(p2**-2.)*((((((-k3*p3)-(k2*p2))-(k1*p1))**2))*((sBeOD0**2)*((1.+ratio)**-2.)))
    aux1=((-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))**-2.)*((k2**-2.)*aux0)
    aux2=(sBeOH0**2)*(((1.+ratio)**-2.)*((((k1*(p1*ratio))+((k2*(p2*ratio))+(k3*(p3*ratio))))**2)))
    aux3=((-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))**-2.)*((k2**-2.)*((p2**-2.)*aux2))
    aux4=((BeOH0*(k1*p1))+((BeOH0*(k2*p2))+((Be0*(k3*p3))+(BeOH0*(k3*p3)))))-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*p3)))
    aux5=(((aux4/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux6=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux7=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux6)
    aux8=(aux7-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux9=((1.+ratio)**-2.)*((((aux8-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1)))
    aux10=(aux5/Be0)-((((aux9/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))))/Be0)
    aux11=((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1))+(((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2))+(k3*(p3*ratio)))
    aux12=(aux11-((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio))))-(k2*p2)
    aux13=((((aux12-(k1*p1))/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux14=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux15=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux14)
    aux16=(aux15-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux17=(((aux16-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1))
    aux18=(((((Be0**-2.)*aux17)/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux19=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux20=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux19)
    aux21=(aux20-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux22=(((aux21-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1))
    aux23=(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(((-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))**-2.)*(k1*(t*aux22)))
    aux24=(BeOH0*(k1*ratio))+(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*(k3*(p3*(t*ratio))))))
    aux25=((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*k1))+aux24)-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*(k2*(p2*t)))))
    aux26=aux25-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*((k1**2)*(p1*t))))
    aux27=(((((aux26-(BeOD0*k1))-(Be0*k1))/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux28=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux29=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux28)
    aux30=(aux29-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux31=(((aux30-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1))
    aux32=(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(((-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))**-2.)*(t*aux31))
    aux33=(BeOH0*(k2*ratio))+(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*(k3*(p3*(t*ratio))))))
    aux34=((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*k2))+aux33)-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*((k2**2)*(p2*t))))
    aux35=aux34-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*(k2*(p1*t)))))
    aux36=(((((aux35-(BeOD0*k2))-(Be0*k2))/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux37=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux38=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux37)
    aux39=(aux38-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux40=(((aux39-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1))
    aux41=((((p2**-2.)*aux40)/(1.+ratio))/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux42=(BeOH0*(k1*(p1*ratio)))+((BeOH0*(k2*(p2*ratio)))+((Be0*(k3*(p3*ratio)))+(BeOH0*(k3*(p3*ratio)))))
    aux43=(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*p1)))+((Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*p2)))+aux42)
    aux44=(aux43-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*(p3*ratio)))))-(BeOD0*(k3*p3))
    aux45=(((aux44-(BeOD0*(k2*p2)))-(Be0*(k2*p2)))-(BeOD0*(k1*p1)))-(Be0*(k1*p1))
    aux46=(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(((-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))**-2.)*(k3*(t*aux45)))
    aux47=(BeOH0*(k3*ratio))+(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*((k3**2)*(p3*(t*ratio)))))
    aux48=((Be0*(k3*ratio))+aux47)-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k3*ratio)))
    aux49=aux48-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k2*(k3*(p2*t)))))
    aux50=aux49-(Be0*((_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t)))*(k1*(k3*(p1*t)))))
    aux51=((((aux50-(BeOD0*k3))/(1.+ratio))/p2)/k2)/(-1.+(_np.exp(((((-k3*p3)-(k2*p2))-(k1*p1))*t))))
    aux52=((sp2**2)*(((((((aux32/(1.+ratio))/p2)/Be0)+(aux36/Be0))-(aux41/Be0))**2)))+((sp3**2)*(((((((aux46/(1.+ratio))/p2)/k2)/Be0)+(aux51/Be0))**2)))
    aux53=((sBe0**2)*((((aux13/Be0)-aux18)**2)))+(((sp1**2)*(((((((aux23/(1.+ratio))/p2)/k2)/Be0)+(aux27/Be0))**2)))+aux52)
    
    output=_np.sqrt((((Be0**-2.)*aux1)+(((Be0**-2.)*aux3)+(((sratio**2)*(aux10**2))+aux53))))
    
    return output
