import labrad

cxn = labrad.connect()
dv = cxn.data_vault


def cd(year, month, day, run, kind='images'):
    """

    Goes to ion_images file
    :param year: year (####)
    :param month: month (int)
    :param day: day (int)
    :param run: run (int)
    :param kind: datavault section
    :return:
    """
    year = str(year)
    month = '%02d' % month
    day = '%02d' % day
    trunk = year + '_' + month + '_' + day
    dv.cd([''])
    if kind == 'images':
        dv.cd([year, month, trunk, 'Be_reaction', 'images_%s' % run])
