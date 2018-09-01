# import tkFileDialog as _tk
# #
# #
# # def select_files():
# #     """
# #     Opens a dialog window to select files
# #     :return: tuple of file paths
# #     """
# #     files = _tk.askopenfilenames
# #
# #     return files


class DataObj(object):
    """
    General Data Class
    .d for data
    .p for parameters
    """
    def __init__(self, data=None, params=None):
        """
        Init function
        data and params usually dictionaries
        :param data: data
        :param params: parameters
        """
        self.d = data
        self.p = params

    def __setitem__(self, key, item):
        self.d[key] = item

    def __getitem__(self, item):
        try:
            return self.d[item]
        except (TypeError, KeyError):
            pass
        return self.p[item]
