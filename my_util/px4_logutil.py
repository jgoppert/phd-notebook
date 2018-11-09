import numpy as np
import csv
from collections import namedtuple
from math import factorial


def dict_to_namedtuple(name, d):
    """
    Convert a dict to a namedtuple, easier for plotting etc. in IPython.
    """
    return namedtuple(name, d.keys())(**d)


def nested_dict_to_nested_namedtuple(name, d):
    """
    Convert a nexted dict to a nested namedtuple, only handles two levels.
    """
    d2 = {}
    for key in d.keys():
        d2[key] = dict_to_namedtuple(key, d[key])
    return dict_to_namedtuple(name, d2)


def px4_log_to_dict(f):
    """
    Convert a px4 csv log file to a nested dictionary.
    """
    reader = csv.DictReader(f)
    d = {}
    for fieldname in reader.fieldnames:
        names = fieldname.split('_')
        msg = names[0]
        field = names[1]
        if msg not in d.keys():
            d[msg] = {}
        d[msg][field] = []
    for row in reader:
        for fieldname in reader.fieldnames:
            names = fieldname.split('_')
            msg = names[0]
            field = names[1]
            val = row[fieldname]
            if val == '':
                val = 0
            try:
                val = float(val)
            except ValueError:
                val = val
            d[msg][field].append(val)
    for fieldname in reader.fieldnames:
        names = fieldname.split('_')
        msg = names[0]
        field = names[1]
        try:
            d[msg][field] = np.array(d[msg][field], dtype='float')
        except:
            pass
    return d


def px4_log_to_namedtuple(f):
    """
    Convert a px4 csv log file to a nested named tuple, easier for plotting.
    """
    return nested_dict_to_nested_namedtuple('px4_log', px4_log_to_dict(f))

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
