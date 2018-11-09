import struct
import numpy as np
from xml.etree.ElementTree import ElementTree
import copy
import matplotlib.pyplot as plt
import scipy.signal
import os
import pickle


def load_calibration(file_name):
    tree = ElementTree(file=file_name)
    d = {}
    for axis in tree.findall('*.UserAxis'):
        d[axis.attrib['Name']] = {
            'values': np.array(axis.attrib['values'].split()).astype(float),
            'max': float(axis.attrib['max'])}
    cal_matrix = np.array([
        d['Fx']['values'], d['Fy']['values'], d['Fz']['values'],
        d['Tx']['values'], d['Ty']['values'], d['Tz']['values']]).T
    max_matrix = np.array([
        d['Fx']['max'], d['Fy']['max'], d['Fz']['max'],
        d['Tx']['max'], d['Ty']['max'], d['Tz']['max']])
    return cal_matrix, max_matrix


def aiscan(daq, fs, tf):
    """
    daq : the daq
    fs : sampling frequency
    tf : final time
    """
    ch_lo = 0
    ch_hi = 7
    v_range = 5
    n_ch = ch_hi - ch_lo + 1
    samples = int(tf*fs)
    setup_str = """
    AISCAN:XFRMODE=BLOCKIO
    AI:CHMODE=DIFF
    AI{{0}}:RANGE=BIP{v_range}V
    AI{{1}}:RANGE=BIP{v_range}V
    AI{{2}}:RANGE=BIP{v_range}V
    AI{{3}}:RANGE=BIP{v_range}V
    AI{{4}}:RANGE=BIP{v_range}V
    AI{{5}}:RANGE=BIP{v_range}V
    AI{{6}}:RANGE=BIP{v_range}V
    AI{{7}}:RANGE=BIP{v_range}V
    AISCAN:LOWCHAN={ch_lo:d}
    AISCAN:HIGHCHAN={ch_hi:d}
    AISCAN:SAMPLES={samples:d}
    AISCAN:RATE={fs:d}
    """.format(**locals())
    for msg in setup_str.split():
        try:
            daq.send_message(msg)
        except IOError:
            raise IOError('send failed for msg: {:s}'.format(msg))
    bytes_in_uint16 = 2
    n_bytes = n_ch*samples*bytes_in_uint16

    daq.send_message('AISCAN:RESET')
    daq.send_message('AISCAN:START')
    packet = daq._ep_in.read(n_bytes, timeout=10000)
    data = np.reshape(struct.unpack('H'*(len(packet)/2), packet),
                      (len(packet)/2/n_ch, n_ch))
    t = np.linspace(0, tf, tf*fs)
    return {
        't': t,
        'tf': tf,
        'ai': data,
        'fs': fs,
    }


def normalize_signal(sig):
    sig2 = sig - sig.mean()
    max = sig2.max()
    return sig2/max


def apply_calibration(data, offsets, cal_matrix):
    x = data['ai'].astype(float)
    x -= offsets
    x[:, :6] = x[:, :6].dot(cal_matrix)*5/2**15
    return {
        't': data['t'],
        'tf': data['tf'],
        'fs': data['fs'],
        'F': x[:, 0:3],
        'M': x[:, 3:6],
        'ch6': x[:, 6],
        'ch7': x[:, 7],
    }


def pwm_to_level(pwm):
    level = np.zeros(len(pwm))
    state = 0
    start = 0
    for i in range(len(pwm)):
        new_state = 1 if pwm[i] > 0 else 0
        # rising edge
        if state == 0 and new_state == 1:
            start = i
        # falling edge
        elif state == 1 and new_state == 0:
            level[start-1:i+2] = i-start
        state = new_state
    return level


def aiscan_or_load_data(filename, daq, fs, tf, offsets, cal_matrix, force_run=False):
    if force_run or (not os.path.isfile(filename)):
        data = apply_calibration(
            aiscan(daq, fs, tf), offsets, cal_matrix)
        with open(filename,  'w') as save_file:
            pickle.dump(data, save_file)
    else:
        print 'loading'
        with open(filename, 'r') as load_file:
            data = pickle.load(load_file)
    return data


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y


def low_pass_data(data, cut_freq):
    d = copy.deepcopy(data)
    for key in ['F', 'M']:
        for dim in range(3):
            d[key][:, dim] = butter_lowpass_filter(
                data[key][:, dim], cut_freq, data['fs'])
    return d


def plot_forces(data):
    N2gm = 1000/9.8
    t = data['t']

    for dim in range(3):
        plt.plot(t, data['F'][:, dim]*N2gm)
    plt.ylabel('gm')
    plt.xlabel('t, sec')
    plt.title('force/g history')
    plt.legend(r'$F_x$ $F_y$ $F_z$'.split(), ncol=3, loc='best')
    plt.grid()


def plot_moments(data):
    N2gm = 1000/9.8
    t = data['t']

    for dim in range(3):
        plt.plot(t, data['M'][:, dim]*N2gm/100)
    plt.ylabel('gm-cm')
    plt.xlabel('t, sec')
    plt.title('moment/g history')
    plt.legend(r'$M_x$ $M_y$ $M_z$'.split(), ncol=3, loc='best')
    plt.grid()


def print_mean_force_moment(data):
    print 'F: {:s} N'.format(data['F'].mean(0))
    print 'M: {:s} N-m'.format(data['M'].mean(0))
