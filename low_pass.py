
import numpy as np
from scipy.signal import butter, filtfilt

def butter_lowpass_filter(data, cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y


def avrx(q, geom):
    # (jm, im, lm) = q.shape
    (jm, im) = q.shape

    nmax = im / 2
    bysn = 1 / np.sin(np.pi / im * np.arange(1, nmax+1))

    drat = geom.dy / geom.dx_j

    sm = 1 - bysn / drat
    smmz = 1 - np.maximum(sm, np.zeros_like(sm))
    smmz = np.insert(smmz, 0, 1, 1)


    f_q = np.fft.rfft(q.m)

    ratios = np.fft.rfftfreq(im, geom.dx_j) * geom.dy
    ratio_mult = np.zeros_like(ratios)
    ratio_mult[ratios <= 0.5] = 1

    f_q_f = f_q * ratio_mult

    q_f = np.fft.irfft(f_q_f) * q.u

    return q_f


def arakawa_1977(q, geom):
    """
    implements the smoothing from Arakawa and Lamb 1977, "THE UCLA GENERAL CIRCULATION MODEL"
    page 248, Longitudinal Averaging of Selected Terms Near the Poles

    The method devised to allow the use of a longer dt in the model is to
    smooth the longitudinal pressure gradient in the momentum equation and
    the longitudinal divergence in the continuity equation with a longitudinal
    averaging operator.
    """
    # (jm, im, lm) = q.shape
    (lm, jm, im) = q.shape
    drat = geom.dy / geom.dx_j


    # compute formula 325
    # a * dw is the latitudinal grid size, aka d*

    nmax = im / 2
    bysn = 1 / np.sin(np.pi / im * np.arange(1, nmax+1))

    sm = 1 - bysn / drat
    smmz = 1 - np.maximum(sm, np.zeros_like(sm))
    smmz = np.insert(smmz, 0, 1, -1)


    f_q = np.fft.rfft(q.m)
    f_q_f = f_q * smmz.m
    q_f = np.fft.irfft(f_q_f) * q.u
    return q_f
