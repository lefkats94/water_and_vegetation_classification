import numpy as np
from scipy.signal import lfilter

def clean_outliers(band_array):
    band_array[band_array < -1000000] = 0
    band_array[band_array > 1000000] = 0
    band_array[np.isnan(band_array)] = 0
    band_array[np.isinf(band_array)] = 0
    return band_array

def calculate_rescaled_array(band_array):
    rescaled_array = (band_array-np.amin(band_array)) / np.amax(band_array)
    return rescaled_array

def moving_average_smooth(Y, width = 5):
    n = len(Y)
    c = lfilter(np.ones(width)/width,1,Y)
    cbegin = np.cumsum(Y[0:width-2])
    cbegin_div = np.arange(1,(width-1),2).astype(float)
    cbegin = cbegin[::2]/cbegin_div
    cend = np.cumsum(Y[(n-1):(n-width+3-2):-1])
    cend_div = np.arange(width-2, 0, -2).astype(float)
    cend = cend[(n)::-2]/cend_div
    c = np.concatenate((cbegin, c[(width-1)::],cend))
    return c

def calculate_ndvi(red_array, nir_array):
    ndvi_array = (nir_array - red_array) / (nir_array + red_array)
    return ndvi_array