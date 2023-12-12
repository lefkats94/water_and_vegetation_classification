import numpy as np
from scipy.signal import lfilter
from osgeo import gdal

def get_metadata(inputs_directory):
    """
    Retrieves metadata from a Sentinel-2 band.

    Parameters:
    - inputs_directory (str): The directory containing the Sentinel-2 bands.

    Returns:
    GDAL Dataset: Metadata information for the specified band.
    """
    random_band = gdal.Open(inputs_directory + "B04.tif")
    return random_band

def read_sentinel2_bands(inputs_directory):
    """
    Reads and returns Sentinel-2 bands.

    Parameters:
    - inputs_directory (str): The directory containing the Sentinel-2 bands.

    Returns:
    Tuple of numpy arrays: (red_array, green_array, nir_array, red_edge1_array, swir_array, swir2_array)
    """
    # 10m resolution bands
    green = gdal.Open(inputs_directory + "B03.tif")
    srcband = green.GetRasterBand(1)
    green_array = srcband.ReadAsArray(0, 0, green.RasterXSize, green.RasterYSize).astype(np.float64)

    red = gdal.Open(inputs_directory + "B04.tif")
    srcband = red.GetRasterBand(1)
    red_array = srcband.ReadAsArray(0, 0, red.RasterXSize, red.RasterYSize).astype(np.float64)

    nir = gdal.Open(inputs_directory + "B08.tif")
    srcband = nir.GetRasterBand(1)
    nir_array = srcband.ReadAsArray(0, 0, nir.RasterXSize, nir.RasterYSize).astype(np.float64)

    # 20m resolution bands
    red_edge1 = gdal.Open(inputs_directory + "B05.tif")
    srcband = red_edge1.GetRasterBand(1)
    red_edge1_array = srcband.ReadAsArray(0, 0, red_edge1.RasterXSize, red_edge1.RasterYSize, nir.RasterXSize, nir.RasterYSize).astype(np.float64)

    swir = gdal.Open(inputs_directory + "B11.tif")
    srcband = swir.GetRasterBand(1)
    swir_array = srcband.ReadAsArray(0, 0, swir.RasterXSize, swir.RasterYSize, nir.RasterXSize, nir.RasterYSize).astype(np.float64)

    swir2 = gdal.Open(inputs_directory + "B12.tif")
    srcband = swir2.GetRasterBand(1)
    swir2_array = srcband.ReadAsArray(0, 0, swir2.RasterXSize, swir2.RasterYSize, nir.RasterXSize, nir.RasterYSize).astype(np.float64)

    return red_array, green_array, nir_array, red_edge1_array, swir_array, swir2_array

def clean_outliers(band_array):
    """
    Cleans outliers in the input band array.

    Parameters:
    - band_array (numpy.ndarray): Input band array.

    Returns:
    numpy.ndarray: Cleaned band array.
    """
    band_array[band_array < -1000000] = 0
    band_array[band_array > 1000000] = 0
    band_array[np.isnan(band_array)] = 0
    band_array[np.isinf(band_array)] = 0
    return band_array

def calculate_rescaled_array(band_array):
    """
    Calculates and returns a rescaled version of the input band array.

    Parameters:
    - band_array (numpy.ndarray): Input band array.

    Returns:
    numpy.ndarray: Rescaled band array.
    """
    rescaled_array = (band_array - np.amin(band_array)) / np.amax(band_array)
    return rescaled_array

def moving_average_smooth(Y, width=5):
    """
    Applies a moving average smoothing to the input array.

    Parameters:
    - Y (numpy.ndarray): Input array.
    - width (int): Width of the moving average window.

    Returns:
    numpy.ndarray: Smoothed array.
    """
    n = len(Y)
    c = lfilter(np.ones(width)/width, 1, Y)
    cbegin = np.cumsum(Y[0:width-2])
    cbegin_div = np.arange(1, (width-1), 2).astype(float)
    cbegin = cbegin[::2] / cbegin_div
    cend = np.cumsum(Y[(n-1):(n-width+3-2):-1])
    cend_div = np.arange(width-2, 0, -2).astype(float)
    cend = cend[(n)::-2] / cend_div
    c = np.concatenate((cbegin, c[(width-1)::], cend))
    return c


def peakdet(v=None, delta=None, x=None):
    
    '''
    PEAKDET Detect peaks in a vector
            [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
            maxima and minima ("peaks") in the vector V.
            MAXTAB and MINTAB consists of two columns. Column 1
            contains indices in V, and column 2 the found values.
          
            With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
            in MAXTAB and MINTAB are replaced with the corresponding
            X-values.
    
            A point is considered a maximum peak if it has the maximal
            value, and was preceded (to the left) by a value lower by
            DELTA.
    '''
    maxtab = []
    mintab = []

    v = np.array(v)
       
    #if args < 3 actually
    if x is None:
        x = range(0, len(v))
        x = np.array(x)
    else:
        if len(v) != len(x):
            print('Input vectors v and x must have same length')
            quit()
    delta = np.array(delta)
    
    if (delta.size)>1:
        print('Input argument DELTA must be a scalar')
        quit()
    if delta <= 0:
        print("Input argument DELTA must be positive")
        quit()
    
    mn = float('Inf')
    mx = -float('Inf')

    mnpos = float('NaN')
    mxpos = float('NaN')

    lookformax = 1

    for i in range (1,len(v)+1):
      temp = v[i-1]
      if temp > mx:
        mx = temp
        mxpos = x[i-1]
      if temp < mn:
        mn = temp
        mnpos = x[i-1]

      
      if lookformax:
        if temp < mx - delta:
            temp_maxtab = []
            temp_maxtab.append(mxpos)
            temp_maxtab.append(mx)
            maxtab.append(temp_maxtab)
            mn = temp
            mnpos = x[i-1]
            lookformax = 0
          
      else:
          if temp > mn + delta:
            temp_mintab = []
            temp_mintab.append(mnpos)
            temp_mintab.append(mn)
            mintab.append(temp_mintab)  
            mx = temp
            mxpos = x[i-1]
            lookformax = 1

    return maxtab, mintab
