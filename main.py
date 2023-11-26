from read_files import read_sentinel2_bands
from help_methods import clean_outliers, calculate_ndvi, calculate_rescaled_array, moving_average_smooth
import conv_peakdet
import numpy as np

if __name__ == "__main__":
    # call the read_ndvi_data function
    bands_directory = input("Enter the NDVIS directory:")

    # correct the given directory
    bands_directory = bands_directory.replace("\\", "/")
    if not bands_directory.endswith("/"):
        bands_directory += "/"
    
    # read the necessary Sentinel-2 bands
    blue_array, red_array, green_array, nir_array, red_edge1_array, red_edge3_array, swir_array, swir2_array, narrow_nir_array = read_sentinel2_bands("C:/Users/lefkats_local/Desktop/test_watermask/")

    # clean the outliers of the swir bands
    swir_array = clean_outliers(swir_array)
    swir2_array = clean_outliers(swir2_array)

    #rescale the swir band
    swir_array = calculate_rescaled_array(swir_array)

    # flat the swir band
    swir_flat = swir_array.flatten()

    # calculate the swir histogram
    swir_hist, bin_edges = np.histogram(swir_flat[~np.isnan(swir_flat)], bins = 255)

    # identify the peaks and valleys in the swir's histogram
    swir_peaks, swir_valleys = conv_peakdet.peakdet(swir_hist, np.mean(swir_hist)/3)

    # smooth the histogram
    swir_hist_smoothed = moving_average_smooth(swir_hist)

    # identify the peaks and valleys in the swir's smoothed histogram
    swir_smoothed_peaks, swir_smoothed_valleys = conv_peakdet.peakdet(swir_hist_smoothed, np.mean(swir_hist_smoothed)/3)

    # calculate the ndvi
    ndvi_array = calculate_ndvi(red_array, nir_array)
    print(swir_smoothed_valleys)