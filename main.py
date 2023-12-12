from methods import (
    get_metadata,
    clean_outliers,
    peakdet,
    moving_average_smooth,
    read_sentinel2_bands
)
import numpy as np
from osgeo import gdal
from PIL import Image

class WaterMaskClassifier:
    def __init__(self, bands_directory):
        self.bands_directory = bands_directory.replace("\\", "/")
        if not self.bands_directory.endswith("/"):
            self.bands_directory += "/"
        self.red_array, self.green_array, self.nir_array, self.red_edge1_array, self.swir_array, self.swir2_array = self.preprocess()

    def preprocess(self):
        """
        Read Sentinel-2 bands and preprocess SWIR arrays by cleaning outliers.

        Returns:
        Tuple of numpy arrays: (red_array, green_array, nir_array, red_edge1_array, swir_array, swir2_array)
        """
        red_array, green_array, nir_array, red_edge1_array, swir_array, swir2_array = read_sentinel2_bands(self.bands_directory)
        swir_array = clean_outliers(swir_array)
        swir2_array = clean_outliers(swir2_array)
        return red_array, green_array, nir_array, red_edge1_array, swir_array, swir2_array

    def swir_thresholds(self):
        """
        Detects and returns the positions of the first and second deep valleys in a Short-Wave Infrared (SWIR) array.

        Returns:
        Tuple of two values: (first_deep_valley, second_deep_valley)
        """
        swir_flat = self.swir_array.flatten()
        swir_hist, bin_edges = np.histogram(swir_flat[~np.isnan(swir_flat)], bins=255)
        swir_hist_smoothed = moving_average_smooth(swir_hist)
        swir_smoothed_peaks, swir_smoothed_valleys = peakdet(swir_hist_smoothed, np.mean(swir_hist_smoothed)/3)
        swir_valleys_smoothed_only = [x[0] for x in swir_smoothed_valleys]
        first_swir_threshold = np.amax(self.swir_array)*(swir_valleys_smoothed_only[0]/255) + np.amin(self.swir_array)
        second_swir_threshold = np.amax(self.swir_array)*(swir_valleys_smoothed_only[1]/255) + np.amin(self.swir_array)
        return first_swir_threshold, second_swir_threshold

    def ndvi_threshold(self):
        """
        Computes and returns the Normalized Difference Vegetation Index (NDVI) threshold from NIR and Red bands.

        Returns:
        Tuple of two values: (ndvi_threshold, ndvi)
        """
        ndvi = np.divide(self.nir_array - self.red_array, self.nir_array + self.red_array, out=np.ones_like(self.nir_array) * -2, where=(self.nir_array + self.red_array) != 0)
        ndvi_flattened = ndvi.flatten()
        ndvi_hist, _ = np.histogram(ndvi_flattened[~np.isnan(ndvi_flattened)], bins=201, range=[-1.0, 1.0])
        ndvi_peaks, ndvi_valleys = peakdet(ndvi_hist, np.mean(ndvi_hist)/10)
        ndvi_valleys_only = [x[0] for x in ndvi_valleys]
        ndvi_valleys_only_ranged = [(x - 100)/100.0 for x in ndvi_valleys_only]

        if len(ndvi_valleys_only_ranged) > 0:
            for ndvi_threshold in (x for x in ndvi_valleys_only_ranged if x > 0.3):
                break
        else:
            ndvi_threshold = -1

        return ndvi_threshold, ndvi

    def floating_vegetation(self):
        """
        Computes and returns a binary mask indicating floating vegetation areas.

        Returns:
        numpy.ndarray: Binary mask indicating floating vegetation areas.
        """
        floating_index = np.divide(self.red_edge1_array, self.swir_array, out=np.ones_like(self.swir_array) * -2, where=(self.red_edge1_array + self.swir_array) != 0)
        ndwi = np.divide(self.nir_array - self.swir_array, self.nir_array + self.swir_array, out=np.ones_like(self.green_array) * -2, where=(self.green_array + self.nir_array) != 0)
        floating_mask = np.full(floating_index.shape, False)

        condition_mask = (0.6 < floating_index) & (floating_index < 1.5) & (0.2 < ndwi) & (ndwi < 0.45) & (100 < self.swir2_array) & (self.swir2_array < 900)
        floating_mask[condition_mask] = True
        return floating_mask

    def watermask_calculation(self, ndvi_thres, floating_mask):
        """
        Computes and returns a classified water mask.

        Parameters:
        - ndvi_thres (float): NDVI threshold.
        - floating_mask (numpy.ndarray): Binary mask indicating floating vegetation areas.

        Returns:
        numpy.ndarray: Classified water mask.
        """
        first_swir_threshold, second_swir_threshold = self.swir_thresholds()
        classified_area = np.full(floating_mask.shape, 0)
        classified_area[self.swir_array < first_swir_threshold] = 1
        classified_area[(self.swir_array >= first_swir_threshold) & (self.swir_array <= second_swir_threshold) & (self.ndvi > ndvi_thres)] = 2
        classified_area[floating_mask] = 3
        return classified_area

    def create_outputs(self):
        """
        Creates output files including a GeoTIFF and a PNG image for the classified water mask.
        """
        inputs_directory = self.bands_directory
        inputs_directory = inputs_directory.replace("\\", "/")
        if not inputs_directory.endswith("/"):
            inputs_directory += "/"

        x, y = self.classified_area.shape

        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(inputs_directory + "classification.tif", y, x, 1, gdal.GDT_Float32)
        dataset.GetRasterBand(1).WriteArray(self.classified_area)
        dataset.GetRasterBand(1).SetNoDataValue(-1)

        random_band = get_metadata(inputs_directory)

        geotrans = random_band.GetGeoTransform()
        proj = random_band.GetProjection()
        dataset.SetGeoTransform(geotrans)
        dataset.SetProjection(proj)

        dataset = None

        rgb = np.zeros((x, y, 3))

        rgb[:,:,0] = np.where((self.classified_area == 0), 255, rgb[:,:,0])
        rgb[:,:,1] = np.where((self.classified_area == 0), 128, rgb[:,:,1])
        rgb[:,:,2] = np.where((self.classified_area == 0), 0, rgb[:,:,2])

        rgb[:,:,0] = np.where((self.classified_area == 1), 0, rgb[:,:,0])
        rgb[:,:,1] = np.where((self.classified_area == 1), 0, rgb[:,:,1])
        rgb[:,:,2] = np.where((self.classified_area == 1), 255, rgb[:,:,2])

        rgb[:,:,0] = np.where((self.classified_area == 2), 0, rgb[:,:,0])
        rgb[:,:,1] = np.where((self.classified_area == 2), 102, rgb[:,:,1])
        rgb[:,:,2] = np.where((self.classified_area == 2), 0, rgb[:,:,2])

        rgb[:,:,0] = np.where((self.classified_area == 3), 0, rgb[:,:,0])
        rgb[:,:,1] = np.where((self.classified_area == 3), 255, rgb[:,:,1])
        rgb[:,:,2] = np.where((self.classified_area == 3), 0, rgb[:,:,2])

        im = Image.fromarray(rgb.astype("uint8"))
        im.save(inputs_directory + "classification.png")

if __name__ == "__main__":
    print("""
    Sentinel-2 Bands Required:
    --------------------------

    - B03 (Green Band, 10m resolution)
    - B04 (Red Band, 10m resolution)
    - B08 (NIR Band, 10m resolution)
    - B05 (Red Edge 1 Band, 20m resolution)
    - B11 (SWIR Band, 20m resolution)
    - B12 (SWIR 2 Band, 20m resolution)

    Note: Ensure these Sentinel-2 bands (in .tif format) are available in the specified input directory.
    """)
    inputs_dir = input("insert the inputs directory:")
    water_mask_classifier = WaterMaskClassifier(inputs_dir)
    ndvi_thres, water_mask_classifier.ndvi = water_mask_classifier.ndvi_threshold()
    water_mask_classifier.floating_mask = water_mask_classifier.floating_vegetation()
    water_mask_classifier.classified_area = water_mask_classifier.watermask_calculation(ndvi_thres, water_mask_classifier.floating_mask)
    water_mask_classifier.create_outputs()
