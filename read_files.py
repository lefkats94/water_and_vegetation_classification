from osgeo import gdal
import numpy as np

def read_sentinel2_bands(inputs_directory):
    path = inputs_directory
    
    # 10m resolution bands
    blue=gdal.Open(path + "B02.tif")
    srcband = blue.GetRasterBand(1)
    blue_array = srcband.ReadAsArray(0, 0, blue.RasterXSize, blue.RasterYSize).astype(np.float)

    green=gdal.Open(path + "B03.tif")
    srcband = green.GetRasterBand(1)
    green_array = srcband.ReadAsArray(0, 0, green.RasterXSize, green.RasterYSize).astype(np.float)

    red=gdal.Open(path + "B04.tif")
    srcband = red.GetRasterBand(1)
    red_array = srcband.ReadAsArray(0, 0, red.RasterXSize, red.RasterYSize).astype(np.float)

    nir=gdal.Open(path + "B08.tif")
    srcband = nir.GetRasterBand(1)
    nir_array = srcband.ReadAsArray(0, 0, nir.RasterXSize, nir.RasterYSize).astype(np.float)

    # 20m resolution bands
    red_edge1 = gdal.Open(path + "B05.tif")
    srcband = red_edge1.GetRasterBand(1)
    red_edge1_array = srcband.ReadAsArray(0, 0, red_edge1.RasterXSize, red_edge1.RasterYSize).astype(np.float)

    red_edge3 = gdal.Open(path + "B07.tif")
    srcband = red_edge3.GetRasterBand(1)
    red_edge3_array = srcband.ReadAsArray(0, 0, red_edge3.RasterXSize, red_edge3.RasterYSize).astype(np.float)

    swir = gdal.Open(path + "B11.tif")
    srcband = swir.GetRasterBand(1)
    swir_array = srcband.ReadAsArray(0, 0, swir.RasterXSize, swir.RasterYSize).astype(np.float)

    swir2 = gdal.Open(path + "B12.tif")
    srcband = swir2.GetRasterBand(1)
    swir2_array = srcband.ReadAsArray(0, 0, swir2.RasterXSize, swir2.RasterYSize).astype(np.float)

    narrow_nir = gdal.Open(path + "B8A.tif")
    srcband = narrow_nir.GetRasterBand(1)
    narrow_nir_array = srcband.ReadAsArray(0, 0, narrow_nir.RasterXSize, narrow_nir.RasterYSize).astype(np.float)

    # resize 20m bands to 10m
    red_edge1_array = srcband.ReadAsArray(0, 0, red_edge1.RasterXSize, red_edge1.RasterYSize, blue.RasterXSize, blue.RasterYSize).astype(np.float)
    red_edge3_array = srcband.ReadAsArray(0, 0, red_edge3.RasterXSize, red_edge3.RasterYSize, blue.RasterXSize, blue.RasterYSize).astype(np.float)
    swir_array = srcband.ReadAsArray(0, 0, swir.RasterXSize, swir.RasterYSize, blue.RasterXSize, blue.RasterYSize).astype(np.float)
    swir2_array = srcband.ReadAsArray(0, 0, swir2.RasterXSize, swir2.RasterYSize, blue.RasterXSize, blue.RasterYSize).astype(np.float)
    narrow_nir_array = srcband.ReadAsArray(0, 0, narrow_nir.RasterXSize, narrow_nir.RasterYSize, blue.RasterXSize, blue.RasterYSize).astype(np.float)

    return blue_array, red_array, green_array, nir_array, red_edge1_array, red_edge3_array, swir_array, swir2_array, narrow_nir_array