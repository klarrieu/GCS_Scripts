from file_functions import *
'''
inputs:
-
-model depth, velocity rasters (bankfull)
-thalweg shapefile
-XS shapefile

aggregate metrics for each reach (at bankfull):
-slope
-sinuosity
-average depth
-average width

if GCS landform attributes in XS shapefile:
-characteristic variables by landform?
'''

def agg_metrics(d_ras, v_ras, thalweg_shp):
    '''
    calculates slope, sinuosity, average depth, average width

    Args:
        d_ras: depth raster
        v_ras: velocity raster
        thalweg_shp: thalweg shapefile
    '''

    # slice XS attribute table by reach
    for

    # calculate mean depth


    # calculate mean width: (average XS area)/(XS length)