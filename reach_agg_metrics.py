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

def agg_metrics(d_ras, v_ras, thalweg_shp, reach_breaks):
    '''calculates slope, sinuosity, average depth, average width'''

    # slice XS attribute table into

    # calculate mean depth


    # calculate mean width: (average XS area)/(XS length)