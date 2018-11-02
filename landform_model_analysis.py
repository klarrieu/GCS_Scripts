from file_functions import *
from classify_landforms_GUI import *
from extract_channel_dims_GUI import *
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Tkinter import *
import logging
arcpy.CheckOutExtension('Spatial')

# run extract channel dims and classify landforms on all flt tuflow velocity outputs between base and flood flow



def landform_stratified_velocities(station_lines, detrended_DEM, vel_ras_list, buffer_size='', rm_up_length=0, rm_down_length=0):
    '''
    Computes and plots landform stratified velocities

    Args:
        vel_ras_list: list of velocity raster files
        detrended DEM: name of raster file for detrended DEM
        station lines (str): name of shapefile for station XS lines
        buffer size (float): if ETGeoWizards was used to make station lines, it will automatically get buffer size.
        rm_up_length (float): distance upstream to remove from analysis (if there are edge effects)
        rm_down_length (float): distance downstream to remove from analysis (if there are edge effects)

    Returns:
        A plot of landform stratified average velocity vs discharge
    '''

    tables = extract_channel_data(station_lines=station_lines, detrended_DEM=detrended_DEM, wetted_polygons_list=vel_ras_list,
                         buffer_size=buffer_size, rm_up_length=rm_up_length, rm_down_length=rm_down_length)

    # get output tables from extract_channel_data
    # tables = [vel_ras.replace('.flt', '_XS_joined_table.csv') for vel_ras in vel_ras_list]

    main_classify_landforms(tables=tables, w_field='W', z_field='Z', dist_field='dist_down', make_plots=False)

    mu_XS_list = [vel_ras.replace('.flt', '_XS.shp') for vel_ras in vel_ras_list]

    for mu_poly, vel_ras in zip(mu_XS_list, vel_ras_list):
        arcpy.sa.ZonalStatisticsAsTable(mu_poly, 'code', vel_ras, os.path.join(os.path.dirname(vel_ras), 'v_avg_%s' % os.path.basename(vel_ras).replace('.flt','')), statistics_type='MEAN')


    # read in zonal stats tables

    # make plot with line for each code, points are discharge and corresponding velocity
    '''
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot()
    plt.show()
    '''

if __name__ == '__main__':

    # initialize logger
    init_logger(__file__)

    #get velocity raster from final timestep for each tuflow run

    vel_ras_list = get_all_files('G:\\Ken_RB\\c01\\Tuflow\\results\\003\\', suffix='V_002_00.flt')+get_all_files('G:\\Ken_RB\\c01\\Tuflow\\results\\003\\', suffix='V_002_00.hdr')
    # copy rasters to folder for analysis
    for vel_ras in vel_ras_list:
        shutil.copyfile(vel_ras, 'G:\\Ken_RB\\c01\\mu_vel\\' + os.path.basename(vel_ras).replace('.', 'pt', 1).replace('RB_', '').replace('_003_V_002_00', ''))
    vel_ras_list = get_all_files('G:\\Ken_RB\\c01\\mu_vel\\', suffix='.flt')


    landform_stratified_velocities(station_lines='G:\\Ken_RB\\c01\\stationing\\stations_100.shp',
                                   detrended_DEM='G:\\Ken_RB\\c01\\DEM\\Detrending\\556\\detrended_DEM.tif',
                                   vel_ras_list=vel_ras_list,
                                   rm_up_length=55,
                                   rm_down_length=55)

    '''
    root = Tk()
    root.wm_title('Landform Stratified Velocities')
    '''
