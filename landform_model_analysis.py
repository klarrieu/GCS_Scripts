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

    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}
    color_code = {-2: 'black', -1: 'blue', 0: 'grey', 1: 'orange', 2: 'red'}

    tables = extract_channel_data(station_lines=station_lines, detrended_DEM=detrended_DEM, wetted_polygons_list=vel_ras_list,
                                  buffer_size=buffer_size, rm_up_length=rm_up_length, rm_down_length=rm_down_length)

    # get output tables from extract_channel_data
    tables = [vel_ras.replace('.flt', '_XS_joined_table.csv') for vel_ras in vel_ras_list]

    main_classify_landforms(tables=tables, w_field='W', z_field='Z', dist_field='dist_down', make_plots=False)

    mu_XS_list = [vel_ras.replace('.flt', '_XS.shp') for vel_ras in vel_ras_list]


    logging.info('Getting mean velocities by landform and discharge...')
    series = []  # fill with [discharge, mean_vel, code] lists
    for mu_poly, vel_ras in zip(mu_XS_list, vel_ras_list):

        # get mean velocity series
        mean_vel = arcpy.sa.ZonalStatisticsAsTable(mu_poly, 'code', vel_ras, os.path.join(os.path.dirname(vel_ras), 'v_avg_%s.dbf' % os.path.basename(vel_ras).replace('.flt','')), statistics_type='MEAN')
        mean_vel = arcpy.TableToExcel_conversion(mean_vel, str(mean_vel).replace('.dbf', '.xls'))
        df = pd.read_excel(str(mean_vel))
        q = os.path.basename(vel_ras).replace('cms.flt','').replace('pt', '.')
        vals = [[float(q), float(v), int(code)]
                for v, code in zip(df['MEAN'].tolist(), df['code'].tolist()) if code != -9999]
        series.extend(vals)

        # get 95th percentile velocity series
        # use con to get set of velocity pixels

    logging.info('OK')

    series = sorted(series)

    # make plot with line for each code, points are discharge and corresponding velocity
    fig, ax = plt.subplots(1, 1, sharex=True)
    fig.suptitle('Landform Stratified Velocity')

    for code in [-2, -1, 0, 1, 2]:
        qs, vs = map(list, zip(*[x[:2] for x in series if x[2] == code]))
        ax.semilogx(qs, vs, label=code_dict[code], color=color_code[code], marker='o', markersize=5)

    ax.set(ylabel=r'Mean Velocity $(m/s)$')
    ax.grid()
    '''
    ax2.set(xlabel=r'Discharge $(cms)$', ylabel=r'95th Percentile Velocity $(m/s)$')
    ax2.grid()
    '''
    ax.legend()
    plt.show()


if __name__ == '__main__':

    # initialize logger
    init_logger(__file__)

    #get velocity raster from final timestep for each tuflow run
    '''
    vel_ras_list = get_all_files('G:\\Ken_RB\\c01\\Tuflow\\results\\003\\', suffix='V_002_00.flt')+get_all_files('G:\\Ken_RB\\c01\\Tuflow\\results\\003\\', suffix='V_002_00.hdr')
    # copy rasters to folder for analysis
    for vel_ras in vel_ras_list:
        shutil.copyfile(vel_ras, 'G:\\Ken_RB\\c01\\mu_vel\\' + os.path.basename(vel_ras).replace('.', 'pt', 1).replace('RB_', '').replace('_003_V_002_00', ''))
    '''
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
