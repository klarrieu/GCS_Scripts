'''
#input: x,y,z coordinate table from centerline-station intersections
#output: plot of fitted elevation profile
         plot of residual elevations after fitting
         plot of standardized residual elevations
         detrended DEM
'''

import arcpy
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
import matplotlib.pyplot as plt
from file_functions import *
import logging


def station_coords(centerline, station_lines, DEM):
    '''Computes points of intersection between centerline and station lines.
        Outputs an .xls file containing intersection XYZ coordinates.'''
    outdir = os.path.dirname(DEM) + '\\'
    check_use([centerline,
               station_lines,
               DEM,
               outdir + 'intersection_pts.shp',
               outdir + 'intersection_pts_sp.shp',
               outdir + 'intersection_coords.shp',
               outdir + 'intersection_coords.xls'
               ])

    # get intersection multipoint feature
    logging.info('Getting centerline/station line intersection points...')
    intersections = arcpy.Intersect_analysis([centerline, station_lines],
                                             outdir + 'intersection_pts.shp',
                                             output_type='POINT'
                                             )
    # convert multipart feature to singlepart features
    intersections = arcpy.MultipartToSinglepart_management(intersections,
                                                           outdir + 'intersection_pts_sp.shp')
    logging.info('OK')

    logging.info('Extracting XYZ coordinates at intersection points...')
    # extract Z coordinates to points attribute table for each point from DEM
    intersections = arcpy.sa.ExtractValuesToPoints(intersections,
                                                   DEM.split('.aux')[0],
                                                   outdir + 'intersection_coords.shp'
                                                   )
    # add XY coordinates to attribute table
    intersections = arcpy.AddXY_management(outdir + 'intersection_coords.shp')
    logging.info('OK')
    logging.info('Exporting intersection attribute table to spreadsheet...')
    # export attribute table to spreadsheet
    table = arcpy.TableToExcel_conversion(outdir + 'intersection_coords.shp',
                                          outdir + 'intersection_coords.xls'
                                          )
    logging.info('OK')

    # delete intermediate files
    logging.info('Deleting intermediate files...')
    del_files = [outdir + 'intersection_pts.shp', outdir + 'intersection_pts_sp.shp', outdir + 'intersection_coords.shp']
    for f in del_files:
        arcpy.Delete_management(f)
    logging.info('OK')

    return table.getOutput(0)


def trend_fit(intersection_coords, station_lines, slope_break_indices=[], regression='linear', make_plot=False):
    '''Fits longitudinal elevation profile to piecewise linear or quadratic function.
        Slope breaks are made at defined indices.'''

    # get units
    units = arcpy.Describe(station_lines).spatialReference.linearUnitName
    if units == 'Meter':
        units = 'm'
    elif units == 'US Foot':
        units = 'ft'

    ###############################
    # read in coordinates to dataframe

    df = pd.read_excel(intersection_coords)
    if 'dist_down' in df.columns.tolist():
        coords = df[['ORIG_FID', 'dist_down', 'POINT_X', 'POINT_Y', 'RASTERVALU']]
    else:
        coords = df[['ORIG_FID', 'LOCATION', 'POINT_X', 'POINT_Y', 'RASTERVALU']]

    # rename columns
    coords.columns = ['station', 'distance downstream', 'x', 'y', 'z']
    # get spacing between stations
    spacing = abs(coords['distance downstream'][1] - coords['distance downstream'][0])
    # index and sort values going downstream
    coords.index = map(int, coords['distance downstream'] / spacing)
    coords = coords.sort_values(by=['distance downstream'])

    # station is numbered going downstream
    station = coords.index.tolist()
    # distance downstream
    dist = coords['distance downstream'].tolist()
    # coordinates
    x = coords.x.tolist()
    y = coords.y.tolist()
    z = coords.z.tolist()

    ###############################
    # use slope_breaks to make a regression for each reach
    slope_break_dist = [spacing * x_i for x_i in slope_break_indices]

    # fit_params contains lists of [slope, intercept] pairs for each fitted line segment
    fit_params = []
    d_sections = split_list(dist, slope_break_indices)
    z_sections = split_list(z, slope_break_indices)

    for i, (d_section, z_section) in enumerate(zip(d_sections, z_sections)):
        if regression == 'linear':
            m, b = np.polyfit(d_section, z_section, deg=1)
            fit_params.append([m, b])
            logging.info('Reach %i length: %i%s' % (i + 1, int(spacing * (len(d_section) - 1)), units))
            logging.info('Reach %i slope: %f' % (i + 1, abs(m)))
        elif regression == 'quadratic':
            a, b, c = np.polyfit(d_section, z_section, deg=2)
            fit_params.append([a, b, c])
            quadratic = lambda x_i: a*x_i**2 + b*x_i + c
            mean_slope = (quadratic(d_section[-1]) - quadratic(d_section[0]))/(d_section[-1] - d_section[0])
            logging.info('Reach %i length: %i%s' % (i + 1, int(spacing * (len(d_section) - 1)), units))
            logging.info('Reach %i mean slope: %f' % (i + 1, abs(mean_slope)))

    ###############################
    # calculate residuals
    z_fits = []
    for i, d_section in enumerate(d_sections):
        if regression == 'linear':
            m, b = fit_params[i]
            z_fits.extend([m * x_i + b for x_i in d_section])
        elif regression == 'quadratic':
            a, b, c = fit_params[i]
            z_fits.extend([a*x_i**2 + b*x_i + c for x_i in d_section])
    z_res = [z_val - z_fit for z_val, z_fit in zip(z, z_fits)]

    # calculate r squared goodness of fit
    ss_res = np.sum(np.asarray(z_res) ** 2)
    ss_tot = np.sum((np.asarray(z) - np.mean(z)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    logging.info('r^2: %f' % r_squared)

    ###############################
    # standardize the residuals

    '''
    #using mean and std dev for each reach
    zs = []
    #split up z_res into section for each domain
    z_res_sections = split_list(z_res, slope_break_indices)
    for i, z_res_section in enumerate(z_res_sections):
        z_res_section_mean = np.mean(z_res_section)
        z_res_section_sigma = np.std(z_res_section)
        zs.extend([(z_res_section_val - z_res_section_mean)*1.0/z_res_section_sigma for z_res_section_val in z_res_section])
    '''

    # using mean and std dev for entire dataset
    z_mean = np.mean(z_res)
    z_sigma = np.std(z_res)
    zs = [(z_val - z_mean) * 1.0 / z_sigma for z_val in z_res]

    ###############################
    # save fitted coordinate data for detrending

    # save a .csv containing x,y,z_fit table
    xyz_fit = pd.DataFrame({'x': x, 'y': y, 'z_fit': z_fits})
    xyz_fit.to_csv(os.path.dirname(intersection_coords) + '/xyz_fit.csv', index=False)
    logging.info('Saved fitted coordinates to %s' % (os.path.dirname(intersection_coords) + '\\xyz_fit.csv'))

    if make_plot == True:
        ################################
        # plot the elevation profile with fitted line segments
        plt.figure(figsize=(24, 16))
        plot1 = plt.subplot(3, 1, 1)
        plt.title('Longitudinal Bed Elevation')
        plt.xlabel('Distance Downstream (%s)' % units)
        plt.ylabel('Centerline Bed Elevation (%s)' % units)
        plt.grid()
        plt.plot(dist, z, label='longitudinal profile')
        plt.xlim(dist[0], dist[-1])
        for i, d_section in enumerate(d_sections):
            if regression == 'linear':
                m, b = fit_params[i]
                plt.plot(d_section, [m * x_i + b for x_i in d_section], 'r', label='%s regression' % regression)
            elif regression == 'quadratic':
                a, b, c = fit_params[i]
                plt.plot(d_section, [a*x_i**2 + b*x_i + c for x_i in d_section], 'r', label='%s regression' % regression)
        for x_val in slope_break_dist:
            plt.axvline(x=x_val - 0.5, linestyle='--')
        plt.legend()

        ###############################
        # plot the residuals

        plot2 = plt.subplot(3, 1, 2)
        plt.title('Elevation Residuals')
        plt.xlabel('Distance Downstream (%s)' % units)
        plt.ylabel('Elevation Residuals (%s)' % units)
        plt.grid()
        plt.axhline(y=0, linestyle='--', color='black')
        plt.plot(dist, z_res)
        plt.xlim(dist[0], dist[-1])
        for x_val in slope_break_dist:
            plt.axvline(x=x_val - 0.5, linestyle='--')

        ###############################
        # plot the standardized residuals

        plot3 = plt.subplot(3, 1, 3)
        plt.title('Standardized Elevation Residuals')
        plt.xlabel('Distance Downstream (%s)' % units)
        plt.ylabel(r'Standardized Elevation Residuals ($\sigma$)')
        plt.grid()
        plt.axhline(y=0, linestyle='--', color='black')
        plt.plot(dist, zs)
        for x_val in slope_break_dist:
            plt.axvline(x=x_val - 0.5, linestyle='--')

        ###############################
        # save the plots
        plot_file = os.path.dirname(intersection_coords) + '\\longitudinal_profile.png'
        plt.savefig(plot_file)
        plt.close()
        logging.info('Saved longitudinal profile plot: %s' % plot_file)

    # delete intermediate files

    logging.info('Deleting intermediate files...')
    del_files = [intersection_coords]
    for f in del_files:
        arcpy.Delete_management(f)
    logging.info('OK')

    return os.path.dirname(intersection_coords) + '\\xyz_fit.csv'


def detrend_DEM(fit_table, DEM):
    '''Uses table fitted coordinates to detrend the DEM'''

    outdir = os.path.dirname(fit_table) + '/'
    check_use([fit_table,
               DEM,
               fit_table.replace('.csv', '.lyr'),
               fit_table.replace('.csv', '.shp'),
               outdir + 'thiessen_polygons.shp',
               outdir + 'z_fit_raster.tif',
               outdir + 'detrended_DEM.tif'
               ])

    logging.info('Making points from fitted coordinates...')
    # create feature class from fitted values
    points = arcpy.MakeXYEventLayer_management(fit_table,
                                               out_layer='points_layer',
                                               in_x_field='x',
                                               in_y_field='y',
                                               in_z_field='z_fit',
                                               spatial_reference=arcpy.Describe(DEM.split('.aux')[0]).SpatialReference
                                               )
    points = arcpy.SaveToLayerFile_management('points_layer',
                                              fit_table.replace('.csv', '.lyr')
                                              )
    points = arcpy.FeatureClassToFeatureClass_conversion(points,
                                                         outdir,
                                                         os.path.basename(fit_table).replace('.csv', '.shp'),
                                                         )

    logging.info('OK')
    logging.info('Creating Thiessen polygons...')
    # create thiessen polygons
    thiessen = arcpy.CreateThiessenPolygons_analysis(points,
                                                     outdir + 'thiessen_polygons.shp',
                                                     fields_to_copy='ALL'
                                                     )
    logging.info('OK')

    logging.info('Creating detrended DEM raster...')
    # convert polygons to raster
    z_fit_raster = arcpy.PolygonToRaster_conversion(thiessen,
                                                    'z_fit',
                                                    outdir + 'z_fit_raster.tif',
                                                    cellsize=float(arcpy.GetRasterProperties_management(DEM.split('.aux')[0],
                                                                                                        'CELLSIZEX').getOutput(
                                                        0))
                                                    )

    # subtract from DEM raster
    detrended_DEM = arcpy.Raster(DEM.split('.aux')[0]) - arcpy.Raster(z_fit_raster)
    detrended_DEM.save(outdir + 'detrended_DEM.tif')

    logging.info('OK')

    logging.info('Deleting intermediate files...')
    del_files = [thiessen, fit_table, points, fit_table.replace('.csv', '.lyr'), z_fit_raster]
    for f in del_files:
        arcpy.Delete_management(f)
    logging.info('OK')

    logging.info('Saved detrended DEM: %s' % detrended_DEM)
    return str(detrended_DEM)


def channel_buffer(centerline, buffer_size=50):
    '''Creates a buffer polygon around the centerline.
        buffer_size specifies distance around centerline to make polygon'''
    outdir = os.path.dirname(centerline) + '/'
    check_use([centerline,
               outdir + 'channel_buffer.shp'
               ])
    logging.info('Creating buffer...')
    polygon = arcpy.Buffer_analysis(centerline,
                                    outdir + 'channel_buffer.shp',
                                    buffer_size
                                    )
    logging.info('OK')
    return str(polygon)


def clip_raster(raster, polygon):
    '''Clips raster to extent of polygon'''
    outdir = os.path.dirname(raster) + '/'
    out_name = 'clipped_' + os.path.basename(raster)
    check_use([raster,
               polygon,
               outdir + out_name
               ])
    logging.info('Clipping raster...')
    clipped_raster = arcpy.Clip_management(raster,
                                           '#',
                                           outdir + out_name,
                                           polygon,
                                           '0',
                                           'ClippingGeometry'
                                           )
    logging.info('OK')
    return clipped_raster


@err_info
@spatial_license
def main_det(DEM, centerline, station_lines, slope_break_indices=[], regression='linear'):
    xyz_table = station_coords(centerline, station_lines, DEM)
    xyz_fit_table = trend_fit(xyz_table, station_lines, slope_break_indices=slope_break_indices, regression=regression, make_plot=True)
    det = detrend_DEM(xyz_fit_table, DEM)
    # det = clip_raster(det, channel_shp)

    return det


#######################################
if __name__ == '__main__':
    # input files
    idir = 'G:\\Ken_RB\\c01\\fit_test\\'
    centerline = idir + 'clipped_centerline.shp'
    station_lines = idir + 'stations_100.shp'
    DEM = idir + 'RB_DEM.tif'

    init_logger(__file__)

    xyz_table = station_coords(centerline, station_lines, DEM)
    # xyz_fit_table = trend_fit(xyz_table, station_lines, slope_break_indices=[556], regression='linear', make_plot=True)
    xyz_fit_table = trend_fit(xyz_table, station_lines, slope_break_indices=[], regression='quadratic', make_plot=True)
    detrended_DEM = detrend_DEM(xyz_fit_table, DEM)
    # buff = channel_buffer(centerline, buffer_size = 50)
    # clipped_detrended_DEM = clip_raster(detrended_DEM, buff)
