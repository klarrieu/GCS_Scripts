import arcpy
import arcpy.sa
from arcpy import env
import os
from file_functions import *
import textwrap
import pandas as pd
import logging

arcpy.env.overwriteOutput = True


@err_info
@spatial_license
def extract_channel_data(station_lines, detrended_DEM, wetted_rasters_list, buffer_size='', rm_up_length=0,
                         rm_down_length=0, reach_breaks=''):
    '''Creates wetted polygon XSs for each wetted polygon, then extracts W,Z and other attributes

    Args:
        station_lines: shapefile containing station XS lines
        detrended_DEM: the detrended DEM raster
        wetted_rasters_list: a list containing wetted polygon/raster filenames
        buffer_size: size of rectangular buffer around each XS. Default is 1/2 XS spacing
    
    Returns:
        polygon_XS: 
        table for each wetted polygon containing width, detrended bed elevation, mean velocity, *reach number
    '''

    check_use([station_lines, detrended_DEM] + wetted_rasters_list)
    check_use(str(station_lines).replace('.shp', '_buffered.shp'))
    check_use([name.replace('.shp', '_XS.shp') for name in wetted_rasters_list])
    check_use([name.replace('.shp', '_XS_z_table.dbf') for name in wetted_rasters_list])
    check_use([name.replace('.shp', '_XS_z_table.xls') for name in wetted_rasters_list])
    check_use([name.replace('.shp', '_XS_attribute_table.xls') for name in wetted_rasters_list])
    check_use([name.replace('.shp', '_XS_joined_table.csv') for name in wetted_rasters_list])

    tables = []

    logging.info('Extracting channel data...')

    if buffer_size == '':

        # get spacing from station lines
        rows = []
        with arcpy.da.UpdateCursor(station_lines, ['dist_down']) as cursor:
            while len(rows) < 2:
                for row in cursor:
                    rows.append(row[0])
        spacing = abs(rows[1] - rows[0])

        # set XS buffer size to 1/2 spacing
        buffer_size = spacing * 1.0 / 2
        logging.info('Using buffer size of %f' % buffer_size)

    else:
        buffer_size = float(buffer_size)

    # create rectangular buffer around station lines
    try:
        logging.info('Buffering XS lines...')
        buff = arcpy.Buffer_analysis(station_lines,
                                     str(station_lines).replace('.shp', '_buffered.shp'),
                                     buffer_size,
                                     line_end_type='FLAT'
                                     )
        logging.info('OK')
        for wetted_raster in wetted_rasters_list:

            # if given an .flt raster (e.g. from Tuflow output), convert it to a polygon
            if wetted_raster.endswith('.flt'):
                logging.info('Converting flt to polygon: %s' % wetted_raster)
                wetted_polygon = flt_to_poly(wetted_raster)
                logging.info('OK')
            # clip rectangles to extent of each wetted polygon
            logging.info('Clipping XSs to wetted area...')
            xs = arcpy.Clip_analysis(buff,
                                     wetted_polygon,
                                     str(wetted_polygon).replace('.shp', '_XS.shp')
                                     )
            logging.info('OK')
            # add reach number attribute (all 1 if reach_breaks='') (short integer type)

            arcpy.AddField_management(xs, 'Reach', 'SHORT')
            if reach_breaks != '':
                codeblock = textwrap.dedent('''
                def get_reach(dist_down):
                    reach_num = 1
                    for brk in %s:
                        if dist_down < brk:
                            return reach_num
                        else:
                            reach_num+=1
                    return reach_num
                ''' % reach_breaks)
            else:
                codeblock = textwrap.dedent('''
                def get_reach(dist_down):
                    return 1
                ''')
            logging.info('Calculating cross-section averaged widths...')
            arcpy.CalculateField_management(xs, 'Reach',
                                            'get_reach(!dist_down!)',
                                            'PYTHON',
                                            codeblock)

            # create area and width attributes for each set of xs's
            arcpy.AddGeometryAttributes_management(xs,
                                                   'AREA'
                                                   )
            arcpy.AddField_management(xs,
                                      'W',
                                      'FLOAT'
                                      )
            arcpy.CalculateField_management(xs,
                                            'W',
                                            '!POLY_AREA!*1.0/(2*%f)' % buffer_size,
                                            'PYTHON'
                                            )
            logging.info('OK')
            # make table of detrended DEM stats for each xs
            logging.info('Calculating cross-section averaged detrended bed elevations...')
            z_table = arcpy.sa.ZonalStatisticsAsTable(xs,
                                                      'FID',
                                                      detrended_DEM,
                                                      str(xs).replace('.shp', '_z_table.dbf'),
                                                      statistics_type='MEAN'
                                                      )
            # convert z_table from dbf to xls
            z_table = arcpy.TableToExcel_conversion(z_table,
                                                    str(z_table).replace('.dbf', '.xls')
                                                    )

            logging.info('OK')
            # granted we're given velocity rasters, let's calculate mean XS velocities bud
            # also assuming we have depth rasters with identical file names as velocity rasters but with _d suffix
            logging.info('Getting mean cross-section velocities...')

            depth_raster = wetted_raster.replace('.flt', '_d.flt')

            # (1) velocity * depth
            vd = arcpy.Raster(wetted_raster) * arcpy.Raster(depth_raster)

            # (2) sum of depths in each xs poly
            d_sums = arcpy.sa.ZonalStatistics(xs,
                                              'FID',
                                              depth_raster,
                                              statistics_type='SUM'
                                              )

            # (3) raster (1) divided by raster (2) yields mean velocity when summed over each xs
            v_weighted = vd / d_sums

            v_table = arcpy.sa.ZonalStatisticsAsTable(xs,
                                                      'FID',
                                                      v_weighted,
                                                      str(xs).replace('.shp', '_v_table.dbf'),
                                                      statistics_type='SUM'
                                                      )
            v_table = arcpy.TableToExcel_conversion(v_table,
                                                    str(v_table).replace('.dbf', '.xls')
                                                    )
            logging.info('OK')

            # export attribute table for analysis
            logging.info('Exporting attributes to spreadsheet...')
            att_table = arcpy.TableToExcel_conversion(xs,
                                                      str(xs).replace('.shp', '_attribute_table.xls')
                                                      )

            # join attribute and z tables on FID
            att_df = pd.read_excel(str(att_table))
            z_df = pd.read_excel(str(z_table))
            z_df = z_df.rename(columns={'MEAN': 'Z'})
            joined_df = att_df.merge(z_df,
                                     'outer',
                                     left_on='FID',
                                     right_on='FID_'
                                     )

            # join attributes and z joined table to v table on FID
            v_df = pd.read_excel(str(v_table))
            v_df = v_df.rename(columns={'SUM': 'V'})
            # dropping these columns because they will be re-added by next join
            # this also guarantees areas are at velocity raster resolution instead of DEM resolution
            joined_df = joined_df.drop(['OID', 'COUNT', 'AREA'], axis=1)
            joined_df = joined_df.merge(v_df,
                                        'outer',
                                        left_on='FID',
                                        right_on='FID_'
                                        )
            logging.info('OK')

            # sort the table so data is longitudinally sequential
            if 'dist_down' in joined_df.columns.tolist():
                joined_df = joined_df.sort_values(by=['dist_down'])
                joined_df = joined_df[(joined_df['dist_down'] >= rm_up_length)]
                joined_df = joined_df[joined_df['dist_down'] <= (joined_df['dist_down'].tolist()[-1] - rm_down_length)]
            else:
                logging.warning('No station lines attribute field named dist_down: sorting by original FID for station lines. Ensure output data is longitudinally sequential.')
                if rm_up_length != 0 or rm_down_length != 0:
                    logging.warning('Cannot remove up/downstream sections.')
                joined_df = joined_df.sort_values(by=['ORIG_FID'])

            out_table = str(xs).replace('.shp', '_joined_table.csv')
            joined_df.to_csv(out_table, index=False)
            logging.info('output table: %s' % out_table)

            tables.append(out_table)

            # delete intermediate files
            logging.info('Deleting intermediate files...')
            del_suffixes = ['_z_table.dbf', '_v_table.dbf', '_z_table.xls', '_v_table.xls', '_attribute_table.xls']
            del_files = [str(xs).replace('.shp', suffix) for suffix in del_suffixes]
            for f in del_files:
                arcpy.Delete_management(f)
            logging.info('OK')

    except Exception, e:
        logging.exception(e)
        raise Exception(e)

    logging.info('Finished.')
    return tables


if __name__ == '__main__':
    '''
    #test data
    station_lines = r'G:\Kenny_Rainbow_Basin\DEM\stations_25.shp'
    detrended_DEM = r'G:\Kenny_Rainbow_Basin\DEM\detrending\0_breaks\detrended_DEM.tif'
    wetted_polygons_list = [r'G:\Kenny_Rainbow_Basin\DEM\detrending\0_breaks\wetted_polygon_%i.shp'%i for i in range(1,5)]
    buffer_size = 2.5
    tables = extract_channel_data(station_lines, detrended_DEM, wetted_polygons_list, buffer_size)
    '''

    # initialize logger
    init_logger(__file__)

    # make the GUI window
    root = Tk()
    root.wm_title('Extract Channel Dimensions Data')

    # specify relevant directories/files

    L1 = Label(root, text='Station Lines:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='Detrended DEM:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.grid(row=1, column=2)
    b2 = Button(root, text='Browse', command=lambda: browse(root, E2, select='file', ftypes=[('Raster', '*.tif'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b2.grid(sticky=W, row=1, column=3)

    L3 = Label(root, text='Wetted Polygons:')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.grid(row=2, column=2)
    b3 = Button(root, text='Browse', command=lambda: browse(root, E3, select='files', ftypes=[('Shapefile', '*.shp'),
                                                                                              ('Float Raster', '*.flt'),
                                                                                              ('All files', '*')]
                                                            )
                )
    b3.grid(sticky=W, row=2, column=3)

    L4 = Label(root, text='XS Buffer Size:')
    L4.grid(sticky=E, row=8, column=1)
    E4 = Entry(root, bd=5)
    E4.grid(row=8, column=2)

    L5 = Label(root, text='Remove Upstream Length:')
    L5.grid(sticky=E, row=9, column=1)
    E5 = Entry(root, bd=5)
    E5.insert(END, '0')
    E5.grid(row=9, column=2)

    L6 = Label(root, text='Remove Downstream Length:')
    L6.grid(sticky=E, row=10, column=1)
    E6 = Entry(root, bd=5)
    E6.insert(END, '0')
    E6.grid(row=10, column=2)

    L7 = Label(root, text='Reach Breaks (distance downstream):')
    L7.grid(sticky=E, row=11, column=1)
    E7 = Entry(root, bd=5)
    E7.grid(row=11, column=2)

    b = Button(root, text='   Run    ',
               command=lambda: extract_channel_data(station_lines=E1.get(),
                                                    detrended_DEM=E2.get(),
                                                    wetted_rasters_list=list(root.tk.splitlist(E3.get())),
                                                    buffer_size=E4.get(),
                                                    rm_up_length=float(E5.get()),
                                                    rm_down_length=float(E6.get()),
                                                    reach_breaks=map(int, E7.get().split(',')) if E7.get() != '' else ''
                                                    )
               )
    b.grid(sticky=W, row=12, column=2)
    root.grid_rowconfigure(12, minsize=80)

    root.mainloop()
