import arcpy
from arcpy import env
import os
from file_functions import *
import pandas as pd
import logging

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')


def extract_channel_data(station_lines, detrended_DEM, wetted_polygons_list, buffer_size='', rm_up_length=0, rm_down_length=0, reach_breaks=''):
    '''Creates wetted polygon XSs for each wetted polygon, then extracts W,Z and other attributes

    Args:
        station_lines: shapefile containing station XS lines
        detrended_DEM: the detrended DEM raster
        wetted_polygons_list: a list containing wetted polygon filenames
        buffer_size: size of rectangular buffer around each XS. Default is 1/2 XS spacing
    
    Returns:
        polygon_XS: 
        table for each wetted polygon containing width, detrended bed elevation, *reach number
    '''

    check_use([station_lines, detrended_DEM] + wetted_polygons_list)
    check_use(str(station_lines).replace('.shp', '_buffered.shp'))
    check_use([name.replace('.shp', '_XS.shp') for name in wetted_polygons_list])
    check_use([name.replace('.shp', '_XS_z_table.dbf') for name in wetted_polygons_list])
    check_use([name.replace('.shp', '_XS_z_table.xls') for name in wetted_polygons_list])
    check_use([name.replace('.shp', '_XS_attribute_table.xls') for name in wetted_polygons_list])
    check_use([name.replace('.shp', '_XS_joined_table.csv') for name in wetted_polygons_list])

    tables = []

    logging.info('Extracting channel data...')

    if buffer_size == '':

        # get spacing from station lines
        rows = []
        with arcpy.da.UpdateCursor(station_lines, ['ET_STATION']) as cursor:
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
        buff = arcpy.Buffer_analysis(station_lines,
                                     str(station_lines).replace('.shp', '_buffered.shp'),
                                     buffer_size,
                                     line_end_type='FLAT'
                                     )

        for wetted_polygon in wetted_polygons_list:

            # if given an .flt raster (e.g. from Tuflow output), convert it to a polygon
            if wetted_polygon.endswith('.flt'):
                logging.info('Converting flt to polygon: %s' % wetted_polygon)
                wetted_polygon = flt_to_poly(wetted_polygon)

            # clip rectangles to extent of each wetted polygon
            xs = arcpy.Clip_analysis(buff,
                                     wetted_polygon,
                                     str(wetted_polygon).replace('.shp', '_XS.shp')
                                     )

            # add reach number attribute (all 1 if reach_breaks='') (short integer type)
            arcpy.AddField_management(xs, 'Reach', 'SHORT')
            codeblock = '''
            def get_reach(dist_down):
                reach_num = 1
                for break in %s:
                    if dist_down < reach_num:
                        return reach_num
                    else:
                        reach_num+=1
                return reach_num
            ''' % reach_breaks
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

            # make table of detrended DEM stats for each xs
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

            # export attribute table for analysis
            att_table = arcpy.TableToExcel_conversion(xs,
                                                      str(xs).replace('.shp', '_attribute_table.xls')
                                                      )

            # join attribute and z tables on FID
            att_df = pd.read_excel(str(att_table))
            z_df = pd.read_excel(str(z_table))
            z_df = z_df.rename(columns={'MEAN': 'Z'})
            joined_df = att_df.join(z_df,
                                    'FID',
                                    'outer'
                                    )

            # sort the table so data is longitudinally sequential
            if 'dist_down' in joined_df.columns.tolist():
                joined_df = joined_df.sort_values(by=['dist_down'])
                joined_df = joined_df[(joined_df['dist_down'] >= rm_up_length)]
                joined_df = joined_df[joined_df['dist_down'] <= (joined_df['dist_down'].tolist()[-1] - rm_down_length)]
            elif 'ET_STATION' in joined_df.columns.tolist():
                joined_df = joined_df.sort_values(by=['ET_STATION'])
                joined_df = joined_df[(joined_df['ET_STATION'] >= rm_up_length)]
                joined_df = joined_df[joined_df['ET_STATION'] <= (joined_df['ET_STATION'].tolist()[-1] - rm_down_length)]
            else:
                logging.warning(
                    'No station lines attribute field named dist_down or ET_STATION: sorting by original FID for station lines. Ensure output data is longitudinally sequential.')
                if rm_up_length != 0 or rm_down_length != 0:
                    logging.warning('Cannot remove up/downstream sections.')
                joined_df = joined_df.sort_values(by=['ORIG_FID'])

            out_table = str(xs).replace('.shp', '_joined_table.csv')
            joined_df.to_csv(out_table, index=False)
            logging.info('output table: %s' % out_table)

            tables.append(out_table)

    except Exception, e:
        logging.exception(e)
        raise Exception(e)

    logging.info('OK')
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

    L7 = Label(root, text='Reach Breaks:')
    L7.grid(sticky=E, row=11, column=1)
    E7 = Entry(root, bd=5)
    E7.grid(row=11, column=2)

    b = Button(root, text='   Run    ',
               command=lambda: extract_channel_data(station_lines=E1.get(),
                                                    detrended_DEM=E2.get(),
                                                    wetted_polygons_list=list(root.tk.splitlist(E3.get())),
                                                    buffer_size=E4.get(),
                                                    rm_up_length=float(E5.get()),
                                                    rm_down_length=float(E6.get()),
                                                    reach_breaks=map(int, E7.get().split(',')) if E7.get() != '' else ''
                                                    )
               )
    b.grid(sticky=W, row=12, column=2)
    root.grid_rowconfigure(12, minsize=80)

    root.mainloop()

arcpy.CheckInExtension('Spatial')
