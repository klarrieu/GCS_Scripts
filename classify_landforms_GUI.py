import os
from file_functions import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import arcpy

arcpy.env.overwriteOutput = True
from Tkinter import *
import logging


def clean_in_table(table, w_field='Width', z_field='Z', v_field='V', dist_field='dist_down'):
    '''Renames columns corresponding to W, Z, and dist_down, if they are not already columns'''
    check_use(table)
    df = pd.read_csv(table)

    for old_name, replacement_name in [(w_field, 'W'), (z_field, 'Z'), (v_field, 'V'), (dist_field, 'dist_down')]:
        if replacement_name not in df.columns.tolist():
            if old_name in df.columns.tolist():
                df.rename(columns={old_name: replacement_name})
                logging.info('Renamed %s to %s in %s' % (old_name, replacement_name, table))
            else:
                logging.exception('Cannot find column named %s or %s in %s' % (old_name, replacement_name, table))

    df.to_csv(table, index=False)

    return df


def standardize(table, fields):
    '''Makes standardized version of field in csv table by subtracting each value by mean and dividing by standard deviation.'''
    check_use(table)
    df = pd.read_csv(table)

    if type(fields) == list:
        for f in fields:
            df = standardize(table, f)
        return df

    new_field = fields + '_s'
    df[new_field] = (df[fields] - np.mean(df[fields])) * 1.0 / np.std(df[fields])
    df.to_csv(table, index=False)

    return df


def std_covar_series(table, field_1, field_2):
    '''Computes the covariance series between standardized variables field_1 and field_2 (GCS series)'''
    check_use(table)
    df = pd.read_csv(table)
    df['%s_%s' % (field_1, field_2)] = df[field_1] * df[field_2]
    df.to_csv(table, index=False)

    return df


def landforms(table, zs_field='Z_s', ws_field='W_s', na=-9999):
    '''Classifies each row by corresponding landform type:
        oversized, nozzle, constricted pool, wide bar, normal channel
        Adds columns to input table'''
    check_use(table)
    df = pd.read_csv(table)

    df['normal'] = [zs * ws if (abs(zs * ws) <= 0.5) else na for zs, ws in zip(df[zs_field], df[ws_field])]
    df['wide_bar'] = [zs * ws if (abs(zs * ws) > 0.5 and zs > 0 and ws > 0) else na for zs, ws in
                      zip(df[zs_field], df[ws_field])]
    df['const_pool'] = [zs * ws if (abs(zs * ws) > 0.5 and zs < 0 and ws < 0) else na for zs, ws in
                        zip(df[zs_field], df[ws_field])]
    df['nozzle'] = [zs * ws if (abs(zs * ws) > 0.5 and zs > 0 and ws < 0) else na for zs, ws in
                    zip(df[zs_field], df[ws_field])]
    df['oversized'] = [zs * ws if (abs(zs * ws) > 0.5 and zs < 0 and ws > 0) else na for zs, ws in
                       zip(df[zs_field], df[ws_field])]

    df['code'] = [-2 if df['oversized'][i] != na
                  else -1 if df['const_pool'][i] != na
                  else 0 if df['normal'][i] != na
                  else 1 if df['wide_bar'][i] != na
                  else 2 if df['nozzle'][i] != na
                  else na
                  for i in range(len(df))
                  ]

    df.to_csv(table, index=False)
    return df


def landform_polygons(table, wetted_XS_polygon):
    '''Adds landform code and Cov(Z_s, W_s) from table to attributes for corresponding wetted polygon'''
    check_use([table, wetted_XS_polygon])
    df = pd.read_csv(table)
    join_on = 'FID'

    logging.info('Classifying by landform: %s...' % wetted_XS_polygon)

    # add column names from table as attribute fields for wetted_polygon
    for col in df.columns.tolist():
        current_fields = [field.name for field in arcpy.ListFields(wetted_XS_polygon)]
        if (col != join_on) and (col not in current_fields) and (col + '_' not in current_fields):
            if col == 'code':  # make 'code' an integer field so it can be referenced by zonal stats
                arcpy.AddField_management(wetted_XS_polygon,
                                          col,
                                          'SHORT'
                                          )
            else:
                arcpy.AddField_management(wetted_XS_polygon,
                                          col,
                                          'FLOAT'
                                          )
    # fill each new column with values using join_on column name
    rows = arcpy.UpdateCursor(wetted_XS_polygon)
    for row in rows:
        join_on_val = row.getValue(join_on)
        for col in df.columns.tolist():
            if col != join_on:
                try:
                    new_val = float(df.loc[df[join_on] == join_on_val][col].tolist()[0])
                    if col == 'OID':
                        row.setValue('OID_', new_val)
                    else:
                        row.setValue(col, new_val)
                # index will be out of range if we removed upstream and downstream sections
                except IndexError:
                    # set 'code' to -9999 instead of 0 for removed sections so it is not mistaken for normal channel
                    row.setValue('code', -9999)
                    pass

        rows.updateRow(row)

    logging.info('OK')
    a = arcpy.da.FeatureClassToNumPyArray(wetted_XS_polygon, '*')
    # remove 'Shape' list from the array so we can put it in a dataframe
    df = pd.DataFrame(a[[x for x in a.dtype.names if x != 'Shape']])
    return df


# *** could generalize to do all GCS plots
def GCS_plot(table, units='m'):
    '''Takes a table with longitudinally sequential W_s, Z_s, and W_s_Z_s as columns and plots them'''
    check_use(table)
    df = pd.read_csv(table)

    if 'dist_down' in df.columns.tolist():
        x = df.dist_down.tolist()
    else:
        logging.warning('No column named dist_down. Plotting with index along x axis.')
        x = df.index.tolist()

    try:
        ws = df.W_s.tolist()
        zs = df.Z_s.tolist()
        zs_ws = df.W_s_Z_s.tolist()
    except Exception, e:
        logging.exception(e)
        logging.info('Could not find columns Z_s, W_s, and W_s_Z_s with GCS data.')
        raise Exception(e)

    fig = plt.figure(figsize=(24, 12))
    fig.suptitle(table)

    # Ws plot
    plot1 = plt.subplot(3, 1, 1)
    plt.xlabel('Distance downstream (%s)' % units)
    plt.ylabel('Zs')
    plt.grid()
    plt.axhline(y=0, linestyle='--', color='black')
    plt.plot(x, zs)

    # Zs plot
    plot2 = plt.subplot(3, 1, 2)
    plt.xlabel('Distance downstream (%s)' % units)
    plt.ylabel('Ws')
    plt.grid()
    plt.axhline(y=0, linestyle='--', color='black')
    plt.plot(x, ws)

    # Zs*Ws plot
    plot3 = plt.subplot(3, 1, 3)
    plt.xlabel('Distance downstream (%s)' % units)
    plt.ylabel('Zs*Ws')
    plt.grid()
    plt.axhline(y=0, linestyle='--', color='black')
    plt.plot(x, zs_ws)

    plt.show()


@err_info
def main_classify_landforms(tables, w_field, z_field, v_field, dist_field, make_plots=False):
    '''Classifies river segments as normal, wide bar, constricted pool, oversized, or nozzle

    Args:
        tables: a list of attribute table filenames for each set of wetted polygon rectangular XS's
        w_field: name of the attribute table field corresponding to width
        z_field: name of the attribute table field corresponding to detrended bed elevation
        v_field: name of the attribute table field corresponding to velocity
        dist_field: name of the attribute table field corresponding to distance downstream
        
    Returns:
        For each input table:
            a .csv containing dist_down, W, Z, W_s, Z_s, W_s*Z_s, and landform classification/code fields
            adds these computed values to attribute tables of corresponding wetted polygon rectangular XS's
    '''
    out_polys = []
    fields = [w_field, z_field, v_field]
    std_fields = [field + '_s' for field in fields]
    std_pairs = list(itertools.combinations(std_fields, 2))

    for table in tables:
        clean_in_table(table, w_field=w_field, z_field=z_field, v_field=v_field, dist_field=dist_field)
        standardize(table, fields=fields)
        for std_pair in std_pairs:
            std_covar_series(table, std_pair[0], std_pair[1])
        df = landforms(table)
        df = landform_polygons(table, table.replace('_joined_table.csv', '.shp'))
        out_polys.append(table.replace('_joined_table.csv', '.shp'))
        if make_plots == True:
            GCS_plot(table)

    logging.info('Finished.')

    return out_polys


if __name__ == '__main__':

    # initialize logger
    init_logger(__file__)

    # make the GUI window
    root = Tk()
    root.wm_title('Classify Landforms')

    # specify relevant directories/files

    L1 = Label(root, text='Wetted XS Attribute Tables:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse',
                command=lambda: browse(root, E1, select='files', ftypes=[('Comma-delimited text', '.csv'),
                                                                         ('All files', '*')]
                                       )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='Width Field:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, 'W')
    E2.grid(row=1, column=2)

    L3 = Label(root, text='Detrended Elevation Field:')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.insert(END, 'Z')
    E3.grid(row=2, column=2)

    L4 = Label(root, text='Velocity Field:')
    L4.grid(sticky=E, row=3, column=1)
    E4 = Entry(root, bd=5)
    E4.insert(END, 'V')
    E4.grid(row=3, column=2)

    L5 = Label(root, text='Distance Downstream Field:')
    L5.grid(sticky=E, row=4, column=1)
    E5 = Entry(root, bd=5)
    E5.insert(END, 'dist_down')
    E5.grid(row=4, column=2)

    b = Button(root, text='   Run    ',
               command=lambda: main_classify_landforms(tables=list(root.tk.splitlist(E1.get())),
                                                       w_field=E2.get(),
                                                       z_field=E3.get(),
                                                       v_field=E4.get(),
                                                       dist_field=E5.get()
                                                       )
               )
    b.grid(sticky=W, row=9, column=2)
    root.grid_rowconfigure(9, minsize=80)

    root.mainloop()
