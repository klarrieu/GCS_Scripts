import arcpy
arcpy.env.overwriteOutput = True
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from file_functions import *
from Tkinter import *

# add DEM resolution and units parameter

def hyps_analysis(detrended_DEM, bin_size, clipping_polygon='', show_plots=True, save_plots=True):
    '''Creates hypsograph CDF and PDF'''
    check_use(detrended_DEM)
    # take Detrended DEM and clip if clipping polygon given
    dem = detrended_DEM
    if clipping_polygon != '':
        check_use([clipping_polygon, detrended_DEM.replace('.tif', '_hypsclip.tif')])
        logging.info('Clipping DEM...')
        dem = arcpy.Clip_management(dem,
                                    '#',
                                    dem.replace('.tif', '_hypsclip.tif'),
                                    clipping_polygon,
                                    '0',
                                    'ClippingGeometry'
                                    )
        logging.info('OK')

    # find min and max Zd and create bin values
    min_zd = float(arcpy.GetRasterProperties_management(dem, 'MINIMUM').getOutput(0))
    max_zd = float(arcpy.GetRasterProperties_management(dem, 'MAXIMUM').getOutput(0))
    bin_vals = np.arange(min_zd, max_zd, bin_size)

    logging.info('Getting area for each Zd bin...')
    areas = []
    for bin in bin_vals:
        # raster from dem with value of 1 if values <= bin value
        rast = arcpy.sa.Con(dem, 1, '', 'VALUE <= %s' % bin)
        # convert to polygon
        poly = arcpy.RasterToPolygon_conversion(rast,
                                            detrended_DEM.replace('.tif','_hyps.shp'),
                                            'NO_SIMPLIFY'
                                            )
        # get area
        arcpy.AddGeometryAttributes_management(poly,
                                               'AREA'
                                               )
        att_table = arcpy.TableToExcel_conversion(poly,
                                                  str(poly).replace('.shp', '_attribute_table.xls')
                                                  )
        att_df = pd.read_excel(str(att_table))
        area = sum(att_df['POLY_AREA'].tolist())
        areas.append(area)
    logging.info('OK')

    # make CDF plot
    fig = plt.figure()
    cdf = plt.subplot(2, 1, 1)
    plt.title(r'$CDF$')
    plt.xlabel(r'$Z_d$' + ' ' + r'$(m)$')
    plt.ylabel(r'$Area$' + ' ' + r'$(m^2)$')
    plt.grid()
    plt.plot(bin_vals, areas, marker='o')

    # make PDF plot
    pdf = plt.subplot(2, 1, 2)
    plt.title(r'$PDF$')
    plt.xlabel(r'$Z_d$' + ' ' + r'$(m)$')
    plt.ylabel(r'$\Delta Area$' + ' ' + r'$(m^2)$')
    plt.grid()
    delta_as = [a2 - a1 for a1, a2 in zip(areas, areas[1:])]
    plt.plot(bin_vals[1:], delta_as, marker='o')

    if save_plots == True:
        fig.savefig('hypsograph.png')

    if show_plots == True:
        plt.show()

    return [cdf, pdf]


if __name__ == '__main__':

    # initialize logger
    init_logger(__file__)

    # make the GUI window
    root = Tk()
    root.wm_title('Hypsograph Analysis')

    # specify relevant directories/files

    L1 = Label(root, text='Detrended DEM:')
    L1.grid(sticky=E, row=1, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=1, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='file', ftypes=[('Raster', '*.tif'),
                                                                                             ('All files', '*')]
                                                            )
                )

    b1.grid(sticky=W, row=1, column=3)

    L2 = Label(root, text='Clipping Polygon (optional):')
    L2.grid(sticky=E, row=2, column=1)
    E2 = Entry(root, bd=5)
    E2.grid(row=2, column=2)
    b2 = Button(root, text='Browse', command=lambda: browse(root, E2, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                              ('All files', '*')]
                                                            )
                )
    b2.grid(sticky=W, row=2, column=3)

    L3 = Label(root, text='Bin size')
    L3.grid(sticky=E, row=8, column=1)
    E3 = Entry(root, bd=5)
    E3.grid(row=8, column=2)

    b = Button(root, text='   Run    ', command=lambda: hyps_analysis(detrended_DEM=E1.get(),
                                                                      clipping_polygon=E2.get(),
                                                                      bin_size=float(E3.get())
                                                                      )
               )
    b.grid(sticky=W, row=11, column=2)
    root.grid_rowconfigure(11, minsize=80)

    root.mainloop()