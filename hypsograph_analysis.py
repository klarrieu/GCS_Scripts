import arcpy
import numpy as np
from Tkinter import *

def hyps_analysis(detrended_DEM, bin_size, clipping_polygon='', make_plots=True)

    # take Detrended DEM and clip if clipping polygon given
    dem = detrended_DEM
    if clipping_polygon != '':
        dem = arcpy.Clip_analysis(detrended_DEM,
                                  clipping_polygon,
                                  detrended_DEM.replace('.tif','_hypsclip.tif')
                                  )

    # find min and max Zd and create bin values
    min_zd = arcpy.GetRasterProperties_management(dem, 'MINIMUM')
    max_zd = arcpy.GetRasterProperties_management(dem, 'MAXIMUM')
    bin_vals = np.linspace(min_zd, max_zd, bin_size)

    # Get area of raster with values up to each bin value


    # remove upstream and downstream sections

    # make CDF plot

    # make PDF plot

    return


if __name__ == '__main__':

    # make the GUI window
    root = Tk()
    root.wm_title('Hypsograph Analysis')

    # specify relevant directories/files

    L1 = Label(root, text='Detrended DEM:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='file', ftypes=[('Raster', '*.tif'),
                                                                                             ('All files', '*')]
                                                            )
                )

    b2.grid(sticky=W, row=1, column=3)

    L3 = Label(root, text='Clipping Polygon (optional):')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.grid(row=2, column=2)
    b3 = Button(root, text='Browse', command=lambda: browse(root, E3, select='files', ftypes=[('Shapefile', '*.shp'),
                                                                                              ('All files', '*')]
                                                            )
                )
    b3.grid(sticky=W, row=2, column=3)

    L4 = Label(root, text='Bin size')
    L4.grid(sticky=E, row=8, column=1)
    E4 = Entry(root, bd=5)
    E4.grid(row=8, column=2)

    L5 = Label(root, text='Remove Upstream Length:')
    L5.grid(sticky=E, row=9, column=1)
    E5 = Entry(root, bd=5)
    E5.insert(END, '0')
    E5.grid(row=9, column=2)

    L6 = Label(root, text='Remove Downstream Length')
    L6.grid(sticky=E, row=10, column=1)
    E6 = Entry(root, bd=5)
    E6.insert(END, '0')
    E6.grid(row=10, column=2)

    b = Button(root, text='   Run    ', command=lambda: hyps_analysis(station_lines=E1.get(),
                                                                             detrended_DEM=E2.get(),
                                                                             wetted_polygons_list=list(
                                                                                 root.tk.splitlist(E3.get())),
                                                                             buffer_size=E4.get(),
                                                                             rm_up_length = float(E5.get()),
                                                                             rm_down_length = float(E6.get())
                                                                             )
               )
    b.grid(sticky=W, row=11, column=2)
    root.grid_rowconfigure(11, minsize=80)

    root.mainloop()