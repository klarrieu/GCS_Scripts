import arcpy
arcpy.env.overwriteOutput = True
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from file_functions import *
from Tkinter import *

# input tuflow output rasters, rating curve params

def stage_hyps_analysis(flow_polys, a, b, show_plots=True, save_plots=True):
    '''Creates hypsograph CDF and PDF'''
    logging.info('Getting area for each stage...')
    areas = []
    discharges = []

    for stage_poly in flow_polys:
        poly = stage_poly
        if stage_poly.endswith('.flt'):
            # convert raster to polygon if necessary
            poly = flt_to_poly(stage_poly)

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

        discharge = float(os.path.basename(poly).replace('.shp', '').replace('pt','.'))
        discharges.append(discharge)
        print(str(discharge) + ' cms')
        print(area)

    logging.info('OK')

    stages = [(Q*1.0/a)**(1.0/b) for Q in discharges]

    # sort by ascending stage
    sorted_pairs = sorted(zip(stages, areas))
    stages = [pair[0] for pair in sorted_pairs]
    areas = [pair[1] for pair in sorted_pairs]

    # make CDF plot
    fig = plt.figure()
    cdf = plt.subplot(2, 1, 1)
    plt.title(r'$CDF$')
    plt.xlabel(r'$Stage$' + ' ' + r'$(m)$')
    plt.ylabel(r'$Area$' + ' ' + r'$(m^2)$')
    plt.grid()
    plt.plot(stages, areas, marker='o')

    # make PDF plot
    pdf = plt.subplot(2, 1, 2)
    plt.title(r'$PDF$')
    plt.xlabel(r'$Stage$' + ' ' + r'$(m)$')
    plt.ylabel(r'$\Delta Area$' + ' ' + r'$(m^2)$')
    plt.grid()
    delta_as = [a2 - a1 for a1, a2 in zip(areas, areas[1:])]
    plt.plot(stages[1:], delta_as, marker='o')

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

    L1 = Label(root, text='Flow Polygons:')
    L1.grid(sticky=E, row=1, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=1, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='files', ftypes=[('Float Raster', '*.flt'),
                                                                                              ('Shapefiles', '*.shp'),
                                                                                              ('All files', '*')]
                                                            )
                )

    b1.grid(sticky=W, row=1, column=3)

    L2 = Label(root, text='Rating curve a:')
    L2.grid(sticky=E, row=2, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, '18.001499')
    E2.grid(row=2, column=2)

    L3 = Label(root, text='Rating curve b:')
    L3.grid(sticky=E, row=8, column=1)
    E3 = Entry(root, bd=5)
    E3.insert(END, '2.489')
    E3.grid(row=8, column=2)

    b = Button(root, text='   Run    ', command=lambda: stage_hyps_analysis(flow_polys=list(root.tk.splitlist(E1.get())),
                                                                            a=float(E2.get()),
                                                                            b=float(E3.get())
                                                                            )
               )
    b.grid(sticky=W, row=11, column=2)
    root.grid_rowconfigure(11, minsize=80)

    root.mainloop()