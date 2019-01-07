from file_functions import *
init_logger(__file__)
arcpy.env.overwriteOutput = True


def make_thalweg(d_ras, v_ras, source, smooth_distance):
    '''
    Creates thalweg shapefile.
    Calculated by calculating conveyance (depth*velocity^2) raster, inverting, then calculating least-cost path

    Args:
        d_ras: depth raster
        v_ras: velocity raster
        source: polygon shapefile for beginning area of thalweg
        smooth_distance: tolerance for PAEK smoothing algorithm

    Returns:
        thalweg shapefile (saved to directory of d_ras)
    '''


    odir = os.path.dirname(d_ras)+'\\'

    depth = arcpy.Raster(d_ras)
    velocity = arcpy.Raster(v_ras)
    logging.info('Calculating inverse conveyance raster...')
    conveyance = depth*velocity**2
    max_conveyance = float(arcpy.GetRasterProperties_management(conveyance, 'MAXIMUM').getOutput(0))
    inv_conveyance = max_conveyance - conveyance
    inv_conveyance.save(odir + 'inv_conveyance.tif')
    logging.info('OK')
    logging.info('Filling sinks in inverse conveyance raster...')
    fill_inv_conveyance = arcpy.sa.Fill(inv_conveyance)
    fill_inv_conveyance.save(odir + 'fill_inv_conveyance.tif')
    logging.info('OK')
    logging.info('Calculating flow direction...')
    flow_direction = arcpy.sa.FlowDirection(fill_inv_conveyance)
    flow_direction.save(odir + 'flow_dir.tif')
    logging.info('OK')

    # least cost path
    logging.info('Calculating least cost path...')
    lc_ras = arcpy.sa.CostPath(source, fill_inv_conveyance, flow_direction, path_type='BEST_SINGLE', destination_field='Id')
    lc_ras.save(odir+'\\lc_ras.tif')
    logging.info('OK')
    logging.info('Converting path to thalweg polyline...')
    thalweg = arcpy.RasterToPolyline_conversion(lc_ras, odir + 'rough_thalweg.shp', simplify='NO_SIMPLIFY')
    logging.info('OK')

    # smooth line
    logging.info('Smoothing thalweg...')
    if smooth_distance != 0:
        thalweg = arcpy.cartography.SmoothLine(thalweg, odir + 'thalweg.shp', algorithm='PAEK', tolerance=smooth_distance)

    logging.info('OK.')
    logging.info('Finished: %s' % thalweg)

    return thalweg


if __name__ == '__main__':

    root = Tk()
    root.wm_title('Create thalweg')

    L1 = Label(root, text='Depth raster:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse',
                command=lambda: browse(root, E1, select='file', ftypes=[('Float Raster', '*.flt'),
                                                                        ('All files', '*')]
                                       )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='Velocity raster:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E2.grid(row=1, column=2)
    b2 = Button(root, text='Browse',
                command=lambda: browse(root, E2, select='file', ftypes=[('Float Raster', '*.flt'),
                                                                        ('All files', '*')]
                                       )
                )
    b2.grid(sticky=W, row=1, column=3)

    L3 = Label(root, text='Flow Source Polygon:')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E3.grid(row=2, column=2)
    b3 = Button(root, text='Browse',
                command=lambda: browse(root, E3, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                        ('All files', '*')]
                                       )
                )
    b3.grid(sticky=W, row=2, column=3)

    L4 = Label(root, text='Smoothing Distance:')
    L4.grid(sticky=E, row=3, column=1)
    E4 = Entry(root, bd=5)
    E4.insert(END, 20)
    E4.grid(row=3, column=2)


    b = Button(root, text='   Run    ',
               command=lambda: make_thalweg(d_ras=E1.get(), v_ras=E2.get(), source=E3.get(), smooth_distance=float(E4.get()))
               )
    b.grid(sticky=W, row=4, column=2)
    root.grid_rowconfigure(4, minsize=80)

    root.mainloop()
