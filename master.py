'''
This combines the scripts for every GCS data processing/analysis step into one compact GUI.

Uses a GCSHandler class to store input data as class attributes, store processing/analysis functions as class methods.

Use an excel notebook to store input data for things like:
-landform classification thresholds?
-reach breaks?
-tuflow discharges?

The central GUI employs a tab for each processing step.

'''
import tkinter as tk
from tkinter import ttk
from file_functions import *
import LiDAR_processing_GUI as lp
import create_centerline_GUI as cc
import create_station_lines as sl
import DEM_Detrending as dd
import extract_channel_dims_GUI as ec
import classify_landforms_GUI as cl
import GCS_analysis as gcsa


class GCS_GUI(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()

        # initialize attribute files
        self.dem = StringVar()
        self.det_dem = StringVar()
        self.centerline = StringVar()
        self.station_lines = StringVar()
        # ***etc. If one of these files is produced when a function runs, or entered by another browse window, update corresponding ttk.Entry widgets accordingly

        # set header
        self.master.title('River GCS Toolkit')
        self.master.iconbitmap('river.ico')

        self.bg_color = 'dark khaki'
        self.padding = 5

        ww = 700  # controls window width
        wh = 700  # controls window height
        wx = (self.master.winfo_screenwidth() - ww) / 2  # position relative to screen width and ww
        wy = (self.master.winfo_screenheight() - wh) / 2  # position relative to screen height and wh
        self.master.geometry("%dx%d+%d+%d" % (ww, wh, wx, wy))  # set window height and location

        # set widget styles
        self.style = ttk.Style()
        self.style.configure('TFrame', background=self.bg_color)
        self.style.configure('TButton', relief=RAISED)
        self.style.configure('TLabel', padding=self.padding, background=self.bg_color)
        self.style.configure('TEntry', borderwidth=self.padding, relief=SUNKEN)
        self.style.configure('TCheckbutton', background=self.bg_color)
        self.style.configure('TRadiobutton', background=self.bg_color)

        # initialize tab handler
        self.tab_container = ttk.Notebook(master)

        self.tab_names = ['Process LiDAR Data', 'Create Centerline', 'Create Station Lines', 'Detrend DEM',
                          'Extract GCS Data', 'Classify Landforms', 'GCS Analysis']

        self.tabs = {}
        for tab_name in self.tab_names:
            tab = ttk.Frame(self.tab_container)
            self.tab_container.add(tab, text=tab_name)
            self.tab_container.pack(expand=1, fill="both")
            self.tabs[tab_name] = tab

        # fill each tab with widgets
        ######################################################################

        root = self.tabs['Process LiDAR Data']

        self.l_lasbin = ttk.Label(root, text='LAStools /bin/ directory:')
        self.l_lasbin.grid(sticky=E, row=0, column=1)
        self.e_lasbin = ttk.Entry(root)
        self.e_lasbin.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
        self.e_lasbin.grid(row=0, column=2)
        self.b_lasbin = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_lasbin, select='folder'))
        self.b_lasbin.grid(sticky=W, row=0, column=3)

        self.l_lidardir = ttk.Label(root, text='LiDAR data directory:')
        self.l_lidardir.grid(sticky=E, row=1, column=1)
        self.e_lidardir = ttk.Entry(root)
        self.e_lidardir.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
        self.e_lidardir.grid(row=1, column=2)
        self.b_lidardir = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_lidardir, select='folder'))
        self.b_lidardir.grid(sticky=W, row=1, column=3)

        self.l_ground_shp = ttk.Label(root, text='Ground area .shp file (optional):')
        self.shp_var = StringVar()
        self.l_ground_shp.grid(sticky=E, row=2, column=1)
        self.e_ground_shp = ttk.Entry(root, textvariable=self.shp_var)
        self.e_ground_shp.grid(row=2, column=2)
        self.b_ground_shp = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_ground_shp,
                                                                                   select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                                          ('All files', '*')]
                                                                                    )
                        )
        self.b_ground_shp.grid(sticky=W, row=2, column=3)

        # if no ground shapefile is provided, disable the fine setting and just run on "coarse"
        def trace_choice(*args):
            fine_entries = [self.e_f_step, self.e_f_bulge, self.e_f_spike, self.e_f_dspike, self.e_f_offset]
            if self.shp_var.get() == '':
                for widget in fine_entries:
                    widget.config(state=DISABLED)
            else:
                for widget in fine_entries:
                    widget.config(state='normal')

        self.shp_var.trace('w', trace_choice)

        # specify lasground_new parameters

        root.grid_rowconfigure(5, minsize=80)

        self.l_coarse_class = ttk.Label(root, text='standard/coarse classification parameters:')
        self.l_coarse_class.grid(row=5, column=0, columnspan=2)

        self.l_c_step = ttk.Label(root, text='step size:')
        self.l_c_step.grid(sticky=E, row=6)
        self.e_c_step = ttk.Entry(root)
        self.e_c_step.grid(row=6, column=1)

        self.l_c_bulge = ttk.Label(root, text='bulge:')
        self.l_c_bulge.grid(sticky=E, row=7)
        self.e_c_bulge = ttk.Entry(root)
        self.e_c_bulge.grid(row=7, column=1)

        self.l_c_spike = ttk.Label(root, text='spike:')
        self.l_c_spike.grid(sticky=E, row=8)
        self.e_c_spike = ttk.Entry(root)
        self.e_c_spike.grid(row=8, column=1)

        self.l_c_dspike = ttk.Label(root, text='down spike:')
        self.l_c_dspike.grid(sticky=E, row=9)
        self.e_c_dspike = ttk.Entry(root)
        self.e_c_dspike.grid(row=9, column=1)

        self.l_c_offset = ttk.Label(root, text='offset:')
        self.l_c_offset.grid(sticky=E, row=10)
        self.e_c_offset = ttk.Entry(root)
        self.e_c_offset.grid(row=10, column=1)

        self.l_fine_class = ttk.Label(root, text='fine classification parameters (in ground area):')
        self.l_fine_class.grid(row=5, column=2, columnspan=2)

        self.l_f_step = ttk.Label(root, text='step size:')
        self.l_f_step.grid(sticky=E, row=6, column=2)
        self.e_f_step = ttk.Entry(root, state=DISABLED)
        self.e_f_step.grid(row=6, column=3)

        self.l_f_bulge = ttk.Label(root, text='bulge:')
        self.l_f_bulge.grid(sticky=E, row=7, column=2)
        self.e_f_bulge = ttk.Entry(root, state=DISABLED)
        self.e_f_bulge.grid(row=7, column=3)

        self.l_f_spike = ttk.Label(root, text='spike:')
        self.l_f_spike.grid(sticky=E, row=8, column=2)
        self.e_f_spike = ttk.Entry(root, state=DISABLED)
        self.e_f_spike.grid(row=8, column=3)

        self.l_f_dspike = ttk.Label(root, text='down spike:')
        self.l_f_dspike.grid(sticky=E, row=9, column=2)
        self.e_f_dspike = ttk.Entry(root, state=DISABLED)
        self.e_f_dspike.grid(row=9, column=3)

        self.l_f_offset = ttk.Label(root, text='offset:')
        self.l_f_offset.grid(sticky=E, row=10, column=2)
        self.e_f_offset = ttk.Entry(root, state=DISABLED)
        self.e_f_offset.grid(row=10, column=3)

        # specify units
        self.l_lidar_units = ttk.Label(root, text='Units')
        self.l_lidar_units.grid(sticky=W, row=11, column=2)
        root.grid_rowconfigure(11, minsize=80)
        self.lidar_units = StringVar()
        self.r_lidar_meters = ttk.Radiobutton(root, text='Meters', variable=self.lidar_units, value=' ')
        self.r_lidar_meters.grid(sticky=E, row=12, column=1)
        self.r_lidar_feet = ttk.Radiobutton(root, text='US Feet', variable=self.lidar_units, value=' -feet -elevation_feet ')
        self.r_lidar_feet.grid(row=12, column=2)
        self.lidar_units.set(' ')

        # specify number of cores
        self.l_lidar_cores = ttk.Label(root, text='Number of cores for processing')
        self.l_lidar_cores.grid(sticky=E, row=13, column=1, columnspan=2)
        root.grid_rowconfigure(13, minsize=80)
        self.core_num = IntVar()
        self.r1_lidar = ttk.Radiobutton(root, text='1', variable=self.core_num, value=1)
        self.r1_lidar.grid(sticky=E, row=14, column=1)
        self.r2_lidar = ttk.Radiobutton(root, text='2', variable=self.core_num, value=2)
        self.r2_lidar.grid(row=14, column=2)
        self.r4_lidar = ttk.Radiobutton(root, text='4', variable=self.core_num, value=4)
        self.r4_lidar.grid(sticky=W, row=14, column=3)
        self.r8_lidar = ttk.Radiobutton(root, text='8', variable=self.core_num, value=8)
        self.r8_lidar.grid(sticky=E, row=15, column=1)
        self.r16_lidar = ttk.Radiobutton(root, text='16', variable=self.core_num, value=16)
        self.r16_lidar.grid(row=15, column=2)
        self.r32_lidar = ttk.Radiobutton(root, text='32', variable=self.core_num, value=32)
        self.r32_lidar.grid(sticky=W, row=15, column=3)
        self.core_num.set(16)

        self.l_keep_orig_lidar = ttk.Label(root, text='Keep original ground/veg points: ')
        self.l_keep_orig_lidar.grid(sticky=E, row=16, column=1)
        self.keep_orig_lidar = BooleanVar()
        self.c_keep_orig_lidar = ttk.Checkbutton(root, variable=self.keep_orig_lidar)
        self.c_keep_orig_lidar.grid(sticky=W, row=16, column=2)
        self.keep_orig_lidar.set(True)

        # make 'Run' ttk.Button in GUI to call the process_lidar() function
        self.b_lidar_run = ttk.Button(root, text='    Run    ',
                                      command=lambda: lp.process_lidar(lastoolsdir=self.e_lasbin.get(),
                                                                       lidardir=self.e_lidardir.get(),
                                                                       ground_poly=self.e_ground_shp.get(),
                                                                       cores=self.core_num.get(),
                                                                       units_code=self.lidar_units.get()[1:-1],
                                                                       keep_orig_pts=self.keep_orig_lidar.get(),
                                                                       coarse_step=self.e_c_step.get(),
                                                                       coarse_bulge=self.e_c_bulge.get(),
                                                                       coarse_spike=self.e_c_spike.get(),
                                                                       coarse_down_spike=self.e_c_dspike.get(),
                                                                       coarse_offset=self.e_c_offset.get(),
                                                                       fine_step=self.e_f_step.get(),
                                                                       fine_bulge=self.e_f_bulge.get(),
                                                                       fine_spike=self.e_f_spike.get(),
                                                                       fine_down_spike=self.e_f_dspike.get(),
                                                                       fine_offset=self.e_f_offset.get()
                                                                       )
                                      )

        self.b_lidar_run.grid(sticky=W, row=17, column=2)
        root.grid_rowconfigure(17, minsize=80)
        
        #########################################################################
        
        root = self.tabs['Create Centerline']

        self.l_cc_dem = ttk.Label(root, text='DEM:')
        self.l_cc_dem.grid(sticky=E, row=0, column=1)
        self.e_cc_dem = ttk.Entry(root)
        self.e_cc_dem.grid(row=0, column=2)
        self.b_cc_dem = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_cc_dem, select='file',
                                                                               ftypes=[('Raster', '*.tif'),
                                                                                       ('Raster .aux', '*.aux, *.aux.xml'),
                                                                                       ('All files', '*')]
                                                                               )
                                   )
        self.b_cc_dem.grid(sticky=W, row=0, column=3)

        self.l_cc_channel = ttk.Label(root, text='Channel Polygon:')
        self.l_cc_channel.grid(sticky=E, row=1, column=1)
        self.e_cc_channel = ttk.Entry(root)
        self.e_cc_channel.grid(row=1, column=2)
        self.b_cc_channel = ttk.Button(root, text='Browse',
                                       command=lambda: browse(root, self.e_cc_channel, select='file',
                                                              ftypes=[('Shapefile', '*.shp'),
                                                                      ('All files', '*')]
                                                              )
                                       )
        self.b_cc_channel.grid(sticky=W, row=1, column=3)

        self.l_flow_src = ttk.Label(root, text='Flow Source Polygon:')
        self.l_flow_src.grid(sticky=E, row=2, column=1)
        self.e_flow_src = ttk.Entry(root)
        self.e_flow_src.grid(row=2, column=2)
        self.b_flow_src = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_flow_src, select='file',
                                                                                 ftypes=[('Shapefile', '*.shp'),
                                                                                         ('All files', '*')]
                                                                                 )
                                     )
        self.b_flow_src.grid(sticky=W, row=2, column=3)

        self.l_smooth_dist = ttk.Label(root, text='Smoothing Distance:')
        self.l_smooth_dist.grid(sticky=E, row=8, column=1)
        self.e_smooth_dist = ttk.Entry(root)
        self.e_smooth_dist.insert(END, 20)
        self.e_smooth_dist.grid(row=8, column=2)

        self.b_cc_run = ttk.Button(root, text='   Run    ',
                                   command=lambda: cc.make_centerline(DEM=self.e_cc_dem.get(),
                                                                      channel=self.e_cc_channel.get(),
                                                                      source=self.e_flow_src.get(),
                                                                      smooth_distance=float(self.e_smooth_dist.get())
                                                                      )
                                   )
        self.b_cc_run.grid(sticky=W, row=9, column=2)
        root.grid_rowconfigure(9, minsize=80)

        #########################################################################

        root = self.tabs['Create Station Lines']

        self.l_sl_centerline = ttk.Label(root, text='Centerline:')
        self.l_sl_centerline.grid(sticky=E, row=0, column=1)
        self.e_sl_centerline = ttk.Entry(root)
        self.e_sl_centerline.grid(row=0, column=2)
        self.b_sl_centerline = ttk.Button(root, text='Browse',
                                          command=lambda: browse(root, self.e_sl_centerline, select='file',
                                                                 ftypes=[('Shapefile', '*.shp'),
                                                                         ('All files', '*')]
                                                                 )
                                          )
        self.b_sl_centerline.grid(sticky=W, row=0, column=3)

        self.l_xs_spacing = ttk.Label(root, text='Spacing:')
        self.l_xs_spacing.grid(sticky=E, row=1, column=1)
        self.e_xs_spacing = ttk.Entry(root)
        self.e_xs_spacing.insert(END, 3)
        self.e_xs_spacing.grid(row=1, column=2)

        self.l_xs_length = ttk.Label(root, text='XS Length:')
        self.l_xs_length.grid(sticky=E, row=2, column=1)
        self.e_xs_length = ttk.Entry(root)
        self.e_xs_length.insert(END, 100)
        self.e_xs_length.grid(row=2, column=2)

        self.b_sl_run = ttk.Button(root, text='   Run    ',
                                   command=lambda: sl.create_station_lines(line_shp=self.e_sl_centerline.get(),
                                                                           spacing=float(self.e_xs_spacing.get()),
                                                                           xs_length=float(self.e_xs_length.get()))
                                   )
        self.b_sl_run.grid(sticky=W, row=3, column=2)
        root.grid_rowconfigure(3, minsize=80)

        #########################################################################

        root = self.tabs['Detrend DEM']

        self.l_dd_dem = ttk.Label(root, text='DEM: ')
        self.l_dd_dem.grid(sticky=E, row=0, column=0)
        self.e_dd_dem = ttk.Entry(root)
        self.e_dd_dem.grid(row=0, column=1)
        self.b_dd_dem = ttk.Button(root, text='Browse',
                                   command=lambda: browse(root, self.e_dd_dem, select='file',
                                                          ftypes=[('Raster', '*.tif'),
                                                                  ('Raster .aux', '*.aux, *.aux.xml'),
                                                                  ('All files', '*')]
                                                          )
                                   )
        self.b_dd_dem.grid(sticky=W, row=0, column=2)

        self.l_dd_centerline = ttk.Label(root, text='Centerline: ')
        self.l_dd_centerline.grid(sticky=E, row=1, column=0)
        self.e_dd_centerline = ttk.Entry(root)
        self.e_dd_centerline.grid(row=1, column=1)
        self.b_dd_centerline = ttk.Button(root, text='Browse',
                                          command=lambda: browse(root, self.e_dd_centerline, select='file',
                                                                 ftypes=[('Shapefile', '*.shp'),
                                                                         ('All files', '*')]
                                                                 )
                                          )
        self.b_dd_centerline.grid(sticky=W, row=1, column=2)

        self.l_dd_station_lines = ttk.Label(root, text='Station XS Lines: ')
        self.l_dd_station_lines.grid(sticky=E, row=2, column=0)
        self.e_dd_station_lines = ttk.Entry(root)
        self.e_dd_station_lines.grid(row=2, column=1)
        self.b_dd_station_lines = ttk.Button(root, text='Browse',
                                             command=lambda: browse(root, self.e_dd_station_lines, select='file',
                                                                    ftypes=[('Shapefile', '*.shp'),
                                                                            ('All files', '*')]
                                                                    )
                                             )
        self.b_dd_station_lines.grid(sticky=W, row=2, column=2)

        self.l_dd_slope_breaks = ttk.Label(root, text='Slope Breaks (distance downstream): ')
        self.l_dd_slope_breaks.grid(sticky=E, row=3, column=0)
        self.e_dd_slope_breaks = ttk.Entry(root)
        self.e_dd_slope_breaks.grid(row=3, column=1)

        self.l_dd_regression = ttk.Label(root, text='Regression: ')
        self.l_dd_regression.grid(sticky=E, row=4, column=0)
        self.cb_dd_regression = ttk.Combobox(root, values=('linear', 'quadratic'), state='readonly')
        self.cb_dd_regression.current(0)
        self.cb_dd_regression.grid(row=4, column=1)

        self.b_dd_main = ttk.Button(root, text='    Run    ',
                                    command=lambda: dd.main_det(self.e_dd_dem.get(),
                                                                self.e_dd_centerline.get(),
                                                                self.e_dd_station_lines.get(),
                                                                map(float, self.e_dd_slope_breaks.get().split(',')) if self.e_dd_slope_breaks.get() != '' else '',
                                                                regression=self.cb_dd_regression.get()
                                                                )
                                    )
        self.b_dd_main.grid(sticky=W, row=5, column=1, columnspan=2)

        #########################################################################

        root = self.tabs['Extract GCS Data']

        self.l_ec_station_lines = ttk.Label(root, text='Station Lines:')
        self.l_ec_station_lines.grid(sticky=E, row=0, column=1)
        self.e_ec_station_lines = ttk.Entry(root)
        self.e_ec_station_lines.grid(row=0, column=2)
        self.b_ec_station_lines = ttk.Button(root, text='Browse',
                                             command=lambda: browse(root, self.e_ec_station_lines, select='file',
                                                                    ftypes=[('Shapefile', '*.shp'),
                                                                            ('All files', '*')]
                                                                    )
                                             )
        self.b_ec_station_lines.grid(sticky=W, row=0, column=3)

        self.l_ec_det_dem = ttk.Label(root, text='Detrended DEM:')
        self.l_ec_det_dem.grid(sticky=E, row=1, column=1)
        self.e_ec_det_dem = ttk.Entry(root)
        self.e_ec_det_dem.grid(row=1, column=2)
        self.b_ec_det_dem = ttk.Button(root, text='Browse',
                                       command=lambda: browse(root, self.e_ec_det_dem, select='file',
                                                              ftypes=[('Raster', '*.tif'),
                                                                      ('All files', '*')]
                                                              )
                                       )
        self.b_ec_det_dem.grid(sticky=W, row=1, column=3)

        self.l_ec_wetted_polys = ttk.Label(root, text='Wetted Polygons/Velocity Rasters:')
        self.l_ec_wetted_polys.grid(sticky=E, row=2, column=1)
        self.e_ec_wetted_polys = ttk.Entry(root)
        self.e_ec_wetted_polys.grid(row=2, column=2)
        self.b_ec_wetted_polys = ttk.Button(root, text='Browse',
                                            command=lambda: browse(root, self.e_ec_wetted_polys, select='files',
                                                                   ftypes=[('Shapefile', '*.shp'),
                                                                           ('Float Raster', '*.flt'),
                                                                           ('All files', '*')]
                                                                   )
                                            )
        self.b_ec_wetted_polys.grid(sticky=W, row=2, column=3)

        self.l_xs_buffer = ttk.Label(root, text='XS Buffer Size:')
        self.l_xs_buffer.grid(sticky=E, row=8, column=1)
        self.e_xs_buffer = ttk.Entry(root)
        self.e_xs_buffer.grid(row=8, column=2)

        self.l_rm_up_length = ttk.Label(root, text='Remove Upstream Length:')
        self.l_rm_up_length.grid(sticky=E, row=9, column=1)
        self.e_rm_up_length = ttk.Entry(root)
        self.e_rm_up_length.insert(END, '0')
        self.e_rm_up_length.grid(row=9, column=2)

        self.l_rm_down_length = ttk.Label(root, text='Remove Downstream Length:')
        self.l_rm_down_length.grid(sticky=E, row=10, column=1)
        self.e_rm_down_length = ttk.Entry(root)
        self.e_rm_down_length.insert(END, '0')
        self.e_rm_down_length.grid(row=10, column=2)

        self.l_ec_reach_breaks = ttk.Label(root, text='Reach Breaks (distance downstream):')
        self.l_ec_reach_breaks.grid(sticky=E, row=11, column=1)
        self.e_ec_reach_breaks = ttk.Entry(root)
        self.e_ec_reach_breaks.grid(row=11, column=2)

        self.b_ec_run = ttk.Button(root, text='   Run    ',
                                   command=lambda: ec.extract_channel_data(station_lines=self.e_ec_station_lines.get(),
                                                                           detrended_DEM=self.e_ec_det_dem.get(),
                                                                           wetted_rasters_list=list(root.tk.splitlist(self.e_ec_wetted_polys.get())),
                                                                           buffer_size=self.e_xs_buffer.get(),
                                                                           rm_up_length=float(self.e_rm_up_length.get()),
                                                                           rm_down_length=float(self.e_rm_down_length.get()),
                                                                           reach_breaks=map(int, self.e_ec_reach_breaks.get().split(',')) if self.e_ec_reach_breaks.get() != '' else ''
                                                                           )
                                   )
        self.b_ec_run.grid(sticky=W, row=12, column=2)
        root.grid_rowconfigure(12, minsize=80)

        #########################################################################

        root = self.tabs['Classify Landforms']

        self.l_cl_xs_attributes = ttk.Label(root, text='Wetted XS Attribute Tables:')
        self.l_cl_xs_attributes.grid(sticky=E, row=0, column=1)
        self.e_cl_xs_attributes = ttk.Entry(root)
        self.e_cl_xs_attributes.grid(row=0, column=2)
        self.b_cl_xs_attributes = ttk.Button(root, text='Browse',
                                             command=lambda: browse(root, self.e_cl_xs_attributes, select='files',
                                                                    ftypes=[('Comma-delimited text', '.csv'),
                                                                            ('All files', '*')]
                                                                    )
                                             )
        self.b_cl_xs_attributes.grid(sticky=W, row=0, column=3)

        self.l_var1_field = ttk.Label(root, text='Width Field:')
        self.l_var1_field.grid(sticky=E, row=1, column=1)
        self.e_var1_field = ttk.Entry(root)
        self.e_var1_field.insert(END, 'W')
        self.e_var1_field.grid(row=1, column=2)

        self.l_var2_field = ttk.Label(root, text='Detrended Elevation Field:')
        self.l_var2_field.grid(sticky=E, row=2, column=1)
        self.e_var2_field = ttk.Entry(root)
        self.e_var2_field.insert(END, 'Z')
        self.e_var2_field.grid(row=2, column=2)

        self.l_var3_field = ttk.Label(root, text='Velocity Field:')
        self.l_var3_field.grid(sticky=E, row=3, column=1)
        self.e_var3_field = ttk.Entry(root)
        self.e_var3_field.insert(END, 'V')
        self.e_var3_field.grid(row=3, column=2)

        self.l_dist_down_field = ttk.Label(root, text='Distance Downstream Field:')
        self.l_dist_down_field.grid(sticky=E, row=4, column=1)
        self.e_dist_down_field = ttk.Entry(root)
        self.e_dist_down_field.insert(END, 'dist_down')
        self.e_dist_down_field.grid(row=4, column=2)

        self.b_cl_run = ttk.Button(root, text='   Run    ',
                                   command=lambda: cl.main_classify_landforms(tables=list(root.tk.splitlist(self.e_cl_xs_attributes.get())),
                                                                              w_field=self.e_var1_field.get(),
                                                                              z_field=self.e_var2_field.get(),
                                                                              v_field=self.e_var3_field.get(),
                                                                              dist_field=self.e_dist_down_field.get()
                                                                              )
                                   )
        self.b_cl_run.grid(sticky=W, row=9, column=2)
        root.grid_rowconfigure(9, minsize=80)

        #########################################################################

        root = self.tabs['GCS Analysis']

        self.l_gcs_tables = ttk.Label(root, text='GCS Data Tables: ')
        self.l_gcs_tables.grid(sticky=E, row=0, column=1)
        self.e_gcs_tables = ttk.Entry(root)
        self.e_gcs_tables.grid(row=0, column=2)
        self.b_gcs_tables = ttk.Button(root, text='Browse',
                                       command=lambda: browse(root, self.e_gcs_tables, select='files',
                                                              ftypes=[('Comma-delimited text', '.csv'),
                                                                      ('All files', '*')]
                                                              )
                                       )
        self.b_gcs_tables.grid(sticky=W, row=0, column=3)

        self.b_gcsa_run = ttk.Button(root, text='    Run    ',
                                     command=lambda: gcsa.complete_analysis(tables=list(root.tk.splitlist(self.e_gcs_tables.get())))
                                     )
        self.b_gcsa_run.grid(sticky=W, row=2, column=2)


if __name__ == '__main__':

    # initialize logger
    init_logger(__file__)

    GCS_GUI().mainloop()
