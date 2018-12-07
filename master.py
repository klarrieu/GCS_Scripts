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


class GCS_GUI(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()

        # set header
        self.master.title('River GCS App')
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

        print(self.tabs)

        # fill each tab with widgets
        ######################################################################

        root = self.tabs['Process LiDAR Data']

        L1 = ttk.Label(root, text='LAStools /bin/ directory:')
        L1.grid(sticky=E, row=0, column=1)
        E1 = ttk.Entry(root)
        E1.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
        E1.grid(row=0, column=2)
        b1 = ttk.Button(root, text='Browse', command=lambda: browse(root, E1, select='folder'))
        b1.grid(sticky=W, row=0, column=3)

        L2 = ttk.Label(root, text='LiDAR data directory:')
        L2.grid(sticky=E, row=1, column=1)
        E2 = ttk.Entry(root)
        E2.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
        E2.grid(row=1, column=2)
        b2 = ttk.Button(root, text='Browse', command=lambda: browse(root, E2, select='folder'))
        b2.grid(sticky=W, row=1, column=3)

        L3 = ttk.Label(root, text='Ground area .shp file (optional):')
        shp_var = StringVar()
        L3.grid(sticky=E, row=2, column=1)
        E3 = ttk.Entry(root, textvariable=shp_var)
        E3.grid(row=2, column=2)
        b3 = ttk.Button(root, text='Browse', command=lambda: browse(root, E3, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                                 ('All files', '*')]
                                                                )
                        )
        b3.grid(sticky=W, row=2, column=3)

        # if no ground shapefile is provided, disable the fine setting and just run on "coarse"
        def trace_choice(*args):
            if shp_var.get() == '':
                for widget in [E1b, E2b, E3b, E4b, E5b]:
                    widget.config(state=DISABLED)
            else:
                for widget in [E1b, E2b, E3b, E4b, E5b]:
                    widget.config(state='normal')

        shp_var.trace('w', trace_choice)

        # specify lasground_new parameters

        root.grid_rowconfigure(5, minsize=80)

        LC1 = ttk.Label(root, text='standard/coarse classification parameters:')
        LC1.grid(row=5, column=0, columnspan=2)

        L1a = ttk.Label(root, text='step size:')
        L1a.grid(sticky=E, row=6)
        E1a = ttk.Entry(root)
        E1a.grid(row=6, column=1)

        L2a = ttk.Label(root, text='bulge:')
        L2a.grid(sticky=E, row=7)
        E2a = ttk.Entry(root)
        E2a.grid(row=7, column=1)

        L3a = ttk.Label(root, text='spike:')
        L3a.grid(sticky=E, row=8)
        E3a = ttk.Entry(root)
        E3a.grid(row=8, column=1)

        L4a = ttk.Label(root, text='down spike:')
        L4a.grid(sticky=E, row=9)
        E4a = ttk.Entry(root)
        E4a.grid(row=9, column=1)

        L5a = ttk.Label(root, text='offset:')
        L5a.grid(sticky=E, row=10)
        E5a = ttk.Entry(root)
        E5a.grid(row=10, column=1)

        LC2 = ttk.Label(root, text='fine classification parameters (in ground area):')
        LC2.grid(row=5, column=2, columnspan=2)

        L1b = ttk.Label(root, text='step size:')
        L1b.grid(sticky=E, row=6, column=2)
        E1b = ttk.Entry(root, state=DISABLED)
        E1b.grid(row=6, column=3)

        L2b = ttk.Label(root, text='bulge:')
        L2b.grid(sticky=E, row=7, column=2)
        E2b = ttk.Entry(root, state=DISABLED)
        E2b.grid(row=7, column=3)

        L3b = ttk.Label(root, text='spike:')
        L3b.grid(sticky=E, row=8, column=2)
        E3b = ttk.Entry(root, state=DISABLED)
        E3b.grid(row=8, column=3)

        L4b = ttk.Label(root, text='down spike:')
        L4b.grid(sticky=E, row=9, column=2)
        E4b = ttk.Entry(root, state=DISABLED)
        E4b.grid(row=9, column=3)

        L5b = ttk.Label(root, text='offset:')
        L5b.grid(sticky=E, row=10, column=2)
        E5b = ttk.Entry(root, state=DISABLED)
        E5b.grid(row=10, column=3)

        # specify units
        L5 = ttk.Label(root, text='Units')
        L5.grid(sticky=W, row=11, column=2)
        root.grid_rowconfigure(11, minsize=80)
        unit_var = StringVar()
        R5m = ttk.Radiobutton(root, text='Meters', variable=unit_var, value=' ')
        R5m.grid(sticky=E, row=12, column=1)
        R5f = ttk.Radiobutton(root, text='US Feet', variable=unit_var, value=' -feet -elevation_feet ')
        R5f.grid(row=12, column=2)
        unit_var.set(' ')

        # specify number of cores
        L4 = ttk.Label(root, text='Number of cores for processing')
        L4.grid(sticky=E, row=13, column=1, columnspan=2)
        root.grid_rowconfigure(13, minsize=80)
        core_num = IntVar()
        R1 = ttk.Radiobutton(root, text='1', variable=core_num, value=1)
        R1.grid(sticky=E, row=14, column=1)
        R2 = ttk.Radiobutton(root, text='2', variable=core_num, value=2)
        R2.grid(row=14, column=2)
        R4 = ttk.Radiobutton(root, text='4', variable=core_num, value=4)
        R4.grid(sticky=W, row=14, column=3)
        R8 = ttk.Radiobutton(root, text='8', variable=core_num, value=8)
        R8.grid(sticky=E, row=15, column=1)
        R16 = ttk.Radiobutton(root, text='16', variable=core_num, value=16)
        R16.grid(row=15, column=2)
        R32 = ttk.Radiobutton(root, text='32', variable=core_num, value=32)
        R32.grid(sticky=W, row=15, column=3)
        core_num.set(16)

        L5 = ttk.Label(root, text='Keep original ground/veg points: ')
        L5.grid(sticky=E, row=16, column=1)
        keep_originals = BooleanVar()
        C1 = ttk.Checkbutton(root, variable=keep_originals)
        C1.grid(sticky=W, row=16, column=2)
        keep_originals.set(True)

        # make 'Run' ttk.Button in GUI to call the process_lidar() function
        b = ttk.Button(root, text='    Run    ', command=lambda: lp.process_lidar(lastoolsdir=E1.get(),
                                                                                  lidardir=E2.get(),
                                                                                  ground_poly=E3.get(),
                                                                                  cores=core_num.get(),
                                                                                  units_code=unit_var.get()[1:-1],
                                                                                  keep_orig_pts=keep_originals.get(),
                                                                                  coarse_step=E1a.get(),
                                                                                  coarse_bulge=E2a.get(),
                                                                                  coarse_spike=E3a.get(),
                                                                                  coarse_down_spike=E4a.get(),
                                                                                  coarse_offset=E5a.get(),
                                                                                  fine_step=E1b.get(),
                                                                                  fine_bulge=E2b.get(),
                                                                                  fine_spike=E3b.get(),
                                                                                  fine_down_spike=E4b.get(),
                                                                                  fine_offset=E5b.get()
                                                                                  )
                       )

        b.grid(sticky=W, row=17, column=2)
        root.grid_rowconfigure(17, minsize=80)
        
        #########################################################################
        
        root = self.tabs['Create Centerline']

        L1 = ttk.Label(root, text='DEM:')
        L1.grid(sticky=E, row=0, column=1)
        E1 = ttk.Entry(root)
        E1.grid(row=0, column=2)
        b1 = ttk.Button(root, text='Browse', command=lambda: browse(root, E1, select='file', ftypes=[('Raster', '*.tif'),
                                                                                                 ('All files', '*')]
                                                                )
                    )
        b1.grid(sticky=W, row=0, column=3)

        L2 = ttk.Label(root, text='Channel Polygon:')
        L2.grid(sticky=E, row=1, column=1)
        E2 = ttk.Entry(root)
        E2.grid(row=1, column=2)
        b2 = ttk.Button(root, text='Browse', command=lambda: browse(root, E2, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                                 ('All files', '*')]
                                                                )
                    )
        b2.grid(sticky=W, row=1, column=3)

        L3 = ttk.Label(root, text='Flow Source Polygon:')
        L3.grid(sticky=E, row=2, column=1)
        E3 = ttk.Entry(root)
        E3.grid(row=2, column=2)
        b3 = ttk.Button(root, text='Browse', command=lambda: browse(root, E3, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                                 ('All files', '*')]
                                                                )
                    )
        b3.grid(sticky=W, row=2, column=3)

        L4 = ttk.Label(root, text='Smoothing Distance:')
        L4.grid(sticky=E, row=8, column=1)
        E4 = ttk.Entry(root)
        E4.insert(END, 20)
        E4.grid(row=8, column=2)

        b = ttk.Button(root, text='   Run    ', command=lambda: make_centerline(DEM=E1.get(),
                                                                            channel=E2.get(),
                                                                            source=E3.get(),
                                                                            smooth_distance=float(E4.get())
                                                                            )
                   )
        b.grid(sticky=W, row=9, column=2)
        root.grid_rowconfigure(9, minsize=80)

        #########################################################################

if __name__ == '__main__':

    GCS_GUI().mainloop()
