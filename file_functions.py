import os, time
import arcpy
from arcpy import env
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')
import pandas as pd
from Tkinter import *
import tkFileDialog
import subprocess
import logging
logger = logging.getLogger(__name__)

def init_logger(filename):
    '''Initializes logger'''
    logging.basicConfig(filename = os.path.basename(filename).replace('.py','.log'), filemode = 'w', level=logging.INFO)
    stderrLogger=logging.StreamHandler()
    stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logging.getLogger().addHandler(stderrLogger)
    return

def cmd(command):
    '''Executes command prompt command'''
    try:
        res = subprocess.Popen(command, stderr = subprocess.PIPE)
    except:
        msg = 'Command failed: %s'%command
        logger.error(msg)
        raise Exception(msg)
        return
        
    if res.wait() != 0:
        msg = res.communicate()[1]
        logger.error(msg)
        raise Exception(msg)
    
    return

#opens window in GUI to browse for folder or file
def browse(root, entry, select = 'file', ftypes = [('All files', '*')]):
    '''GUI button command: opens browser window and adds selected file/folder to entry'''
    if select == 'file':
        filename = tkFileDialog.askopenfilename(parent = root,title = 'Choose a file', filetypes = ftypes)
        if filename != None:
            entry.delete(0, END)
            entry.insert(END, filename)

    elif select == 'files':
        files = tkFileDialog.askopenfilenames(parent = root, title = 'Choose files', filetypes = ftypes)
        l = root.tk.splitlist(files)
        entry.delete(0, END)
        entry.insert(END, l)
    
    elif select == 'folder':
        dirname = tkFileDialog.askdirectory(parent = root,initialdir = entry.get(), title = 'Choose a directory')
        if len(dirname ) > 0:
            entry.delete(0, END)
            entry.insert(END, dirname+'/')


def check_use(filepath):
    """Checks if a file or list of files is in use by another process
    If the file cannot be opened or there is an associated .lock file, it throws an exception.
    """
    
    if type(filepath) == list:
        for f in filepath:
            check_use(f)
        return
    
    file_object = None
    if os.path.exists(filepath):
        try:
            buffer_size = 8
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'a', buffer_size)
            if file_object:
                for filename in os.listdir(os.path.dirname(filepath)):
                    if filename.startswith(os.path.basename(filepath)) and filename.endswith('.lock'):
                        logger.error('%s is open in another program. Close the file and try again.'%filepath)
                        raise Exception('%s is open in another program. Close the file and try again.'%filepath)
                
        except IOError, message:
            logger.error('%s is open in another program. Close the file and try again.'%filepath)
            raise Exception('%s is open in another program. Close the file and try again.'%filepath)
        
        finally:
            if file_object:
                file_object.close()
    return

def split_list(l, break_pts):
    '''returns list l split up into sublists at break point indices'''
    l_0 = len(l)
    sl = []
    if break_pts == []:
        return [l]
    else:
        for brk in break_pts:
            delta_l = l_0-len(l)
            sl.append(l[:brk-delta_l])
            l = l[brk-delta_l:]
        sl.append(l)
    return sl

def split_reaches(l, new_reach_pts):
    '''splits l into sections where new_reach_pts contains the starting indices for each slice'''
    new_reach_pts = sorted(new_reach_pts)
    sl = [l[i1:i2] for i1,i2 in zip(new_reach_pts,new_reach_pts[1:])]
    last_index = new_reach_pts[-1]
    sl.append(l[last_index:])
    return sl

class DF(pd.DataFrame):
    '''pandas DataFrame class with an additional title attribute'''
    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False, title=None):
        pd.DataFrame.__init__(self, data, index, columns, dtype, copy)
        self.title = title

    def show(self):
        if self.title != None:
            print(self.title)
        print(self)

def flt_to_poly(flt):
    '''Converts .flt raster to a single polygon covering area that is not null'''
    ras = arcpy.Raster(flt)
    #make integer raster
    int_raster = arcpy.sa.Con(arcpy.sa.IsNull(ras) == False, 1)
    #convert to polygon
    poly = arcpy.RasterToPolygon_conversion(int_raster,
                                            flt.replace('.flt','.shp'),
                                            'NO_SIMPLIFY'
                                            )

    return poly.getOutput(0)


# Test
if __name__ == '__main__':
    files = [r"G:\Kenny_Rainbow_Basin\DEM\detrending\station_coords.xls",
             r"G:\Kenny_Rainbow_Basin\DEM\detrending\xyz_fit.csv"]

