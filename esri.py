"""
Class for reading and manipulating ESRI ASCII-format rasters.
"""

import os
import numpy as np
from weakref import WeakKeyDictionary

class NonNegative(object):

    """
    Decriptor for non-negative integer attributes.
    """
    
    def __init__(self, default):
        self.default = default
        self.data = WeakKeyDictionary()
    
    def __get__(self, instance, owner):
        return self.data.get(instance, self.default)
        
    def __set__(self, instance, value):
        if not isinstance(value, int):
            raise ValueError('must be an integer value: {}'.format(value))
        if value < 0:
            raise ValueError('negative values cannot be assigned: {}'.format(value))
        self.data[instance] = value
        

class AsciiGrid(object):

    """
    ESRI ASCII-format raster class.
    
    Example:
        dem = AsciiGrid('dem.asc')
        print dem.ncols, dem.nrows
        dem.write('dem_plus_1.asc')
    
    Attributes:
        ncols (int):        Number of columns in AsciiGrid (read-only).
        nrows (int):        Number of rows in AsciiGrid (read-only).
        xllcorner (int):    X-coordinate of lower-left grid corner.
        yllcorner (int):    Y-coordinate of lower-left grid corner.
        cellsize (int):     Grid cell size (in metres).
        no_data_val (int):  No-data value (default = -9999).
        data (np.ndarray):  Raster array with no-data values represented as numpy.nan's.
    
    TODO:
        * See https://github.com/olemb/dbfread/blob/master/dbfread/dbf.py for example!
        * Options in write() to specify xllcorner vs. xllcenter
        * Implement __repr__() and __str__() methods.
        * Error-checking on filepaths for read, write methods.
        * Implement some way to compare grids (i.e. __eq__, same_extents(), etc...).
        * Implement a 'change_no_data_val' method:
          --> Can probably do this with a @property or something similar? So just 'step in' when
              the no_data_val attribute is changed and change it in the data data too :).
        * Replace add_grid() method with __add__().
        * Implement a way to access data via own coordinate system???
    """
    
    # version number
    version = '0.1'

    # descriptor attributes
    xllcorner = NonNegative(0)
    yllcorner = NonNegative(0)
    cellsize  = NonNegative(1)
    
    # property attributes
    @property
    def ncols(self):
        return self.data.shape[1]
    
    @property
    def nrows(self):
        return self.data.shape[0]
        
    @property
    def no_data_val(self):
        return self._no_data_val
        
    @no_data_val.setter
    def no_data_val(self, value):
        if not isinstance(value, int):
            raise ValueError('no_data_val must take an integer value: {}'.format(value))
        self._no_data_val = value
        
    @property
    def data(self):
        return self._data
        
    @data.setter
    def data(self, arr):
        if not isinstance(arr, np.ndarray):
            raise TypeError('data attribute must be a numpy.ndarray object')
        self._data = np.copy(arr)
    
    
    def __init__(self, xllcorner=0, yllcorner=0, cellsize=1, 
                 no_data_val=-9999, data=np.zeros((10,10))):
        """
        Initialiser: Create an AsciiGrid raster object from kwargs and a numpy array. 
        Defaults to a 10x10 grid of zeros.
        """
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.no_data_val = no_data_val
        self.data = np.where(data == no_data_val, np.nan, data)
        
    
    # TODO: def __repr__(self):
    # TODO: def __str(self)__:
    
        # """Pretty-print ascii-grid object."""

        # print 'xllcorner: {s.xllcorner:d}'.format(s=self)
        # print 'yllcorner: {s.yllcorner:d}'.format(s=self)
    
    
    def read(self, filename):
    
        """Reads an ascii grid from the filename (deprecated; maintained for backwards compatibility)."""
        
        return AsciiGrid.fromfile(filename)
        
        # try:
            # fh = open(filename, 'rb')
        # except IOError:
            # print 'can\'t find or read ASCII raster file.'
        # else:
            # # read header data
            # self.ncols = int(fh.readline().split()[1])
            # self.nrows = int(fh.readline().split()[1])
            # self.xllcorner = int(float(fh.readline().split()[1]))
            # self.yllcorner = int(float(fh.readline().split()[1]))
            # self.cellsize = int(float(fh.readline().split()[1]))
            # self.no_data_val = float(fh.readline().split()[1])
            
            # # read raster data
            # arr = np.fromfile(fh, dtype=float, count=-1, sep=' ')
            # self.data = np.reshape(arr, (self.nrows, self.ncols))
            
            # # close handle
            # fh.close()
            
            # return self
            
    
    @classmethod
    def fromfile(cls, filename):
    
        """Reads an ESRI ASCII-format raster from file."""
        
        # read file headers, data array
        try:
            fh = open(filename, 'r')
        except IOError:
            print 'File not found/error reading file'
        else:
            hdrs = [fh.readline().split() for i in xrange(6)]
            arr = np.fromfile(fh, dtype=float, count=-1, sep=' ')
            fh.close()
            
        # parse headers
        d_hd = dict([(tup[0].lower(),int(tup[1])) for tup in hdrs])
        try:
            ncols = d_hd.pop('ncols')
            nrows = d_hd.pop('nrows')
            cellsize = d_hd.pop('cellsize')
            no_data_val = d_hd.pop('nodata_value',-9999)
            if 'xllcorner' in d_hd:
                xllcorner = d_hd.pop('xllcorner')
            else:
                xllcorner = int(d_hd.pop('xllcenter') - cellsize*0.5)
            if 'yllcorner' in d_hd:
                yllcorner = d_hd.pop('yllcorner')
            else:
                yllcorner = int(d_hd.pop('yllcenter') - cellsize*0.5)
        except KeyError:
            print 'Missing header entries in file'
            raise
        if d_hd:
            raise TypeError('Unexpected header entries found in file: {}'.format(d_hd.keys()))
        
        # check data, create instance
        try:
            data = np.reshape(arr, (nrows,ncols))
        except ValueError:
            print 'NROWS, NCOLS inconsistent with raster array size'
            raise
        else:
            return cls(xllcorner=xllcorner, yllcorner=yllcorner, cellsize=cellsize, 
                       no_data_val=no_data_val, data=data)
    
    
    def write(self, filename, precision=2):
    
        """Writes the raster to ascii-format file specified in first argument (filename)."""
        
        # # check filename
        # if not os.path.isdir(filename):
            # raise IOError('filepath not recognised: {}'.format(filename))
        
        # check precision
        if not isinstance(precision, int):
            raise ValueError('precision argument takes an integer value: {}'.format(precision))
        fmt = '%1.'+str(int(precision))+'f'
        
        # write header
        header = ''.join(['ncols         ',str(self.ncols),'\n',
                          'nrows         ',str(self.nrows),'\n',
                          'xllcorner     ',str(self.xllcorner),'\n',
                          'yllcorner     ',str(self.yllcorner),'\n',
                          'cellsize      ',str(self.cellsize),'\n',
                          'NODATA_value  ',str(self.no_data_val)])
        
        # write grid data
        np.savetxt(filename, self.data, delimiter=' ', header=header, fmt=fmt, comments='')
        
    
    def get_headers(self):
    
        """Returns headers as a dict."""
        
        h = {k:v for k,v in vars(self).iteritems() if k != 'data'}
        return h
        
        
    def add_grid(self, newAsciiGrid):
    # def __add__(self, newAsciiGrid):
    
        """Adds an AsciiGrid raster."""
        
        # TODO: check same xllcorner, yllcorner, shape, etc...
        
        # error checking
        # if newAsciiGrid.shape != self.data.shape:
            # print 'Error: grid must have same shape as Ascii raster grid.'
            # raise Exception
    
        try:
            self.data += newAsciiGrid.data
        except:
            print 'error: cannot add grids'
        
        
    def clip(self, xmin=None, xmax=None, ymin=None, ymax=None):
    
        """Clips the raster to specified limits (xmin, xmax, ymin & ymax)."""
        
        # sanity-check clipping limits:
        #  1. must be numeric
        try:
            int(xmin); int(xmax); int(ymin); int(ymax)
        except ValueError:
            print 'Error: non-numerical value specified for clipping limits'
        
        #  2. must align with grid (i.e. be a multiple of cellsize removed from llcorner)
        for arg in [xmin-self.xllcorner, xmax-self.xllcorner, ymin-self.yllcorner, ymax-self.yllcorner]:
            if arg%self.cellsize != 0:
                print 'Error: clipping limits must align with raster origin/cellsize'
                return
        
        #  3. must have max > min for x and y
        if xmin >= xmax or ymin >= ymax: 
            print 'Error: check clipping limits - maxs must be greater than mins'
            return
        
        # rationalise clipping limits if they lie outside raster limits
        gxmax, gymax = self.xllcorner + self.ncols*self.cellsize, self.yllcorner + self.nrows*self.cellsize
        if xmin < self.xllcorner or ymin < self.yllcorner or xmax > (gxmax) or ymax > (gymax):
            print 'Warning: some clipping limits lie outside raster extents - defaulting these to raster extents.'
            xmin, ymin = max(xmin, self.xllcorner), max(ymin, self.yllcorner)
            xmax, ymax = min(xmax, gxmax), min(ymax, gymax)
        
        # calculate new indices, clip raster
        ixmin, iymin = (xmin - self.xllcorner)/self.cellsize, (ymin - self.yllcorner)/self.cellsize
        ixmax, iymax = (xmax - self.xllcorner)/self.cellsize, (ymax - self.yllcorner)/self.cellsize
        ytop = self.data.shape[0]
        self.data = self.data[ytop-iymax:ytop-iymin, ixmin:ixmax]
        
        # reset llcorner coords, ncols, nrows
        if xmin > self.xllcorner: self.xllcorner = xmin
        if ymin > self.yllcorner: self.yllcorner = ymin
        self.ncols = (xmax - xmin)/self.cellsize
        self.nrows = (ymax - ymin)/self.cellsize
    
    
    def write_to_xyz(self, filename):
        
        """Writes AsciiGrid object to an .xyz format file."""
        
        # correct filename extension if necessary
        if not filename.endswith('.xyz'): filename.append('.xyz')
        
        # write points to xyz
        xl,yl,cs = self.xllcorner, self.yllcorner, self.cellsize
        it = np.nditer(self.data, flags=['multi_index'], op_flags=['readonly'])
        with open(filename, 'wb') as fh:
            while not it.finished:
                iy,ix = it.multi_index
                fh.write(''.join([str(yl+iy*cs),' ',str(xl+ix*cs),' ',str(it[0]),'\n']))
                it.iternext()
                
                
    def value_at(self, x, y):
   
        """Return value at (x,y), specified in grid coordinates."""
       
        # check inputs
        if not (isinstance(x,int) and isinstance(y,int)):
            raise ValueError('x, y args must be integers: {},{}'.format(x,y))
           
        # calculate local coords: subtract xll/yll and divide by cellsize
        xl,yl = x-self.xllcorner, y-self.yllcorner
        xl /= self.cellsize
        yl /= self.cellsize
       
        # check if outside raster
        try:
            return self.data[self.nrows-yl-1, xl-1]
        except IndexError:
            print 'Coordinates lie outside raster area'
                
                
    
# test client
if __name__ == "__main__":
    
    grid = AsciiGrid()
    data = np.vstack([np.arange(1,6), np.arange(6,11), np.arange(11,16),
                      np.arange(16,21), np.arange(21,26)])
    grid2 = AsciiGrid(xllcorner=10, yllcorner=10, cellsize=5, data=data)
    grid3 = AsciiGrid.fromfile('example.asc')
    grid4 = AsciiGrid().read('example2.asc')

