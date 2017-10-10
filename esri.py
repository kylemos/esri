"""
Class for handling ESRI ASCII-format rasters.
"""
from __future__ import print_function
import os.path as pth
import numpy as np
from weakref import WeakKeyDictionary

class NonNegative(object):
    
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
    
    Attributes:
        ncols (int):        Number of columns in AsciiGrid (read-only).
        nrows (int):        Number of rows in AsciiGrid (read-only).
        xllcorner (int):    X-coordinate of lower-left grid corner.
        yllcorner (int):    Y-coordinate of lower-left grid corner.
        cellsize (int):     Grid cell size (in metres).
        no_data_val (int):  No-data value (default = -9999).
        data (np.ndarray):  Raster array with no-data values represented as numpy.nan's.
    
    TODO:
        * Use decorators to enable argument-checking on methods? Write as both function & class to learn!
        * Create a context manager and __iter__ method? i.e. to do with AsciiGrid as x, loop over x.cols???
        * See https://github.com/olemb/dbfread/blob/master/dbfread/dbf.py for example!
        * Option in write() method to specify xllcorner vs. xllcenter
        * Option in to_xyz() method to specify precision?
        * Implement some way to compare grids (i.e. __eq__ etc...).
        * Implement a 'change_no_data_val' method:
          --> Can probably do this with a @property or something similar? So just 'step in' when
              the no_data_val attribute is changed and change it in the data attribute too?
        * Replace add_grid() method with __add__()?
        * REPLACE value_at(x,y) METHOD TO __getitem__() ?
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
        # (protected attribute)
        return self._no_data_val
        
    @no_data_val.setter
    def no_data_val(self, value):
        if not isinstance(value, int):
            raise ValueError('no_data_val must be an integer value: {}'.format(value))
        self._no_data_val = value
        
    @property
    def data(self):
        return self.__data
        
    @data.setter
    def data(self, arr):
        if not isinstance(arr, np.ndarray):
            raise TypeError('data attribute must be a numpy.ndarray object')
        self.__data = np.copy(arr)
    
    
    def __init__(self, xllcorner=0, yllcorner=0, cellsize=1, 
                 no_data_val=-9999, data=np.zeros((10,10))):
                 
        """
        Creates an AsciiGrid raster object from kwargs and a numpy array. 
        Defaults to a 10x10 grid of zeros.
        """
        
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.no_data_val = no_data_val
        self.data = np.where(data == no_data_val, np.nan, data)
        
    
    def __repr__(self):
        
        return '{}'.format(self.__class__.__name__)
        
    
    def __str__(self):
        
        s1 = '{0}x{1} AsciiGrid'.format(self.ncols, self.nrows)
        s2 = 'xll: {0}, yll: {1}'.format(self.xllcorner, self.yllcorner)
        s3 = 'cellsize: {0}'.format(self.cellsize)
        return ', '.join([s1,s2,s3])
        
    
    def read(self, filename):
    
        """Reads an ascii grid from the filename (deprecated; maintained for backwards compatibility)."""
        
        return AsciiGrid.fromfile(filename)
        
    
    @classmethod
    def fromfile(cls, filename):
    
        """Reads an ESRI ASCII-format raster from file."""
        
        # read file headers, data array
        try:
            fh = open(filename, 'r')
        except IOError:
            print('File not found/error reading file')
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
            print('Missing headers in file')
            raise
        if d_hd:
            raise TypeError('Unexpected header entries found in file: {}'.format(d_hd.keys()))
        
        # check data, create instance
        try:
            data = np.reshape(arr, (nrows,ncols))
        except ValueError:
            print('nrows, ncols inconsistent with raster array shape')
            raise
        else:
            return cls(xllcorner=xllcorner, yllcorner=yllcorner, cellsize=cellsize, 
                       no_data_val=no_data_val, data=data)
            
    
    def write(self, filename, precision=2):
    
        """Writes the raster to ESRI ASCII-format file specified in first argument."""
        
        # check filename arg
        if not pth.isdir(pth.dirname(filename)):
            raise IOError('filepath not recognised: {}'.format(filename))
        
        # check precision arg
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
        
        # write header, raster to file
        np.savetxt(filename, self.data, delimiter=' ', header=header, fmt=fmt, comments='')
        
    
    def add_grid(self, newAsciiGrid):
    
        """Adds an AsciiGrid raster."""
        
        # TODO: check same xllcorner, yllcorner, shape, etc...
        
        # error checking
        # if newAsciiGrid.shape != self.data.shape:
            # print('Error: grid must have same shape as Ascii raster grid.')
            # raise Exception
    
        # try:
        self.data += newAsciiGrid.data
        # except <some_error>:
            # print('error: cannot add grids')
        
    
    def clip(self, xmin=None, xmax=None, ymin=None, ymax=None):
    
        """Clips the raster to specified limits (xmin, xmax, ymin & ymax)."""
        
        # sanity-check clipping limits:
        # integer args
        if not all(isinstance(arg, int) for arg in [xmin,xmax,ymin,ymax]):
            raise ValueError('Non-numerical value specified in clipping limits')
        
        # align with grid (i.e. be a multiple of cellsize removed from llcorner)
        for arg in [xmin-self.xllcorner, xmax-self.xllcorner, ymin-self.yllcorner, ymax-self.yllcorner]:
            if arg % self.cellsize != 0:
                raise  ValueError('Clipping limits must align with raster origin & cellsize')
        
        # max > min for both x and y
        if xmin >= xmax or ymin >= ymax:
            raise ValueError('Inconsistent clipping limits: max must be greater than min')
        
        # rationalise clipping limits if they lie outside raster limits
        gxmax, gymax = self.xllcorner + self.ncols*self.cellsize, self.yllcorner + self.nrows*self.cellsize
        if xmin < self.xllcorner or ymin < self.yllcorner or xmax > (gxmax) or ymax > (gymax):
            print('Warning: some clipping limits lie outside raster extents -> defaulted to raster extents.')
            xmin, ymin = max(xmin, self.xllcorner), max(ymin, self.yllcorner)
            xmax, ymax = min(xmax, gxmax), min(ymax, gymax)
        
        # calculate new indices, clip raster
        ixmin, iymin = (xmin - self.xllcorner)/self.cellsize, (ymin - self.yllcorner)/self.cellsize
        ixmax, iymax = (xmax - self.xllcorner)/self.cellsize, (ymax - self.yllcorner)/self.cellsize
        ytop = self.data.shape[0]
        self.data = self.data[ytop-iymax:ytop-iymin, ixmin:ixmax]
        
        # reset llcorner
        if xmin > self.xllcorner: self.xllcorner = xmin
        if ymin > self.yllcorner: self.yllcorner = ymin
    
    
    def to_xyz(self, filename):
        
        """Writes AsciiGrid object to an .xyz format file."""
        
        # correct file extension if required
        if not filename.endswith('.xyz'): 
            filename.append('.xyz')
        
        # write points to file
        xl,yl,cs = self.xllcorner, self.yllcorner, self.cellsize
        it = np.nditer(self.data, flags=['multi_index'], op_flags=['readonly'])
        with open(filename, 'w') as fh:
            while not it.finished:
                iy,ix = it.multi_index
                fh.write(' '.join([str(yl+iy*cs),str(xl+ix*cs),str(it[0]),'\n']))
                it.iternext()
                
    
    def value_at(self, x, y, interpolate=True):
   
        """
        Return value at (x,y), where x & y are specified in grid coordinates.
       
        Keyword arguments:
            interpolate - Set True to return bilinear interpolation in x & y
                          between 4 nearest data nodes.                       
        """
       
        # check inputs
        if not (isinstance(x,int) and isinstance(y,int)):
            raise ValueError('x, y args must be integers: {},{}'.format(x,y))
           
        # calculate local coords: subtract xll/yll and divide by cellsize
        xo,yo = x-self.xllcorner, y-self.yllcorner
        xl = int(math.floor(xo/self.cellsize))
        yl = int(math.floor(yo/self.cellsize))
       
        # calculate fractional
        if interpolate:
            xf,yf = xo-xl, yo-yl
       
        # check if outside raster area
        try:
            ll = self.data[self.nrows-(yl+1), xl  ]
            lr = self.data[self.nrows-(yl+1), xl+1]
            ul = self.data[self.nrows-yl,     xl  ]
            ur = self.data[self.nrows-yl,     xl+1]
        except IndexError:
            print('Coordinates lie outside raster area')
            return None
       
        # return (and interpolate if requested)
        if interpolate:
            return AsciiGrid.interpolate(ll,lr,ul,ur,xf,yf)
        else:
            return ll
           
        
    @staticmethod
    def interpolate(ll,lr,ul,ur,xf,yf):
   
        """bilinear interpolation..."""
            
    
    # def __eq__(self, other):
    
        # if isinstance(other, self.__class__):
            # eq = (np.array_equal(self.data, other.data) and
                  # self.ncols == other.ncols             and 
                  # self.nrows == other.nrows             and
                  # self.xllcorner == other.xllcorner     and
                  # self.yllcorner == other.yllcorner     and
                  # self.cellsize == other.cellsize       and
                  # self.no_data_val == other.no_data_val)
            # return eq
        # return NotImplemented
            
    
    # def __ne__(self, other):
    
        # return not self == other
                
    
# test client
if __name__ == "__main__":

    pass
    
    # grid = AsciiGrid()
    # data = np.vstack([np.arange(1,6), np.arange(6,11), np.arange(11,16),
                      # np.arange(16,21), np.arange(21,26)])
    # grid2 = AsciiGrid(xllcorner=10, yllcorner=10, cellsize=5, data=data)
    # grid3 = AsciiGrid.fromfile('example.asc')
    # grid4 = AsciiGrid().read('example2.asc')

