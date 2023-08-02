import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
import dask.dataframe as dask
import os
import dask as dsk
#=====================================================
#Data class and functions to load and save the data

class Dataset:
    def __init__(self, files,setname,dataloc, frequency = 'monthly'):
        self.files = files
        self.setname = setname
        self.frequency = frequency
        self.dataloc = dataloc
        
        self.chunks = {'time': 'auto', 'lon': 'auto', 'lat': 'auto'}#"time": 10, "p": 42, "lon": 288, "lat": 144}
        self.loaddata()
    
    
    """
    #load data
    def loaddata(self):
        if self.setname=='MERRA':
            self.loaddata()
        elif self.setname=='MERRA2':
            self.loaddata()
        elif self.setname=='JRA-55':
            self.loaddata()
        else:
            raise Exception("Sorry, notexisting Setname")"""
    
    #load data routine for Merra 
    def loaddata(self):
        print("Loading full dataset "+ self.setname)
        #with dsk.config.set(**{'array.slicing.split_large_chunks': True}):
        with ProgressBar():
            self.ds = xr.open_mfdataset(self.files,parallel=True,engine='netcdf4',
                               ) #combine='by_coords', chunks=self.chunks, 
        #print(self.ds)
        if self.setname=="TOPO":
            rename_dict = {}
            self.ds=self.ds.rename(rename_dict)
            # Check the longitude coordinate values
            lon_values = self.ds['lon'].values
            # Check if the longitude values are in the 0 to 360 range
            if lon_values.max() > 180:
                self.ds = self.ds.assign_coords(lon=(self.ds['lon'] + 180) % 360 - 180)
                self.ds = self.ds.sortby('lon')
        elif self.setname=="MERRA" and self.frequency == 'monthly':
            rename_dict = {'TIME': 'time','XDim': 'lon', 'YDim': 'lat', 'Height': 'p', 'RH': 'RH', 'T': 'T','Var_RH': 'RH_Var', 'Var_T': 'T_Var'}
            self.ds=self.ds.rename(rename_dict)
            #ds.RH.attrs["units"] = '%'  #To rename Units
        elif self.setname=="MERRA2" and self.frequency == 'monthly':
            rename_dict = {'lev': 'p', 'RH': 'RH', 'T': 'T','PS':'PS'}
            self.ds=self.ds.rename(rename_dict)
        elif self.setname=="MERRA2_2D" and self.frequency == 'monthly':
            rename_dict = {}
            self.ds=self.ds.rename(rename_dict)
        elif self.setname=='JRA-55' and self.frequency == 'monthly':
            rename_dict = {'plev': 'p', 'r': 'RH', 't': 'T'}
            self.ds=self.ds.rename(rename_dict)
            self.ds=self.ds.isel(p=slice(None, None, -1))  # Reverse order of pressurelevels
            self.ds=self.ds.isel(lat=slice(None, None, -1))  # Reverse order of pressurelevels
            self.ds['RH']=0.01 *self.ds.RH  #Convert to unitless
            self.ds.RH.attrs["units"] = ''  #To rename Units
            self.ds['p']=0.01 *self.ds.p  #Convert to hPa
            self.ds.p.attrs["units"] = 'hPa'  #To rename Units
            
            # Check the longitude coordinate values
            lon_values = self.ds['lon'].values

            # Check if the longitude values are in the 0 to 360 range
            if lon_values.max() > 180:
                self.ds = self.ds.assign_coords(lon=(self.ds['lon'] + 180) % 360 - 180)
                self.ds = self.ds.sortby('lon')
        elif self.setname=='JRA-55_2D' and self.frequency == 'monthly':
            if 'sp' in list(self.ds.keys()):
                rename_dict = {'sp': 'SP'}
                self.ds=self.ds.rename(rename_dict)
                self.ds['SP']=0.01 *self.ds.SP  #Convert to hPa
                self.ds.SP.attrs["units"] = 'hPa'  #To rename Units
            if '2t' in list(self.ds.keys()):
                rename_dict = {'2t': 'T2M'}
                self.ds=self.ds.rename(rename_dict)
            if '2r' in list(self.ds.keys()):
                rename_dict = {'2r': 'RH2M'}
                self.ds=self.ds.rename(rename_dict)
                self.ds['RH2M']=0.01 *self.ds.RH2M  #Convert to unitless
                self.ds.RH2M.attrs["units"] = ''  #To rename Units
                
            self.ds=self.ds.isel(lat=slice(None, None, -1))  # Reverse order of pressurelevels
            # Check the longitude coordinate values
            lon_values = self.ds['lon'].values
            # Check if the longitude values are in the 0 to 360 range
            if lon_values.max() > 180:
                self.ds = self.ds.assign_coords(lon=(self.ds['lon'] + 180) % 360 - 180)
                self.ds = self.ds.sortby('lon')
             
            
        else:
            raise Exception("Sorry, notexisting Setname or Frequency") 
        
        if self.setname!="TOPO" and self.setname!='JRA-55_2D' and self.setname!='MERRA2_2D':
            self.ds['T']=self.ds.T-273.15
            self.ds.T.attrs["units"] = "deg C"
        
        print(self.ds)
        return self.ds
    
    
    
    
    # =======================================================================
    #Method to get surface values:   
    def get_surface_level(self, ds):
        # Select the highest pressure level where pressure is not NaN
        print("Extracting surface values")
        with ProgressBar():
             ds_surface = dask.compute( ds.where(ds.p == ds.T.notnull().idxmax(dim="p")).max("p"))[0] 
             #ds.where(ds.notnull()).max(dim='p') 
        return ds_surface
    
    def create_surfacevals(self):
        DS_new = self.get_surface_level(self.ds)
        DS_new.to_netcdf(path = self.dataloc+ 'surfacevals.netcdf')
        return DS_new

    def loadsurfacevals(self):
        return xr.open_mfdataset(self.dataloc+ 'surfacevals.netcdf', chunks=self.chunks)
    
    #Load Surfacevalues or create them
    def init_surfacevals(self):
        if os.path.exists(self.dataloc +'surfacevals.netcdf'):
            print("Load existing file for Surfacevals from "+self.dataloc +'surfacevals.netcdf')
            self.ds_surface=self.loadsurfacevals()
        else:
            self.ds_surface=self.create_surfacevals()
        print(self.ds_surface)
    
    def reloadsurfacevals(self):
        if os.path.exists(self.dataloc +'surfacevals.netcdf'):
            os.remove(self.dataloc+ 'surfacevals.netcdf')
        else:
            self.ds_surface=self.create_surfacevals()
    
    

    
    '''
    # =======================================================================
    #Wie geht das?
    #Clean up
    def clean():
        del self.ds
    '''

    

def loadsimpledata(file):
    ds = xr.open_dataset(file,engine='netcdf4',
                               ) #combine='by_coords', chunks=self.chunks, 
        
    return ds
    