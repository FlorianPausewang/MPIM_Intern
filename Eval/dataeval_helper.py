import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
import dask

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import cmocean

import cartopy.crs as ccrs
import cartopy.feature as cf

#=====================================================
#Notes

'''
Names:
da = dataarray
ds = dataset, usually a dataarray works too

'''

#=====================================================
#Select data  ---  Mostly redundant
'''
get_sea_level
get_surface level --> dataload_helper
get_vertical_levels


get_lat_zone
get_double_lat_zone

get_timeframe
get_timepoint

elevation2nan

extract below temp inv
'''

def get_sea_level(ds):
    return ds.isel(p=0)

def get_vertical_levels(ds,p1,p2):
    return ds.sel(p=slice(p1,p2))
    
def get_lat_zone(ds,l1,l2):
    if l1<l2:    
        min_lon=l1
        max_lon=l2
    else:
        min_lon=l2
        max_lon=l1
    #print(min_lon,max_lon)
    return ds.where((ds.lat >= min_lon) & (ds.lat <= max_lon))
    
    
def get_double_lat_zone(ds,l1,l2):
    return ds.where((np.abs(ds.lat) >= l1) & (np.abs(ds.lat) <= l2))
    
    
def get_timeframe(ds,t1,t2):
    return ds.sel(time=slice(t1,t2))

def get_timepoint(ds,t1):
    return ds.sel(time=t1)
    

def elevation2nan(ds,ds_2d):
    return ds.where( (ds.p < ds_2d.SP).compute() ,drop=True)

def extract_below_temp_inv(ds):
    with ProgressBar():
             ds_below_min = dask.compute( ds.where(ds.p >= ds['T'].idxmin(dim='p', keep_attrs=True)  )   )[0]   # ds.T.notnull().idxmax(dim="p")).max("p"))[0] 
    return ds_below_min   

#=====================================================
#Take mean functions over dataarray  ---  Mostly redundant
'''
time_mean
year_mean
long_mean
lat_mean
p_mean
'''


def mean_over_da(da,mean_axis):
    return da.mean(dim=mean_axis)
    
def time_mean(da,  mean_axis=['time']):  
    return mean_over_da(da,mean_axis)

def year_mean(da):
    return time_mean(da.groupby('time.year'))

def lon_mean(da, mean_axis=['lon']):
    return mean_over_da(da,mean_axis)

## use weighted mean for lat_mean()
def lat_weighting(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name="weights"
    weighted_data = data.weighted(weights)
    weighted_data.attrs = data.attrs
    return weighted_data

def lat_mean(da, mean_axis=['lat']):
    return mean_over_da(lat_weighting(da),mean_axis)

def p_mean(da, mean_axis=['p']):
    return mean_over_da(da,mean_axis)


#=====================================================
#Unit conversion

def K2Deg(da):
    da.attrs["units"] = "deg C"
    return (da - 273.15)
    
#=====================================================
#Plot Stuff 

plt.rcParams.update(plt.rcParamsDefault)
'''
create_empty_map


'''

def create_empty_2d_map(cl=True,grid=True,projection=None):
    # Inspired by https://medium.com/@lubomirfranko/climate-data-visualisation-with-python-visualise-climate-data-using-cartopy-and-xarray-cf35a60ca8ee
    # First we specify Coordinate Refference System for Map Projection
    # We will use Mercator, which is a cylindrical, conformal projection. 
    # It has bery large distortion at high latitudes, cannot 
    # fully reach the polar regions.
    if projection is None:
        projection = ccrs.PlateCarree()
    # Specify CRS, that will be used to tell the code, where should our data be plotted
    
    # Now we will create axes object having specific projection 
    fig= plt.figure(figsize=(10,6), dpi=300)#
    ax = plt.axes(projection=projection, frameon=True)
    
    if grid==True:
        # Draw gridlines in degrees over Mercator map
        gl = ax.gridlines( draw_labels=True,
                        linewidth=.5,color='gray', alpha=0.5, linestyle='--') #crs=projection,
        gl.top_labels = False
        gl.right_labels = False
        #gl.xlines = False
        gl.ylocator = mticker.FixedLocator([-90,-70,-50,-30,-10,0,10,30,50,70,90])
        #gl.yformatter = LatitudeFormatter()
        gl.xlabel_style = {"size" : 6}
        gl.ylabel_style = {"size" : 6}
        
    if cl==True:
        #ax.coastlines()
        # To plot borders and coastlines, we can use cartopy feature
        ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=0.5)
    
    ''' For specific extend add
    lon_min = -20
    lon_max = 45
    lat_min = 34
    lat_max = 60
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=crs)
    '''
    
    cbar_kwargs = {'shrink':0.6}#'orientation':'vertical', 'shrink':0.6, "pad" : .05, 'aspect':40
    
    ##### ADD THESE LINES IN SCRIPT#####
    
    #dataset["t2m"].plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs=cbar_kwargs, levels=21)
    
    #plt.title(f"Temperature anomaly over Europe in {dataset.valid_time.dt.strftime('%B %Y').values}")
    #plt.show()
    ################################
    return fig, ax, gl, cbar_kwargs

    
# Define the custom colormaps
cmap_T = plt.get_cmap(cmocean.cm.balance, 30)

colors_RH_change = ['#964B00', 'white', '#006400']
cmap_RH_change = mcolors.LinearSegmentedColormap.from_list('brown_green_cmap', colors_RH_change)

colors_RH_change_gray = ['#964B00', 'lightgray', '#006400']
cmap_RH_change_gray = mcolors.LinearSegmentedColormap.from_list('brown_green_cmap', colors_RH_change_gray)


colors_RH = ['#964B00', '#D2B48C', '#90EE90', '#006400']#['#964B00', '#7F853D', '#5A9E63', '#37B58A', '#14CDB2', '#00E7DB']
cmap_RH= mcolors.LinearSegmentedColormap.from_list('brown_green_cmap', colors_RH)

elevation_colors = ["darkblue", "lightblue", "white", "lightgreen", "darkgreen", "brown"]
cmap_elevation= mcolors.LinearSegmentedColormap.from_list("elevation_cmap", elevation_colors , N=256)

TOP_LIM4PLOT=10













'''
Trash!!!
'''

'''
def get_highest_pressure(ds):
    # Create a mask to identify NaN values in the temperature dimension
    nan_mask = np.isnan(ds.T)
    #print(nan_mask.to_numpy())
    # Apply the mask to the min_pressure_index to ignore NaN values
    ds_masked = ds.where(nan_mask, drop = False)
    print(ds_masked)

    # Find the index of the minimum pressure (maximum height) along the 'p' dimension
    ds_min_pressure_index = ds_masked.p.argmax(dim='p')
    
    print(ds_min_pressure_index.to_numpy())
     
    # Use the index to extract the values corresponding to the highest pressure
    highest_pressure_data = ds.sel(p=ds['p'][ds_min_pressure_index])
    highest_pressure_data.attrs = ds.attrs
    
    del nan_mask, ds_masked, ds_min_pressure_index
    return highest_pressure_data   
    


def get_surface_level(ds, var = 'T'):
    
    mask= ds[var].notnull()
    #print(mask.to_numpy())
    # Create a new dataset with the highest pressure values
    new_dataset = p_mean(ds)
    
    # Extract the longitude, latitude, and time dimensions in for loops
    for lon_val in ds['lon']:
        print(lon_val)
        for lat_val in ds['lat']:
            for time_val in ds['time']:
                isnan=mask.sel(lon=lon_val, lat=lat_val, time=time_val)
                i=0
                while isnan[i]==False:
                    i=i+1
                    if isnan[i]==True:
                        # Access the climate data at the current combination of coordinates and time
                        new_dataset[var].loc[dict(lon=lon_val, lat=lat_val, time=time_val)]  = ds[var].sel(lon=lon_val, lat=lat_val, time=time_val).isel(p=i)

    return new_dataset

def get_surface_level(ds, var = 'T'):                                            ###########################################################################  SOS HELP  
    # Select the highest pressure level where pressure is not NaN
    max_pressure = ds.where(~np.isnan(ds[var])).max(dim='p')
    #max_pressure= ds.argmax(dim='p')#,skipna=True
    return ds.sel(p=max_pressure)

'''