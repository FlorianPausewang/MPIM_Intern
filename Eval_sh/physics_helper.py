import os

import numpy as np
import xarray as xr
import typhon.physics as typ

from konrad import constants

import konrad



def sh2vmr(sh):
    sh=sh
    return sh/(1-sh)

def vmr2sh(vmr):
    vmr=vmr
    return vmr/(1+vmr)

def sh2rh_1D(ds):
    temp=ds.T+273.15
    p=temp['p']*100
    vmr= sh2vmr(ds.SH)
    da=konrad.physics.vmr2relative_humidity(vmr,p,temp)
    return da

def saturation_pressure_3D(temperature):
    r"""Return equilibrium pressure of water with respect to the mixed-phase.

    The equilibrium pressure over water is taken for temperatures above the
    triple point :math:`T_t` the value over ice is taken for temperatures
    below :math:`T_t–23\,\mathrm{K}`.  For intermediate temperatures the
    equilibrium pressure is computed as a combination
    of the values over water and ice according to the IFS documentation:

    .. math::
        e_\mathrm{s} = \begin{cases}
            T > T_t, & e_\mathrm{liq} \\
            T < T_t - 23\,\mathrm{K}, & e_\mathrm{ice} \\
            else, & e_\mathrm{ice}
                + (e_\mathrm{liq} - e_\mathrm{ice})
                \cdot \left(\frac{T - T_t - 23}{23}\right)^2
        \end{cases}

    References:
        IFS Documentation – Cy45r1,
        Operational implementation 5 June 2018,
        Part IV: Physical Processes, Chapter 12, Eq. 12.13,
        https://www.ecmwf.int/node/18714

    Parameters:
        temperature (float or ndarray): Temperature [K].

    See also:
        :func:`~typhon.physics.e_eq_ice_mk`
            Equilibrium pressure of water over ice.
        :func:`~typhon.physics.e_eq_water_mk`
            Equilibrium pressure of water over liquid water.

    Returns:
        float or ndarray: Equilibrium pressure [Pa].
    """

    e_eq_water = typ.e_eq_water_mk(temperature)
    e_eq_ice = typ.e_eq_ice_mk(temperature)
    
    is_water = temperature > constants.triple_point_water

    is_ice = temperature < (constants.triple_point_water - 23.0)
    
    #print(is_ice)
    
    e_eq = (
        e_eq_ice
        + (e_eq_water - e_eq_ice)
        * ((temperature - constants.triple_point_water + 23) / 23) ** 2
    )
    
    #e_eq = e_eq_ice.where(is_ice,e_eq_ice)
    #e_eq = e_eq_water.where(is_water,e_eq_water)
    
    e_eq = e_eq.where(~is_ice).where(~is_water)
    e_eq = e_eq.combine_first(e_eq_ice).combine_first(e_eq_water)
    
    return e_eq


def sh2rh_3D(ds):
    temp=ds.T+273.15
    p=temp['p']*100
    vmr= sh2vmr(ds.SH)
    e=saturation_pressure_3D(temp)
    return vmr*p/e.to_numpy()


def add_rh_calc_in_intervals(ds,dataloc):
    
    dates=ds.time.to_numpy()
    interval=12 #month
    print(len(dates)/interval)
    for i in range(int(len(dates)/interval)):
        print(str(i))
        start_date=dates[interval*i]
        end_date=dates[interval*(i+1)-1]
        
        tmp=sh2rh_3D(ds.sel(time=slice(start_date, end_date)))
        if i==0:
            da = tmp
        else:
            da = xr.concat([da, tmp], dim='time')
    da.attrs["description"] = "RH calculated from SH."
    da.name='RH_calc'
    da.to_netcdf(path = dataloc+ 'rh_calc_vals.netcdf')
    print(da)
    ds["RH_calc"]=da
    
    return ds

def load_vals(Dds):
    return xr.open_mfdataset(Dds.dataloc+ 'rh_calc_vals.netcdf')#, chunks=DS.chunks)

#Load values or create them
def init_rh_calc_vals(Dds):
    if os.path.exists(Dds.dataloc +'rh_calc_vals.netcdf'):
        print("Load existing file for rh_calc_vals from "+Dds.dataloc +'rh_calc_vals.netcdf')
        tmp=load_vals(Dds)
        
        Dds.ds=xr.merge([Dds.ds,tmp])
    else:
        Dds.ds=add_rh_calc_in_intervals(Dds.ds,Dds.dataloc)
    print(Dds.ds)
    return

def add_rh_calc(ds):
    
    rh_calc=sh2rh_3D(ds)
    ds["RH_calc"]=rh_calc
    ds["RH_calc"].attrs["description"] = "RH calculated from SH."
    return ds