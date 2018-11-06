import numpy as np
import netCDF4
from datetime import datetime
import os

def make_frc(rootdir='../project/',
             his_name='../toy/ocean_his.nc',
             frc_name='ocean_frc.nc'):
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)
    
    keep = ['swrad', 'shflux']
    tnames = ['srf_time', 'ssf_time']
    
    new ={'bhflux':{'long_name':'bottom net heat flux',
                    'units': 'watt meter-2',
                    'field':'bottom heat flux, scalar, series',
                    'time':'bhf_time'}}
                
    frc = netCDF4.Dataset(os.path.join(rootdir, frc_name), 'w', format='NETCDF3_CLASSIC')
    frc.Description = 'Forcing for offline simulation'
    cfrcAuthor = 'VRX'
    frc.Created = datetime.now().isoformat()
    frc.type = 'ROMS FRC file'
    
    his = netCDF4.Dataset(os.path.join(rootdir, his_name))
    
    for dim in his[keep[0]].dimensions:
        frc.createDimension(dim, his.dimensions[dim].size)
            
    def write_nc(var, tname):
        frc.createVariable(tname, 'd', (tname,))
        frc[tname][:] = his['ocean_time'][:]
        frc[tname].long_name = "time for "+his[var].long_name
        frc[tname].units = his['ocean_time'].units
        frc[tname].calendar =his['ocean_time'].calendar
        frc[tname].field = tname +", scalar, series"
        frc.createVariable(var, 'f8', (tname, 'eta_rho', 'xi_rho'))
        frc[var][:] = his[var][:]
        frc[var].long_name = his[var].long_name
        try:
            frc[var].units = his[var].units
        except:
            print(var+' has no units')
        try:
            frc[var].field = his[var].field
        except:
            print(var+' has no field')
        frc[var].time = tname  
          
    times = his.dimensions['ocean_time'].size
    for var, tname in zip(keep, tnames):
        if tname not in frc.dimensions.keys():
            frc.createDimension(tname, times)
        write_nc(var, tname)
    for var in new:
        tname = new[var]['time']      
        frc.createDimension(tname, times)
        frc.createVariable(new[var]['time'], 'd', (tname,))
        frc[tname][:] = his['ocean_time'][:]
        frc[tname].long_name = "time for "+new[var]['long_name']
        frc[tname].units = his['ocean_time'].units
        frc[tname].calendar =his['ocean_time'].calendar
        frc[tname].field = tname +", scalar, series"        
        frc.createVariable(var, 'f8', (tname, 'eta_rho', 'xi_rho'))
        frc[var][:] = 0.
        frc[var].long_name = new[var]['long_name']
        frc[var].units = new[var]['units']
        frc[var].field = new[var]['field']
        frc[var].time = new[var]['time']        
    frc.close()
    his.close()
    
if __name__ == '__main__':

    make_frc(rootdir=os.getcwd())
