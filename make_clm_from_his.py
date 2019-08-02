import numpy as np
import netCDF4
from datetime import datetime
import os

def make_clm(rootdir='../project/',
             his_name='ocean_his.nc',
             clm_name='ocean_clm.nc'):
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)
    
    keep = ['zeta', 'ubar', 'vbar', 'u', 'v', 'omega', 'temp', 'salt', 'swrad', 'shflux']
    tnames = ['zeta_time', 'u2d_time', 'v2d_time', 'u3d_time', 'v3d_time', 
                'ocean_time', 'temp_time', 'salt_time', 'srf_time', 'ssf_time']
                
    clm = netCDF4.Dataset(os.path.join(rootdir, clm_name), 'w', format='NETCDF3_CLASSIC')
    clm.Description = 'CLimatology for offline simulation'
    clm.Author = 'VRX'
    clm.Created = datetime.now().isoformat()
    clm.type = 'ROMS CLM file'
    
    his = netCDF4.Dataset(os.path.join(rootdir, his_name))
    
    for dim in his.dimensions.keys():
        clm.createDimension(dim, his.dimensions[dim].size)
        
    
    #########################################
    times = his.dimensions['ocean_time'].size
    for var, tname in zip(keep, tnames):
        if tname not in clm.dimensions.keys():
            clm.createDimension(tname, times)
        clm.createVariable(tname, 'd', (tname,))
        clm[tname][:] = his['ocean_time'][:]
        clm[tname].long_name = "time for "+his[var].long_name
        clm[tname].units = his['ocean_time'].units
        clm[tname].calendar =his['ocean_time'].calendar
        clm[tname].field = tname +", scalar, series"
        clm.createVariable(var, 'f8', his[var].dimensions)
        clm[var][:] = his[var][:]
        clm[var].long_name = his[var].long_name
        try:
            clm[var].units = his[var].units
        except:
            print(var+' has no units')
        try:
            clm[var].field = his[var].field
        except:
            print(var+' has no field')
        clm[var].time = tname           
        
    clm.close()
    his.close()
    
if __name__ == '__main__':

    make_clm(rootdir=os.getcwd())
    

