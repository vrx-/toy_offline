import numpy as np
import netCDF4
from datetime import datetime
import os

def make_frc(rootdir='../project/',
             his_path='ocean_his.nc',
             frc_name='ocean_frc.nc'):
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    
    vars = {
        'swrad': {
            'long_name':'solar shortwave radiation flux', 
            'dims':('srf_time', 'eta_rho', 'xi_rho')
            }, 
        'shflux': {
            'long_name':'surface net heat flux', 
            'dims':('shf_time', 'eta_rho', 'xi_rho')
            },
        'bhflux': {
            'long_name':'bottom net heat flux',
            'units': 'watt meter-2',
            'field':'bottom heat flux, scalar, series', 
            'dims':('bhf_time', 'eta_rho', 'xi_rho')
            },
        'sustr': {
            'long_name':'surface u-momentum stress',
            'units': 'newton meter-2',
            'field':'surface u-momentum stress, scalar, series', 
            'dims':('sms_time', 'eta_u', 'xi_u')
            },
        'svstr': {
            'long_name':'surface v-momentum stress',
            'units': 'newton meter-2',
            'field':'surface v-momentum stress, scalar, series', 
            'dims':('sms_time', 'eta_v', 'xi_v')},
        }
                
    frc = netCDF4.Dataset(os.path.join(rootdir, frc_name), 'w', format='NETCDF3_CLASSIC')
    frc.Description = 'Forcing for offline simulation'
    cfrcAuthor = 'VRX'
    frc.Created = datetime.now().isoformat()
    frc.type = 'ROMS FRC file'
    
    his = netCDF4.Dataset(his_path)
            
    def write_nc(var):
        for dim in vars[var]['dims']:
            if dim not in frc.dimensions.keys():
                try:
                    frc.createDimension(dim, his.dimensions[dim].size)
                except:
                    frc.createDimension(dim, his.dimensions['ocean_time'].size)
                    frc.createVariable(dim, 'd', (dim,))
                    frc[dim][:] = his['ocean_time'][:]
                    frc[dim].long_name = "time for "+ vars[var]['long_name']
                    frc[dim].units = his['ocean_time'].units
                    frc[dim].calendar =his['ocean_time'].calendar
                    frc[dim].field = dim +", scalar, series"
        frc.createVariable(var, 'f8', vars[var]['dims'])
        if var in his.variables.keys():
            frc[var][:] = his[var][:]
            frc[var].long_name = his[var].long_name
            try:
                frc[var].units = his[var].units
            except:
                pass
            try:
                frc[var].field = his[var].field
            except:
                pass
        else:
            frc[var][:] = 0.
            frc[var].long_name = vars[var]['long_name']
            frc[var].units = vars[var]['units']
            frc[var].field = vars[var]['field']
            frc[var].time = vars[var]['dims'][0]  
          

    for var in vars:
        write_nc(var)
        
    frc.close()
    his.close()
    
if __name__ == '__main__':

    make_frc(rootdir=os.getcwd())
