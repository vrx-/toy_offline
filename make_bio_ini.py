import netCDF4 as nc
import pandas as pd
import numpy as np
import scipy.stats as stats
import os
import subprocess
from datetime import datetime
import time


def make_ini(rootdir='../project/',
            cdl_path='CDL/ini_fennel.cdl',
            bio_path= 'biovars.csv',
            clm_path='ocean_clm.nc', 
            output_path='bio_ini.nc',
            Is=None, Js=None, gauss=True, std=3.98942):
            
    command = 'ncgen ' + cdl_path + ' -o ' + output_path
    print(command)
    subprocess.call(command, shell=True)
    while not os.path.exists(output_path):
        time.sleep(1)
    ibio = nc.Dataset(output_path, 'a')
    clm = nc.Dataset(clm_path)
    
    ibio.Author = 'VRX'
    ibio.Created = datetime.now().isoformat()
    ibio.type = 'ROMS INI file'       
    ibio['ocean_time'].units = 'seconds since    1-01-01 00:00:00'
    ibio['ocean_time'].calendar = 'gregorian_proleptic'
    ibio['ocean_time'].field = 'time, scalar, series'
    ibio.variables['ocean_time'][0] = 0.0
    ibio['Vtransform'][:] = 2
    ibio['Vstretching'][:] = 4
    ibio['theta_s'][:] = 3.0
    ibio['theta_b'][:] = 0.0    
    ibio['Tcline'][:] = 25.0
    ibio['hc'][:] = 25.0
    ### from clm
    ibio['zeta'][0,:,:]  = clm['zeta'][0,:,:]
    ibio['ubar'][0,:,:] = clm['ubar'][0,:,:]
    ibio['vbar'][0,:,:] = clm['vbar'][0,:,:]
    ibio['u'][0,:,:,:] = clm['u'][0,:,:,:]
    ibio['v'][0,:,:,:] = clm['v'][0,:,:,:]
    ibio['temp'][0,:,:,:] = clm['temp'][0,:,:,:]
    ibio['salt'][0,:,:,:] = clm['salt'][0,:,:,:]
    
    if gauss is True:
        eta_rho, xi_rho = ibio.dimensions['eta_rho'].size, ibio.dimensions['xi_rho'].size
        fill_matrix = np.tile(10*stats.norm.pdf(np.arange(xi_rho), .5*xi_rho, std), (eta_rho, 1))
    else:
        if all is True:
            fill_matrix = np.ones(ibio['temp'].shape)
        else:
            fill_matrix = np.zeros(ibio['temp'].shape)
            if Is is None:
                fill_matrix[:, :, Js, :] = 1
            elif Js is None:
                fill_matrix[:, :, :, Is] = 1
            else:
                for i,j in zip(Is,Js):
                    fill_matrix[:, :, j, i] = 1

    ### from bio csv 
    df = pd.read_csv('biovars.csv')
    for key in df.keys():
        ibio[key][:] = df[key].values[np.newaxis, :, np.newaxis, np.newaxis] * fill_matrix


    ibio.close()
if __name__ == '__main__':

    make_ini(rootdir=os.getcwd())
