import netCDF4 as nc
import pandas as pd
import numpy as np
import os
import subprocess
from datetime import datetime
import time


def make_nud(rootdir='../project/',
            cdl_path='CDL/nudgcoef.cdl', 
            output_path='ocean_nud.nc'):
            
    command = 'ncgen ' + cdl_path + ' -o ' + output_path
    print(command)
    subprocess.call(command, shell=True)
    while not os.path.exists(output_path):
        time.sleep(1)
    nud = nc.Dataset(output_path, 'a')
    
    nud.Author = 'VRX'
    nud.Created = datetime.now().isoformat()
    nud.type = 'ROMS NUD file'   

    cf = 1440
    nud['M2_NudgeCoef'][:] = cf
    nud['M3_NudgeCoef'][:] = cf
    nud['temp_NudgeCoef'][:] = cf
    nud['salt_NudgeCoef'][:] = cf
    nud['tracer_NudgeCoef'][:] = 0
    nud.close()
    
if __name__ == '__main__':

    make_nud(rootdir=os.getcwd())

