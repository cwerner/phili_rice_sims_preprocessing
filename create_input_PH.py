#!/usr/bin/env python 

import luigi
import os
import numpy as np
import xarray as xr

# source tasks
from general import MakeDirectory 
from task_management import ProcessManagement
from task_climdata import ExtractClimData
from task_etopo1 import CreateElevationLayers
from task_ids import CreateAllBaseLayers
from task_climdata import DownsampleClimateData
from task_ricearea import ComputeRicePercentage

from luigi.contrib.external_program import ExternalProgramTask


LR = 0.5
HR = 0.08333333
BBOX = [116,4,127,22]


class CreateMiscFile(luigi.Task):
    
    inv = luigi.Parameter(default='PH')    

    def requires(self):
        yield MakeDirectory(path='downloads')
        yield MakeDirectory(path='tmp')
        yield ProcessManagement(bbox=BBOX)   # get management
        #yield ExtractClimData(bbox=[116,4,134,22], syear=1990, eyear=2012)   # get management
        yield CreateElevationLayers(bbox=BBOX)
        yield CreateAllBaseLayers(bbox=BBOX)
        yield ComputeRicePercentage(bbox=BBOX) 

    def run(self):
        file_names = ['tmp/' + x for x in ['manamask.nc', 'base_layers.nc', 'cids.nc', 'climids.nc', 'elevation.nc', 'manaids.nc', 'RiceArea_HR.nc']]
        print(file_names)
        files = [xr.open_dataset(x).transpose('lat','lon') for x in file_names]

        out = files[0].copy(deep=True)
        for file in files[1:]:
            vars = file.data_vars
            print(vars)

            for var in vars:
                out[var] = xr.DataArray(file[var].values, coords=out.coords).where(out['region'] == 1) #[('lat', out.coords['lat']), ('lon', out.coords['lon'])], dims=['lat', 'lon'])

        out['simmask'] = out['ricemask'] * out['region']
        out.to_netcdf(self.output().path, format='NETCDF4_CLASSIC')

    def output(self):
        file_path = f'data/{self.inv}_MISC.nc'
        return luigi.LocalTarget(file_path)   

class CreateSiteFile(ExternalProgramTask):
    bbox     = luigi.ListParameter()

    file_path = 'dummy.xml'

    # TODO: install ldndctools into virtual environment
    def program_args(self):
        return ["/Users/cwerner/Dropbox/development/ldndc_raster_tools/DLSC.py", "-r HR", f"-b {BBOX}", "--region=102"]

    def output(self):
        return luigi.LocalTarget(self.file_path)


class ProcessAll(luigi.Task):

    def requires(self):
        yield CreateMiscFile(inv='PH')
    
    def run(self):
        yield DownsampleClimateData(data_nc='data/PH_MISC.nc')
        
        # TODO: fix DLSC/ NLCC install woes
        #yield CreateSiteFile(bbox=BBOX)


if __name__ == '__main__':
    luigi.run(['ProcessAll', '--workers', '2', '--local-scheduler'])
