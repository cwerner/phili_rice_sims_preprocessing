# Ricearea Task
# ============= 

import io
import luigi
from luigi.contrib.external_program import ExternalProgramTask
import urllib
import numpy as np
import xarray as xr
import zipfile
from affine import Affine

from general import add_attributes

HR = 0.08333333

def rebin(arr, new_shape, mode='mean'):
    """Mode: sum or mean"""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    if mode == 'sum':
        return arr.reshape(shape).sum(-1).sum(1)
    return arr.reshape(shape).mean(-1).mean(1) 



class DownloadRiceMap(luigi.Task):
    ricemap_location = "https://drive.google.com/uc?id=0BwTW1oKZU90LTUtqaDQ2XzZpbjg&export=download"
    file_path = "downloads/IRRI_Asia_Rice_Extent_Map"
    tiffile = 'asia-rice-extent.tif'

    def run(self):
        r = urllib.request.urlopen(self.ricemap_location)
        z = zipfile.ZipFile(io.BytesIO(r.read()))
        z.extractall(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path + '/' + self.tiffile)

class PatchRiceMap(luigi.Task):
    bbox     = luigi.ListParameter()

    file_path = 'tmp/IRRI_Asia_Rice_Extent.nc'

    def requires(self):
        yield DownloadRiceMap()

    def run(self):
        da = xr.open_rasterio(self.input()[0].path).squeeze()
        transform = Affine(*da.attrs['transform'][:6])
        nx, ny = da.sizes['x'], da.sizes['y']
        x, y = np.meshgrid(np.arange(nx)+0.5, np.arange(ny)+0.5) * transform

        da = da.rename({'y':'lat'})
        da['lat'].data = da['lat'].data

        da = da.rename({'x':'lon'})
        da['lon'].data = da['lon'].data

        da = add_attributes(da)

        # flip latitude orientation
        da = da.sel(lat=slice(None, None, -1))

        # clipt to bounding box
        da = da.sel(lat=slice(self.bbox[1], self.bbox[3]), lon=slice(self.bbox[0], self.bbox[2]))

        print(np.unique(da.values))

        da.to_dataset(name="ricearea").to_netcdf(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path)

class ResampleRiceMap(luigi.contrib.external_program.ExternalProgramTask):
    file_path = 'tmp/RiceArea_HR_raw.nc'

    bbox     = luigi.ListParameter()
    
    def requires(self):
        yield PatchRiceMap(bbox=self.bbox)

    def program_args(self):
        return ["cdo", "-remapbil,misc/grid.regional_sr", self.input()[0].path, self.output().path]

    def output(self):
        return luigi.LocalTarget(self.file_path)


class ComputeRicePercentage(luigi.Task):
    file_path = 'tmp/RiceArea_HR.nc'

    bbox     = luigi.ListParameter()
    # TODO: clip to bbox
    
    def requires(self):
        yield ResampleRiceMap(bbox=self.bbox)

    def run(self):
        da = xr.open_dataset(self.input()[0].path)['ricearea']

        data = rebin(da.values, (216,132), mode='sum')

        print(data.shape)        
        da_out = xr.DataArray(data, coords=[('lat', np.arange(self.bbox[1]+HR*0.5, self.bbox[3], HR)), 
                                            ('lon', np.arange(self.bbox[0]+HR*0.5, self.bbox[2], HR))],
                                    dims=['lat', 'lon'], name='ricearea')
        
        ds = xr.Dataset()
        ds['ricearea'] = da_out
        ds['ricemask'] = xr.where(da_out > 0, 1, np.nan)
        ds.to_netcdf(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path)