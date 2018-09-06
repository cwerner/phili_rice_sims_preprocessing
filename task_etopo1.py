# ETopo1 Task
# =========== 

import io
import luigi
from luigi.contrib.external_program import ExternalProgramTask
import urllib
import numpy as np
import xarray as xr
import zipfile

from general import add_attributes

class DownloadETopo1(luigi.Task):
    etopo1_location = "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip"
    file_path = "downloads/ETOPO1_Ice_c_geotiff.tif"

    def run(self):
        r = urllib.request.urlopen(self.etopo1_location)
        z = zipfile.ZipFile(io.BytesIO(r.read()))
        z.extractall(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path)

class PatchETopo1(luigi.Task):
    file_path = 'tmp/ETOPO1_Ice_c.nc'

    def requires(self):
        yield DownloadETopo1()

    def run(self):
        da = xr.open_rasterio(self.input()[0].path).squeeze()
        da = da.rename({'y':'lat'})
        da['lat'].data = da['lat'].data

        da = da.rename({'x':'lon'})
        da['lon'].data = da['lon'].data

        da = add_attributes(da)
        
        da.to_dataset(name="relief").to_netcdf(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path)

class ResampleETopo1(luigi.contrib.external_program.ExternalProgramTask):
    file_path = 'tmp/ETOPO1_HR_raw.nc'

    bbox     = luigi.ListParameter()
    # TODO: clip to bbox
    
    def requires(self):
        yield PatchETopo1()

    def program_args(self):
        return ["cdo", "-selname,relief", "-remapbil,r4320x2160", self.input()[0].path, self.output().path]

    def output(self):
        return luigi.LocalTarget(self.file_path)


class ReorientETopo1(luigi.contrib.external_program.ExternalProgramTask):
    file_path = 'tmp/ETOPO1_HR_raw2.nc'

    bbox     = luigi.ListParameter()
    # TODO: clip to bbox
    
    def requires(self):
        yield ResampleETopo1(bbox=self.bbox)

    def program_args(self):
        return ["cdo", "-sellonlatbox,-180,180,-90,90", self.input()[0].path, self.output().path]

    def output(self):
        return luigi.LocalTarget(self.file_path)

class FixGridETopo1(luigi.contrib.external_program.ExternalProgramTask):
    file_path = 'tmp/ETOPO1_HR.nc'

    bbox     = luigi.ListParameter()
    # TODO: clip to bbox
    
    def requires(self):
        yield ReorientETopo1(bbox=self.bbox)

    def program_args(self):
        return ["cdo", "-setgrid,misc/grid.global_hr", self.input()[0].path, self.output().path]

    def output(self):
        return luigi.LocalTarget(self.file_path)


class CreateElevationLayers(luigi.Task):
    """Create elevation layers for final <VN,PH>_MISC.nc file
    """

    file_path = 'tmp/elevation.nc'

    bbox     = luigi.ListParameter(default=int)

    BASE_NC = 'data/climate/misc/elevation_05deg.nc'

    def requires(self):
        yield FixGridETopo1(bbox=self.bbox)

    def run(self):

        da_LR = xr.open_dataset(self.BASE_NC).squeeze()['data']
        da = xr.open_dataset(self.input()[0].path)['relief']


        scale_factor_lon = len(da.lon.values) / len(da_LR.lon.values)
        scale_factor_lat = len(da.lat.values) / len(da_LR.lat.values)

        scale_factor_lon = int(len(da.lon.values) / len(da_LR.lon.values))
        scale_factor_lat = int(len(da.lat.values) / len(da_LR.lat.values))
       
        assert scale_factor_lon == scale_factor_lat
        
        scale_factor = scale_factor_lon

        # scale LR elevation data to HR grid
        a_elevation_LR = np.kron(da_LR.values, np.ones((scale_factor, scale_factor)))

        # build new nc file with elevation, topo_diff vars
        da_elevation = xr.DataArray(da.values, coords=da.coords, dims=da.dims, name='elevation')
        da_topodiff  = xr.DataArray(da.values - a_elevation_LR, coords=da.coords, dims=da.dims, name='topodiff')

        # mask non-land elevation
        da_topodiff  = da_topodiff.where(~da_elevation.isnull())
        da_elevation = da_elevation.where(~da_topodiff.isnull())

        da_elevation.attrs['units'] = 'm'
        da_topodiff.attrs['units'] = 'm'

        ds = xr.Dataset()
        ds['elevation'] = da_elevation
        ds['topodiff'] = da_topodiff

        # clip to bounding box (maybe we should do this further up the chain)
        lon1, lat1, lon2, lat2 = self.bbox
        ds = ds.sel(lat=slice(lat1, lat2), lon=slice(lon1, lon2))
        ds.to_netcdf(self.output().path)

    def output(self):
        return luigi.LocalTarget(self.file_path)    