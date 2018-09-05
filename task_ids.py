# Id Task
# ===============

import luigi

import numpy as np
import xarray as xr

from general import bbox_to_latlon

LR = 0.5
HR = 0.083333333333333


# LR grid
a_lr = np.array(range(720*360)).reshape((360,720))
l_climids = []
for j in range(a_lr.shape[0]):
    for i in range(a_lr.shape[1]):
        l_climids.append(j*1000+i)
a_climids = np.array(l_climids).reshape(a_lr.shape)

# HR grid
sf = int(LR / HR)
a_climids_HR = np.kron(a_climids, np.ones((sf, sf)))
da_climids_HR = xr.DataArray(np.flipud(a_climids_HR), coords=[('lat', np.arange(-90+HR*0.5,90, HR)), ('lon', np.arange(-180+HR*0.5,180, HR))], 
                                                      dims=['lat','lon'], name='climid')


class CreateCid(luigi.Task):
    bbox = luigi.ListParameter()
    res  = luigi.FloatParameter(default=HR)
    path = luigi.Parameter(default='tmp')    

    file_path = 'tmp/cids.nc'

    #def requires(self):
    #    yield FixGridETopo1(bbox=self.bbox)

    def run(self):
        lats_, lons_ = bbox_to_latlon(self.bbox, self.res)

        # create cin variable
        DIM0 = len(lats_)
        DIM1 = len(lons_)
        # TODO: check if ids start with 0 or with 1 ???
        data = np.flipud( (np.arange(DIM0*DIM1) ).reshape((DIM0, DIM1)) )
        da = xr.DataArray(data, coords=(lats_, lons_), dims=('lat', 'lon'), name='cid')
        da.attrs['units'] = 'id'
       
        da.to_dataset().to_netcdf(self.output().path, format='NETCDF4_CLASSIC')

    def output(self):
        return luigi.LocalTarget(self.file_path)    


class CreateClimid(luigi.Task):
    bbox = luigi.ListParameter()
    res  = luigi.FloatParameter(default=HR)
    path = luigi.Parameter(default='tmp')    

    file_path = 'tmp/climids.nc'

    #def requires(self):
    #    yield FixGridETopo1(bbox=self.bbox)

    def run(self):
        lats_, lons_ = bbox_to_latlon(self.bbox, self.res)

        da = da_climids_HR.sel(lat=slice(min(lats_), max(lats_)), lon=slice(min(lons_), max(lons_)))
        da.attrs['units'] = 'id'
       
        da.to_dataset().to_netcdf(self.output().path, format='NETCDF4_CLASSIC')

    def output(self):
        return luigi.LocalTarget(self.file_path)   


class CreateIdLayers(luigi.WrapperTask):
    bbox = luigi.ListParameter()

    def requires(self):
        yield CreateCid(bbox=self.bbox)
        yield CreateClimid(bbox=self.bbox)