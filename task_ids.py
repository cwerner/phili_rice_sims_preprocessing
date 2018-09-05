# Id Task
# ===============

import luigi

import numpy as np
import xarray as xr

from general import bbox_to_latlon, calc_area, add_attributes

LR = 0.5
HR = 0.083333333333333

sf = int(LR/HR)

class CreateBaseLayers(luigi.Task):
    bbox = luigi.ListParameter()
    res  = luigi.FloatParameter(default=HR)
    path = luigi.Parameter(default='tmp')    

    file_path = 'tmp/base_layers.nc'

    def run(self):
        lats_, lons_ = bbox_to_latlon(self.bbox, self.res)

        # create cin variable
        DIM0 = len(lats_)
        DIM1 = len(lons_)
        # TODO: check if ids start with 0 or with 1 ???
        data = np.flipud( (np.arange(DIM0*DIM1) ).reshape((DIM0, DIM1)) )
        da = xr.DataArray(data, coords=(lats_, lons_), dims=('lat', 'lon'), name='cellid')
        da.attrs['units'] = 'id'

        da = add_attributes(da)

        da_area = da.copy(deep=True)
        da_area[:] = 1
        for l in da_area.lat.values:
            da_area.loc[{'lat': l}] = calc_area(l, HR)

        ds = xr.Dataset()
        ds['cellid'] = da
        ds['area_ha'] = da_area       
        ds.to_netcdf(self.output().path, format='NETCDF4_CLASSIC')

    def output(self):
        return luigi.LocalTarget(self.file_path)    


class CreateCid(luigi.Task):
    bbox = luigi.ListParameter()
    res  = luigi.FloatParameter(default=HR)
    path = luigi.Parameter(default='tmp')    

    file_path = 'tmp/cids.nc'

    def run(self):

        a_hr = np.array(range(720*360*sf*sf)).reshape((360*sf,720*sf))
        l_cids = []
        for j in range(a_hr.shape[0]):
            for i in range(a_hr.shape[1]):
                l_cids.append(j*10000+i)
        a_cids_HR = np.array(l_cids).reshape(a_hr.shape)

        da_cids_HR = xr.DataArray(np.flipud(a_cids_HR), coords=[('lat', np.arange(-90+HR*0.5,90, HR)), ('lon', np.arange(-180+HR*0.5,180, HR))], 
                                                        dims=['lat','lon'], name='cid') 


        lats_, lons_ = bbox_to_latlon(self.bbox, self.res)

        da = da_cids_HR.sel(lat=slice(self.bbox[1], self.bbox[3]), lon=slice(self.bbox[0], self.bbox[2]))
        da.attrs['units'] = 'id'

        da = add_attributes(da)

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

        a_lr = np.array(range(720*360)).reshape((360,720))
        l_climids = []
        for j in range(a_lr.shape[0]):
            for i in range(a_lr.shape[1]):
                l_climids.append(j*1000+i)
        a_climids = np.array(l_climids).reshape(a_lr.shape)



        # HR grid
        a_climids_HR = np.kron(a_climids, np.ones((sf, sf)))
        da_climids_HR = xr.DataArray(np.flipud(a_climids_HR), coords=[('lat', np.arange(-90+HR*0.5,90, HR)), ('lon', np.arange(-180+HR*0.5,180, HR))], 
                                                            dims=['lat','lon'], name='climid')

        lats_, lons_ = bbox_to_latlon(self.bbox, self.res)

        da = da_climids_HR.sel(lat=slice(self.bbox[1], self.bbox[3]), lon=slice(self.bbox[0], self.bbox[2]))
        da.attrs['units'] = 'id'

        da = add_attributes(da)

        da.to_dataset().to_netcdf(self.output().path, format='NETCDF4_CLASSIC')

    def output(self):
        return luigi.LocalTarget(self.file_path)   



class CreateAllBaseLayers(luigi.WrapperTask):
    bbox = luigi.ListParameter()

    def requires(self):
        yield CreateBaseLayers(bbox=self.bbox)
        yield CreateCid(bbox=self.bbox)
        yield CreateClimid(bbox=self.bbox)