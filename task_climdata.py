# Management Task
# ===============

import luigi
import gzip
import os
import pandas as pd
import xarray as xr
from io import StringIO

from general import RasterizeShapefile 

RES = 0.5

#./DLSC.py -r HR --bbox="[116,4,127,22]" --country=102
#
#./DLSC.py -r HR --bbox="[116,4,127,22]" --country=102

class ExtractClimData(luigi.Task):
    bbox = luigi.ListParameter()
    syear  = luigi.IntParameter(default=None)
    eyear  = luigi.IntParameter(default=None)

    path = luigi.Parameter(default='tmp')    

    def run(self):
        lats_, lons_ = bbox_to_latlon(self.bbox, RES)

        # create cin variable
        DIM0 = len(lats_)
        DIM1 = len(lons_)
        data = np.flipud( (np.arange(DIM0*DIM1) + 1).reshape((DIM0, DIM1)) )
        da = xr.DataArray(data, coords=(lats_, lons_), dims=('lat', 'lon'), name='cin')
        da.attrs['units'] = 'id'
        
        FILE = os.path.join(self.path, self.fname)
        da.to_dataset().to_netcdf(FILE, format='NETCDF4_CLASSIC')

    def output(self):
        FILE = os.path.join(self.path, self.fname)
        return luigi.LocalTarget( FILE )
    
    def requires(self):
        yield RasterizeShapefile(shp_file=FILE, attr="NONE", bbox=self.bbox, res=RES, name="region", file_path="tmp/manamask.nc")
        yield RasterizeShapefile(shp_file=FILE, attr="MANA_ID", bbox=self.bbox, res=RES, name="province", file_path="tmp/manaids.nc")



def get_cell_info(data_nc):
    """Extract site info from base netcdf file
    """
    data_nc_stacked = data_nc.stack(cell=['lat','lon'])
    return data_nc_stacked.where(data_nc_stacked.region == 1, drop=True)


class DownsampleClimateData(luigi.Task):
    """Downsample climate files (0.5x0.5deg) based on high-res topography using lapse rate and topodiff
    
    requires: <VN,PH>_MISC.nc file
    """

    ADIABATICLR = 0.65    # 0.5 degrees by 100m; lapse rate

    data_nc = luigi.Parameter()

    out_file = luigi.Parameter(default="climate.txt")
    climate_dbpath = luigi.Parameter(default='data/climate/db')

    syear = luigi.IntParameter(default=None)
    eyear = luigi.IntParameter(default=None) 

    path = luigi.Parameter(default='tmp')    

    file_path = 'data/climate/climate_HR.txt' 

    def run(self):
        YEARS = None
        if self.syear is not None and self.eyear is not None:
            YEARS = range(self.syear, self.eyear+1)

        data = get_cell_info(xr.open_dataset(self.data_nc))


        for cnt, c in enumerate(data.cell):
            data_cell = data.sel(cell=c)
            lat, lon = data_cell.cell.item()
            
            sid = data_cell.cid.values
            climate_id = data_cell.climid.values
            #lat = data_cell.lat.values
            #lon = data_cell.lon.values
            elevation = data_cell.elevation.values
            topodiff = data_cell.topodiff.values

            if cnt == 0:
                fout = gzip.open( self.output().path + '.gz', 'wt')

            climate_dbfile = os.path.join(self.climate_dbpath, "climate_%08d.txt.gz" % climate_id)
            climate_lines = gzip.open( climate_dbfile, 'rt').readlines()
            
            # parse climate file and modify temp(s)
            if cnt == 0:
                # write global section (only once)
                globalHeader = f'%global\n\t\ttime     = {self.syear}-1-1/1\n%cuefile = "cuefile.txt"\n\n'

                for lcnt, line in enumerate(climate_lines):
                    if "climate\n" in line:
                        break

                fout.write( "".join( climate_lines[0:lcnt-1] ) )
                fout.write( f'%cuefile = "cuefile.txt\n\n')

            # remove those lines
            climate_lines = climate_lines[lcnt-1:]

            # write all lines until we see the %data flag
            for lcnt, line in enumerate(climate_lines):
                if "id" in line:
                    climate_lines[lcnt] = f"        id = {sid}\n"
                if "elevation" in line:
                    climate_lines[lcnt] = f"        elevation = {elevation}\n"
                if "latitude" in line:
                    climate_lines[lcnt] = f"        latitude  = {lat}\n"
                if "longitude" in line:
                    climate_lines[lcnt] = f"        longitude = {lon}\n"
                if "data\n" in line:
                    break

            # fix elevation line
            #for i, l in enumerate(climate_lines[0:lcnt]):
            #    if "elevation" in l:
            #        climate_lines[i] = "        elevation = %d" % 

            COMMENT = (f"    # Data height-adjusted (TOPODIFF)\n" +
                       f"    #   original climid (0.5deg res):     {climate_id}\n" +
                       f"    #   new id == cid (native res):       {sid}\n" +
                       f"    #   elevation difference to orig [m]: {topodiff}\n" +
                       f"    #   lapse rate corr. [degC/ 100m]:    {self.ADIABATICLR}\n" +
                       f"    #\n" +
                       f"    # WARNING: tavg & tamp are not corrected (for now)\n")

            sectionHeader = climate_lines[0: lcnt+2]
            sectionHeader[-3:-3] = COMMENT.splitlines(True)

            headerLength = len(sectionHeader)-2

            # write all header data to file
            fout.write( "".join( sectionHeader ) )
            # delete header lines, but not column header line
            climate_lines = climate_lines[lcnt+1:]

            df = pd.read_csv( StringIO("".join( climate_lines )), sep='\t', header=0 )
            df.rename(columns={'*': "date", "*.1": "year", '*.2': "doy"}, inplace=True)
            
            # limit years
            if YEARS:
                df = df[(df.year >= YEARS[0]) & (df.year <= YEARS[-1])]

            temp_adjust = self.ADIABATICLR * topodiff * 0.01 

            # adjust temp columns
            if df['tavg'].mean(axis=0) > -99:
                df['tavg'] = df['tavg'] - temp_adjust
                df['tavg'] = df.tavg.round(1)
            
            df['tmin'] = df['tmin'] - temp_adjust
            df['tmax'] = df['tmax'] - temp_adjust
            df['tmin'] = df.tmin.round(1)
            df['tmax'] = df.tmax.round(1)

            out = df.to_csv(None, sep="\t", index=False, header=False) + '\n'
            fout.write( out )

        fout.close()
    
    def output(self):
        return luigi.LocalTarget(self.file_path)  

