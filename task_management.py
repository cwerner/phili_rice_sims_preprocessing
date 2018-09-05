# Management Task
# ===============

import luigi

from general import RasterizeShapefile 

RES = 0.0833333333
FILE = "original_data/PH_rice-agriculture.shp"

class ProcessManagement(luigi.WrapperTask):
    bbox     = luigi.ListParameter()
    
    def requires(self):
        yield RasterizeShapefile(shp_file=FILE, attr="NONE", bbox=self.bbox, res=RES, name="region", file_path="tmp/manamask.nc")
        yield RasterizeShapefile(shp_file=FILE, attr="MANA_ID", bbox=self.bbox, res=RES, name="province", file_path="tmp/manaids.nc")

