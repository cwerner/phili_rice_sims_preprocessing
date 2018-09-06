import xml.dom.minidom as MD
import xml.etree.cElementTree as ET
import glob
import os
import sys
import datetime as dt
import xarray as xr

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
AUTHOR   = 'Christian Werner'
EMAIL    = 'christian.werner@senckenberg.de'
DATE     = str(dt.datetime.now())
DATASET  = 'n.a.'
VERSION  = 'v0.1'
SOURCE   = 'BIK-F'


class BaseXML( object ):
    def __init__(self):
        self.xml     = None             # ET.Element("setups"), define by child
        self.tags = {}
        # --- dict of initial xml tags
        desc = ET.Element("description")
        e = ET.SubElement(desc, "author");  e.text = AUTHOR
        e = ET.SubElement(desc, "email");   e.text = EMAIL        
        e = ET.SubElement(desc, "date");    e.text = DATE
        e = ET.SubElement(desc, "dataset"); e.text = DATASET
        e = ET.SubElement(desc, "version"); e.text = VERSION
        e = ET.SubElement(desc, "source");  e.text = SOURCE
        self.tags['desc'] = desc

    def write(self, ID = None, filename = 'all.xml'):
        strOut = MD.parseString(ET.tostring(self.xml)).toprettyxml()
        # fix special characters
        sc = {'&gt;': '>', '&lt;': '<'}
        for key, val in sc.items():
            strOut = string.replace(strOut, key, val)
        open(filename, 'w').write(strOut)



class SetupXML( BaseXML ):
    def __init__(self, **k):
        BaseXML.__init__(self)
        lat  = str(round(k['lat'], 6))
        lon  = str(round(k['lon'], 6))
        x=str(k['x'])
        y=str(k['y'])
        try:
            elev = str(k['elev'])
        except KeyError:
            #print 'Using default elevation of 50m'
            elev = "50.0"

        if 'id' in k:
            theId = "%d" % k['id']
        else:
            theId = "0"

        if 'model' in k:
            theModel = str(k['model'])
        else:
            theModel = "mobile"

        self.xml = ET.Element( "ldndcsetup" )
        self.xml.append(self.tags['desc'])

        setup    = ET.SubElement(self.xml, "setup", id=theId) #, name="na", model=theModel, active="on")
        
        location = ET.SubElement(setup, "location", elevation=elev, latitude=lat, longitude=lon, \
                                                   slope="0", zone="", aspect="")
        topo     = ET.SubElement(setup, "topology", area='10000', x=x, y=y,z=elev )


        models = ET.SubElement(setup, "models")
        ET.SubElement(models, "model", id="RiceFarmer")
        ET.SubElement(models, "model", id="_MoBiLE")
        ET.SubElement(models, "model", id="IPCC")

        ricefarmer = ET.SubElement(setup, "RiceFarmer", file="%I/regional/PH_arable/PH_arable_ext_awd_hm.rice", id=str(int(k['mana'])))
        
        use =      ET.SubElement(setup, "use")

        ET.SubElement(use, "climate", id=str(k['cid']))
        ET.SubElement(use, "airchemistry", id=str(k['acid']))
        ET.SubElement(use, "site", id=theId) 

        mobile  = ET.SubElement(setup,  "mobile")
        modules = ET.SubElement(mobile, "modulelist")

        mods = [("microclimate:canopyecm"      , "subdaily"),
                ("microclimate:dndc-impl"      , "subdaily"),
                ("watercycle:dndc"             , "subdaily"),
                ("airchemistry:depositiondndc" , "subdaily"),
                ("physiology:photofarquhar"    , "subdaily"),
                ("physiology:plamox"           , "subdaily"),
                ("soilchemistry:metrx"         , "subdaily"),

                ("output:ecosystem:daily"      , ""),
                ("output:ecosystem:yearly"     , ""), 
                ("output:microclimate:daily"   , ""),
                ("output:watercycle:daily"     , ""),
                ("output:watercycle:yearly"    , ""),
                ("output:physiology:daily"     , ""),
                ("output:soilchemistry:daily"  , ""),
                ("output:soilchemistry:yearly" , ""),
                ("output:report:arable"        , "subdaily")]

        for i, timemode in mods:
            if timemode != "":
                x = ET.SubElement(modules, "module", id=i, timemode=timemode)
            else:
                x = ET.SubElement(modules, "module", id=i)
            if "metrx" in i:
                ET.SubElement(x, "options", algae="yes", wetland="yes")


def _decompose(cid, multiplier=1000):

    lMulti = len(str(multiplier)) - 1
    j = int(str(cid)[:-lMulti])
    i = int(str(cid)[-lMulti:])

    return (j, i)



if __name__ == "__main__":

    ds = xr.open_dataset('data/PH_MISC.nc')

    lats = ds.lat.values #sorted(ds.lat.values, reverse=True) #sorted(ds.lat.values, reverse=True)    #FLIP WHEN DAVID FIXED THE MANA IDS
    lons = ds.lon.values

    #print "WARNING !!! cids are wrong and range from lower left to upper right !!!"

    # create a list of cids from site file
    site_ids = []
    tree = ET.ElementTree(file=sys.argv[1])
    for elem in tree.iterfind('site'):
        cid = int(elem.attrib['id'])
        site_ids.append(cid)

    setups = []

    for j in range(len(ds['simmask'])):
        for i in range(len(ds['simmask'][0])):
            if ds['simmask'][j,i] == 1:
                cid = int(ds['cid'][j,i])

                if cid in site_ids:

                    elev = ds['elevation'][j,i].values
                    manaid = ds['province'][j,i].values

                    setup = SetupXML( lat=lats[j], lon=lons[i], elev=elev, id=cid, cid=cid, acid=0, mana=manaid, y=j, x=i) 
                    setups.append(setup)                    



    print(f"Total Nuber of setup items: {len(setups)}...")

    # merge setup chunks into common setup file
    xml = ET.Element( "ldndcsetup" )
    for scnt, setup in enumerate(setups):
        if scnt == 0:
            xml.append(setup.xml.find("description"))
        xml.append(setup.xml.find("setup"))

    strOut = MD.parseString(ET.tostring(xml)).toprettyxml()
    open('setups_HR.xml', 'w' ).write( strOut )

