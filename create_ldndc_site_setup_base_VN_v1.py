#
# create_site_Setup_base.py
# =========================
# 
# Classes, functions and data mapping dictionaries used in the 
# specialized site_setup scripts... 
#
# 2011-03-28; christian.werner@senckenberg.de

import math, string, os
import datetime as dt
import xml.dom.minidom as MD
import xml.etree.cElementTree as ET

# ----------- Default data attributes -----------------------

# ---------------------------------------------------------------------
#  (0)  Default data for this dataset
# ---------------------------------------------------------------------
AUTHOR   = 'Christian Werner'
EMAIL    = 'christian.werner@senckenberg.de'
DATE     = str(dt.datetime.now())
DATASET  = 'Phillipines GHG Emissions'
VERSION  = 'v1.0'
SOURCE   = 'BiK-F'

# --------------------------------- D I C T I O N A R I E S --------------------------------------------
# crude mapping from vegtype to humus type
LitTypes = {'FASY': 'MODER', 
            'PIAB': 'MODER', 
            'PISY': 'RAWHUMUS', 
            'QURO': 'MULL', 
            'BEPE': 'MODER', 
            'EUTE': 'MULL'}


# ------------------------------------------- C L A S S E S --------------------------------------------
class BaseXML( object ):
    def __init__(self, startY, endY):
        self.xml     = None             # ET.Element("setups"), define by child
        self.startY  = startY
        self.endY    = endY
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




class SiteParaXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat = str(k['lat'])
        lon = str(k['lon'])

        if k.has_key('id'):
            theId = "%d" % k['id']
        else:
            theId = "0"

        self.xml = ET.Element( "siteparameters", id=theId, lat=lat, lon=lon )
        self.xml.append(self.tags['desc'])


    def addPar(self, DATA, ID=None):
        par = ET.SubElement(self.xml, "par", name=DATA['name'], value=DATA['value'])


class SpeciesParaXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat = str(k['lat'])
        lon = str(k['lon'])

        if k.has_key('id'):
            theId = "%d" % k['id']
        else:
            theId = "0"

        self.xml = ET.Element( "ldndcspeciesparameters" )
        self.xml.append(self.tags['desc'])

        ET.SubElement(self.xml,  "speciesparameters", id=theId, lat=lat, lon=lon )

    def addSpecies(self, mnemonic, name, PARS, group=None, master='none'):

        if master == 'create':
            # add a basic species tag from which the following derive
            species = ET.Element("species", mnemonic=mnemonic, group=group, name=name)
            self.xml.find('./speciesparameters').append(species) 

        if PARS != []:
            if group is not None:
                species = ET.Element("species", mnemonic=mnemonic, group=group, name=name)
            else:
                species = ET.Element("species", mnemonic=mnemonic, name=name)
            
            for parD in PARS:
                name, value, comment = parD
                if comment != None:
                    ET.SubElement( species, "par", name=name, value=str(value), comment=comment)
                else:
                    ET.SubElement( species, "par", name=name, value=str(value) )

            # check if we have a
            if master=='append':
                self.xml.find('./speciesparameters/species').append(species) 
            elif master=='none':
                self.xml.find('./speciesparameters').append(species)
            else:
                print('addSpecies requires keyword master = [create,append,none]; assuming none')
                self.xml.find('./speciesparameters').append(species)
                



class SetupXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat  = str(k['lat'])
        lon  = str(k['lon'])
        try:
            elev = str(k['elev'])
        except KeyError:
            #print 'Using default elevation of 50m'
            elev = "50.0"

        if k.has_key('id'):
            theId = "%d" % k['id']
        else:
            theId = "0"

        if k.has_key('model'):
            theModel = str(k['model'])
        else:
            theModel = "mobile"

        if k.has_key('area'):
            # convert ha -> m-2
            theArea = "%d" % (k['area'] * 10000)
        else:
            theArea = "10000"

        if k.has_key('region'):
            theRegion = k["region"]
        else:
            theRegion = "undefined"

        if k.has_key('province'):
            theProvince = k["province"]
        else:
            theProvince = "undefined"


        self.xml = ET.Element( "ldndcsetup" )
        self.xml.append(self.tags['desc'])

        setup    = ET.SubElement(self.xml, "setup", id=theId, region=theRegion, province=theProvince) #, name="na", model=theModel, active="on")
        
        location = ET.SubElement(setup, "location", elevation=elev, latitude=lat, longitude=lon, \
                                                   slope="0", zone="", aspect="")
        topo     = ET.SubElement(setup, "topology", area=theArea )

        use =      ET.SubElement(setup, "use")

        ET.SubElement(use, "climate", id=str(k['cid']))
        ET.SubElement(use, "airchemistry", id=str(k['acid']))
        ET.SubElement(use, "event", id=str(k['mana']) )
        ET.SubElement(use, "site", id=theId) 
        #ET.SubElement(use, "siteparameters", id=str(k['mana']) )
        ET.SubElement(use, "speciesparameters", id=theId) 
        # hardcoded
        #str(k['mana']) )


        mobile  = ET.SubElement(setup,  "mobile")
        modules = ET.SubElement(mobile, "modulelist")

        mods = [("microclimate:canopyecm"      , "subdaily"),
                ("microclimate:dndc-impl"      , "subdaily"),
                ("watercycle:dndc"             , "subdaily"),
                ("physiology:photofarquhar"    , "subdaily"),
                ("physiology:farquhardndc"     , "subdaily"),
                ("soilchemistry:metrx"         , "subdaily"),

                ("output:soilchemistry:daily"  , ""),
                ("output:soilchemistry:yearly" , ""),
                ("output:watercycle:daily"     , ""),
                ("output:microclimate:daily"   , ""),
                ("output:physiology:daily"     , ""),
                ("output:agmip"                , ""),
                ("output:report:arable:fertilize", ""),
                ("output:report:arable:manure", "")
                ]

        for i, timemode in mods:
            if timemode != "":
                x = ET.SubElement(modules, "module", id=i, timemode=timemode)
            else:
                x = ET.SubElement(modules, "module", id=i)
            if "metrx" in i:
                ET.SubElement(x, "options", algae="yes", iron="yes", paddy="yes", wetland="yes")


class EventXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat = str(k['lat'])
        lon = str(k['lon'])

        if k.has_key('id'):
            theId = str(k['id'])
        else:
            theId = "0"

        self.xml = ET.Element( "event", id=theId, lat=lat, lon=lon )
        self.xml.append(self.tags['desc'])

    def addFire(self, DATA, ID=None):
        timeStr="%d-%d-%d-24/24" % (DATA['date'].year , DATA['date'].month, DATA['date'].day)
        event = ET.SubElement(self.xml, "event", type='fire', time=timeStr)
        fire  = ET.SubElement(event, "fire",  burnedarea=DATA['burnedarea'])

    def addPlant(self, DATA, ID=None):
        ''' this adds a plant/ PFT to the patch '''
        # args:
        # ptype WOOD | GRASS
        # type EUGL
        # name Eucalypt
        # covercrop NONE

        #timeStartStr = "%d-01-01-24/24" % self.startY
        timeStartStr = "01/24"

        event = ET.SubElement(self.xml, "event", type='plant', time=timeStartStr)
        plant  = ET.SubElement(event, "plant", type=DATA['type'], name=DATA['name'])
        ptype = DATA['ptype']
        DATA.pop('ptype', None)
        DATA.pop('type',  None)
        DATA.pop('name',  None)

        e = ET.SubElement(plant, ptype)
        for key, val in DATA.items():
            e.attrib[key] = val

class ProjectXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat  = str(k['lat'])
        lon  = str(k['lon'])

        if k.has_key('id'):
            theId = "None"
        else:
            theId = "0"

        if k.has_key('slim'):
            slim = True
        else:
            slim = False


        if k.has_key('sourceprefix'):
            thePrefix = k["sourceprefix"]
            theSite   = thePrefix.split('/')
            theSite   = theSite[-1] + "_"
        else:
            thePrefix = "parameters"
            theSite   = ""

        thePrefix += "/"

        days = k['days']
        self.xml = ET.Element( "ldndcproject" )

        self.xml.append(self.tags['desc'])

        timeStr = "%d-01-01/24 -> %d-01-01" % (startY, endY+1)
        schedule = ET.SubElement(self.xml, "schedule", time=timeStr)

        input   = ET.SubElement(self.xml, "input" )
        sources = ET.SubElement(input, "sources", sourceprefix="files/")

        #ET.SubElement(sources, "speciesparameters", sourceprefix=thePrefix, source=theSite+"speciesparams.xml")
        #ET.SubElement(sources, "siteparameters", sourceprefix=thePrefix, source=theSite+"siteparams.xml")

        attribs = ET.SubElement(input, "attributes", use="0")
        ET.SubElement(attribs, "airchemistry", endless="yes")

        output  = ET.SubElement(self.xml, "output" )
        sinks   = ET.SubElement(output, "sinks", sinkprefix="%03r-")

        ET.SubElement(sinks, "metrxdaily", sink="metrx-daily.txt", format="txt")
        ET.SubElement(sinks, "agmip", sink="agmip.txt", format="txt")
        #ET.SubElement(sinks, "metrxdailylayer", sink="metrx-daily-layer.txt", format="txt")
        #ET.SubElement(sinks, "metrxsubdaily", sink="metrx-subdaily.txt", format="txt")


    def addFiles(self, DATA, sourceprefix=None, sinkprefix=None, multi=False, sinkcounter=True):
        ''' add file entry to input section (dri1, dri2, site, etc. )'''
        e = self.xml.findall('input/sources')
        e = e[0]
        if sourceprefix != None:
            e.attrib['sourceprefix'] = sourceprefix

        #e.insert(0, ET.Element("siteparameters",    source=DATA['sitepara']))
        e.insert(0, ET.Element("speciesparameters", source=DATA['speciespara']))
        e.insert(0, ET.Comment("parameter files"))

        e.insert(0, ET.Element("airchemistry", source=DATA['airchem'], format="txt"))
        e.insert(0, ET.Element("event",        source=DATA['event']))
        e.insert(0, ET.Element("climate",      source=DATA['clim'],    format="txt"))
        e.insert(0, ET.Element("setup",        source=DATA['setup']))
        e.insert(0, ET.Element("site",         source=DATA['site']))
        e.insert(0, ET.Comment("control files"))



        e = self.xml.findall('output/sinks')
        e = e[0]
        if sinkprefix != None:
            if sinkcounter==True:
                e.attrib['sinkprefix'] = sinkprefix + "%03r_"
            else:
                e.attrib['sinkprefix'] = sinkprefix


class SiteXML( BaseXML ):
    def __init__(self, startY, endY, **k):
        BaseXML.__init__(self, startY, endY)
        lat   = str(k['lat'])
        lon   = str(k['lon'])
        if k.has_key('id'):
            theId = "%d" % k['id']
        else:
            theId = "0"
        self.xml = ET.Element( "site", id=theId, lat=lat, lon=lon )
        self.xml.append(self.tags['desc'])
        
        # gernal tags
        general = ET.SubElement(self.xml, "general")
        #ET.SubElement(general, "location", lat  =lat,      lon   = lon, \
        #                                   elev ='0',      zone  ='-1')
        #ET.SubElement(general, "terrain",  slope='0.0',    aspect='-99.99')
        #ET.SubElement(general, "climate",  tavg ='-99.99', tamp  ='-99.99', \
        #                                   psum ='-99.99', riref ='-99.99', \
        #                                   cloud='-99.99', wavg  ='-99.99')
        #ET.SubElement(general, "misc",     watertable = '')

        # soil tags
        soil       = ET.SubElement(self.xml, "soil")
        ET.SubElement(soil, "general", usehistory='arable', soil        ='NONE', \
                                       humus     ='NONE', lheight='0.0',  \
                                       corg5     ='-99.99', corg30      ='-99.99')
        layers = ET.SubElement(soil, "layers")

    def addSoilLayer(self, DATA, ID=None, litter=False, accuracy={}):
        ''' this adds a soil layer to the given site (to current if no ID given)'''
        # only calculate hydr. properties if we have a mineral soil layer added
        if litter == False:
            DATA['wcmax'],  DATA['wcmin'] = calcHydaulicProperties( DATA )
        
        soilLayer=ET.Element("layer", depth='-99.99', split ='1', ph  ='-99.99', \
                                      scel ='-99.99', bd    ='-99.99', sks ='-99.99', \
                                      norg ='-99.99', corg  ='-99.99', clay='-99.99', \
                                      wcmax='-99.99', wcmin ='-99.99', sand='-99.99', \
                                      silt ='-99.99', iron  ='-99.99')
        keys = DATA.keys()
        for k in keys:
            digits=2
            if k in accuracy.keys(): 
                digits = accuracy[k]
                if digits == 0:
                    # int
                    soilLayer.attrib[k] = str(int(round(DATA[k], digits)))
                else:
                    soilLayer.attrib[k] = str(round(DATA[k], digits))

        self.xml.find('./soil/layers').append(soilLayer) 

# ------------------------------------------- F U N C T I O N S --------------------------------------------

def calcLitter(litterMass, litname): # mass in t C ha-1
    if litname == 'MULL':
        density = 0.2;  accumulationFactor = 1.5
    elif litname == 'MODER':
        density = 0.25; accumulationFactor = 2.5
    else:
        density = 0.3;  accumulationFactor = 3.5
        
    # explanation:
    #   (tCha-1) > x2 > tBMha-1 > x0.1 > kgBMm-2 > x0.1 > gBMcm-2 > /density > height_cm > *10 > height_mm  
    # littermass (t C ha-1) * 2 (BM conv) * 0.1 * 0.1 / density * 10
    depth = ((litterMass * accumulationFactor * 2 * 0.1 * 0.1) / density) * 10.0
    numberOfLayers = math.floor( depth / 20.0 )
    layerHeight = depth
    if numberOfLayers != 0:
        layerHeight = 20.0 + ((depth % 20.0) / numberOfLayers) 
    return (density, depth, layerHeight)

def calcHeight(TK, N): 
    if TK == -9999:
        TKmm = -9999
    else:
        TKmm = TK * 10 #change the unit from cm to mm

    if TKmm > -9999:
        numberOfLayers = math.floor( TKmm / N )
        if numberOfLayers != 0:
            layerHeight    = N + ((TKmm % N) / numberOfLayers) 
        else:
            layerHeight = TKmm
    else:
        layerHeight = -9999
    return (TKmm, layerHeight)

def calcHydaulicProperties(D):
    ''' Calc hydraulic properties based on et al. (1996) '''
    # shape parameters Woesten et al. (1999) Geoderma
    #
    # OM      (% organic matter)
    # D       (bulk denisty)
    # topsoil 1, subsoil 0
    # C, S,   (clay, silt in %)
    #
    #ThetaS = 0.7919 + 0.001691 * C - 0.29619 * D - 0.000001491 * S*S + 0.0000821 * OM * OM + 0.02427 * C**-1 + 0.01113 * S**-1 + \
    #         0.01472 * math.ln( S ) - 0.0000733 * OM * C - 0.000619 * D * C - 0.001183 * D * OM - 0.0001664 * topsoil * S
    
    # ad-hoc AG Boden
    #
    # Sand, Clay [%], BD [g cm-3], corg [%]

    corg = D['corg'] * 100
    clay = D['clay'] * 100
    sand = D['sand'] * 100
    bd   = float(D['bd'])   

    ThetaR = 0.015 + 0.005 * clay + 0.014 * corg
    ThetaS = 0.81 - 0.283 * bd + 0.001 * clay

    logAlpha = -2.486 + 0.025 * sand - 0.351 * corg - 2.617   * bd - 0.023 * clay
    logN     =  0.053 - 0.009 * sand - 0.013 * clay + 0.00015 * sand**2

    try:
        ALPHA = math.e**logAlpha
    except:
        print( D )
    vGn = math.e**logN
    vGm = 1.0 # (1.0 - (1.0/ vGn)) disabled as we do not use texture classes but real fractions

    FLDcap = ThetaR + (ThetaS- ThetaR) / math.pow( ( 1.0 + math.pow ( ALPHA * 100.0, vGn ) ), vGm )  
    WILTpt = ThetaR + (ThetaS- ThetaR) / math.pow( ( 1.0 + math.pow ( ALPHA * 15800.0, vGn ) ), vGm )
    return FLDcap * 1000, WILTpt * 1000


def calcHydaulicPropertiesSavanna(D):
    ''' Calc hydraulic properties based on Medarado and Lima, Geoderma Regional (2014)'''
    # Sand, Clay [%], BD [g cm-3], corg [%]

    SOC = D['corg'] * 100
    CLAY = D['clay'] * 100
    SILT = D['silt'] * 100
    SAND = D['sand'] * 100
    BD   = float(D['bd'])


    #x1=bulkdensity
    #x2=clay(%)
    #x3 = total sand (%)
    #x4 = silt(%)
    #x5 = organic matter (%)

    # K
    a5_1 = 0.88626397; b5_1 = 0.00786574

    # theta_r
    a1_2 = 0.11268382; b1_2 = -2.31854838
    a3_2 = 0.17891998; b3_2 = -0.41220364
    a5_2 = 0.01329004; b5_2 = 1.02011236

    # n
    a1_4 = 0.71330436;  b1_4 = -1.28657786
    a2_4 = 3.21317122;  b2_4 = -0.29289780
    a4_4 = -0.23122583; b4_4 = 0.12757866
    a5_4 = -0.00511163; b5_4 = 1.68342787

    # alpha
    a1_3 = 6.32532102;  b1_3 = -0.02318223
    a2_3 = -7.49868149; b2_3 = -0.03241971
    a4_3 = -7.26003177; b4_3 = 0.00905719
    a5_3 = 8.51527042;  b5_3 = 0.00322043

    # K = sum(a*x^b)
    K = a5_1 * pow( SOC, b5_1 )

    soilparticledensity = 2.65

    # soil water porosity
    ThetaP = 1.0/BD - 1.0/soilparticledensity
    ThetaS = K * ThetaP

    ThetaR = a1_2 * pow(BD, b1_2) + a3_2 * pow(SAND, b3_2) + a5_2 * pow(SOC, b5_2)
    vGn    = a1_4 * pow(BD, b1_4) + a2_4 * pow(CLAY, b2_4) + a4_4 * pow(SILT, b4_4) + a5_4 * pow(SOC, b5_4)
    vGm    = 1.0  - (1.0 / vGn);
    ALPHA  = a1_3 * pow(BD, b1_3) + a2_3 * pow(CLAY, b2_3) + a4_3 * pow(SILT, b4_3) + a5_3 * pow(SOC, b5_3)

    # formula in kPa !!!
    # FLDcap for 10 kPa
    # WILTpt for 1580 kPa
    # used to be hPa !!!

    FLDcap = ThetaR + (ThetaS- ThetaR) / pow( ( 1.0 + pow ( ALPHA * 10.0, vGn ) ), vGm )
    WILTpt = ThetaR + (ThetaS- ThetaR) / pow( ( 1.0 + pow ( ALPHA * 1580.0, vGn ) ), vGm )
    return FLDcap * 1000, WILTpt * 1000

