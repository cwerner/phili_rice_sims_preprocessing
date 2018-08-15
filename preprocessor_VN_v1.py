# preprocessor.py
#
# descr: preprocessor for aussi regional simulations
#
# Christian Werner (BiK-F)
# christian.werner@senckenberg.de
#
# 2013-06-10 started (only climate file for now)
# 2013-06-20 extend functionality
#             - also parse firedata netcdf file
#             - use scoop to parallelize program if desired
# 2015-02-20 v0.5
# 2015-03-05 v0.6 major rewrite (NEW climdb dump)
# 2015-03-05 v0.7 new hydrology, other fixes

# enable if coords should be created in parallel
#from scoop import futures, shared

import sys, string, os, datetime, errno, gzip, math
import numpy as np
import datetime as dt
from netCDF4 import Dataset
import operator
from optparse import OptionParser
import itertools
from create_ldndc_site_setup_base_VN_v1 import SiteXML, SetupXML, EventXML, ProjectXML, SiteParaXML, SpeciesParaXML
import xml.dom.minidom as MD
import xml.etree.cElementTree as ET
import progressbar


import pandas as pd
from StringIO import StringIO

# ENUMS
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)

Mode = enum('FIXED', 'NATT', 'DYN')             # we have 3 veg modes: FIXED, NATT, DYN

# CONSTS (columns and precision)

# climate
CLIMCOLS = [('years',0), ('jdays',0), \
            ('tas',1), ('tmin',1), ('tmax',1), \
            ('dswrf',1), ('prcp',2), ('rhum',1), ('wind',2)]

# airchem
AIRCHEMCOLS = [('years',0), ('jdays',0), \
               ('n_total',8), ('n_dry',  8), ('n_wet',8), \
               ('noy_wet',8), ('nhx_wet',8), \
               ('noy_dry',8), ('nhx_dry',8)]

# header info airchem file: * * n ndry nwet no3 nh4 no3dry nh4dry co2


# soil
nmap = { "TOTC": ('corg', 0.001, 5),
         "TOTN": ('norg', 0.001, 6),
         "PHAQ": ('ph', 1, 2),
         "BULK": ('bd', 1, 2),
         "CFRAG": ('scel', 0.01, 2),
         "SDTO":  ('sand', 0.01,  2),
         "STPC":  ('silt', 0.01,  2),
         "CLPC":  ('clay', 0.01,  2),
         "TopDep": ('topd', 1, 0),
         "BotDep": ('botd', 1, 0)}

cmap = dict((x[0], x[2]) for x in nmap.values())
cmap['depth']  = 0
cmap['split']  = 0
cmap['wcmin']  = 1
cmap['wcmax']  = 1
cmap['iron']   = 5

# Helper functions -------------------------------------------------------
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = MD.parseString(rough_string)
    str1 = reparsed.toprettyxml(indent="  ")
    str2 = []
    ss = str1.split('\n')
    for s in ss:
        x = "".join(s.split())
        if x != "":
            str2.append(s)
    return '\n'.join( str2 ) + '\n'



# Airchemistry data functions --------------------------------------------------------------------

# Header info (A)
def airchem_headerA(fname, sdate):
    ''' default header for ldndc daily airchemistry file (1)'''
    header = "\n".join(['%global',
                        '\ttime = "%s/1"\n' % sdate]) + '\n'
    return header

# Header info (B)
def airchem_headerB(fname, cdata):
    ''' default header for ldndc daily airchemistry file (2)'''
    header = "\n".join(['',
                        '%airchemistry',
                        '\tname = "%s"' % fname,
                        '\tarchive = "%s"' % cdata['climarchive'],
                        '\tid   = "%d"\n' % cdata['cid'],                        
                        '%data',
                        '* * n ndry nwet no3 nh4 no3dry nh4dry co2']) + '\n'
    return header

# Data formatting and sorting
def arrangeDataAirchem( data, cols, co2=370):
    """ return data dict as list of outlines sorted by cols """
    rdata = []

    # hack to bring nh4 / no3 wet deposition concentration to corret units
    data["noy_wet"] *= 1000.0
    data["nhx_wet"] *= 1000.0

    days = len(data["years"])
    for d in range(days):
        # first two cols are int, rest is rounded according to digits
        a = [ data["years"][d], data["jdays"][d] ] + [ round(data[c][d], di) for c, di in cols[2:] ]
        # use int value of co2 or the appropriate value for this year
        if isinstance( co2, ( int, long ) ):
            a.append( round(co2, 1) )
        else:
            a.append( round(co2[data["years"][d]],1) )
        rdata.append( '\t'.join( map(str, a) ) )
    return rdata



# Data format functions --------------------------------------------------------------------

def translateDataFormat( d ):
    """ translate data from nc soil file to NEU naming and units """ 
    data = []
    size = len(d[d.keys()[0]])
    
    ks = nmap.keys()
    for l in range(size):
        od = {}
        for k in ks:
            name, conv, ignore = nmap[k]
            od[ name ] = d[k][l] * conv
        data.append( od )
    return data



# Classes --------------------------------------------------------------------------------

class NCFile( object ):
    def __init__(self, fname, id_map=None, diskless=False):
        self.fname = fname
        self.nc    = Dataset(fname, mode="r")

        # ignore those variables
        i_vars = ['z', 'lev', 'layers', 'record', 'rec', 'slm', 'crs', 'lat', 'lon', 'time', 'longitude', 'latitude']

        try:
            self.lons  = list(self.nc.variables['lon'])
            self.lats  = list(self.nc.variables['lat'])
        except KeyError:
            self.lons  = list(self.nc.variables['longitude'])
            self.lats  = list(self.nc.variables['latitude'])

        # list of data variables
        self.vars  = [x for x in self.nc.variables.keys() if x not in i_vars]
        #self.soil  = False
        #self.clim  = False
        #self.fire  = False

        self.id_map   = id_map

        # store entire array in memory
        self.diskless = diskless

        self.mem = {}
        for v in self.vars:
            self.mem[v] = self.nc.variables[v][:]

    def _closest(self, l, v):
        return min(enumerate(l), key=lambda x: abs(x[1]-v))

    def _get(self, lat, lon):

        if self.id_map == None:
            # translate to closests grid cell idx or gridded lat/ lon spacing
            idx, glon = self._closest(self.lons, lon)
            jdx, glat = self._closest(self.lats, lat)

            # check latitude orientation
            if self.lats[0] > self.lats[-1]:
                print("WRONG LATITUDE ORIENTATION in self.nc.name")
                exit(1)
            theId = (len(self.lons) * jdx) + idx # calc the continous id based on idx, idy
        else:
            theId, idx, jdx = self.id_map[(lat, lon)]

        data = {}
        for var in self.vars:
            dims = len( self.nc.variables[var].dimensions )
            if dims > 3:
                # this is to catch z level variables in clim data 
                if self.diskless:
                    d = self.mem[var][:,0,jdx,idx]
                else:
                    d = self.nc.variables[var][:,0,jdx,idx]

            else:
                try:
                    if dims == 3:
                        if self.diskless:
                            d = self.mem[var][:,jdx,idx]
                        else:
                            d = self.nc.variables[var][:,jdx,idx]
                    elif dims == 2:
                        if self.diskless:
                            d = self.mem[var][jdx,idx]
                        else:
                            d = self.nc.variables[var][jdx,idx]
                except IndexError:
                    if self.diskless:
                        print( var, idx, jdx, self.nc.mem[var].shape )
                    else:
                        print( var, idx, jdx, self.nc.variables[var].shape )
            data[var] = d
        return (theId, data)


    def _getXY(self, x, y):

        idx = x
        jdx = y

        theId = (len(self.lons) * jdx) + idx # calc the continous id based on idx, idy

        data = {}
        for var in self.vars:
            dims = len( self.nc.variables[var].dimensions )

            if dims > 3:
                # this is to catch z level variables in clim data 
                if self.diskless:
                    d = self.mem[var][:,0,jdx,idx]
                else:
                    d = self.nc.variables[var][:,0,jdx,idx]

            else:
                try:
                    if dims == 3:
                        if self.diskless:
                            d = self.mem[var][:,jdx,idx]
                        else:
                            d = self.nc.variables[var][:,jdx,idx]
                    elif dims == 2:
                        if self.diskless:
                            d = self.mem[var][jdx,idx]
                        else:
                            d = self.nc.variables[var][jdx,idx]
                except IndexError:
                    if self.diskless:
                        print( var, idx, jdx, self.nc.mem[var].shape )
                    else:
                        print( var, idx, jdx, self.nc.variables[var].shape )
            data[var] = d
        return (theId, data)


    def create_id_map( self, coords ):
        D = {}
        if self.lats[0] > self.lats[-1]:
            print( "WRONG LATITUDE ORIENTATION in self.nc.name" )
            exit()

        for c in coords:
            lat, lon = c
            idx, glon = self._closest(self.lons, lon)
            jdx, glat = self._closest(self.lats, lat)

            theId = (len(self.lons) * jdx) + idx # calc the continous id based on idx, idy

            # data structure: cellid, lon_position, lat_position
            D[c] = (theId, idx , jdx)

        return D

    def close(self): 
        self.nc.close()



class NDepFile( NCFile ):
    def __init__(self, fname, year, id_map=None):
        super(NDepFile, self).__init__( fname, id_map )
        self.year = year
        #self.clim = True
        self.exists = set()

    def get(self, lat, lon):
        theId, data = super(NDepFile, self)._get(lat, lon)
        vars = data.keys()
        data['years'] = np.array([self.year]*len(list( data[vars[0]] )))
        data['jdays'] = range(1, len( data[vars[0]] )+1)

        return (theId, data)



class SoilFile( NCFile ):
    def __init__(self, fname, id_map=None):
        super(SoilFile, self).__init__( fname, id_map )
        self.soil = True

    def get(self, lat, lon):
        theId, data = super(SoilFile, self)._get(lat, lon)
        vars = data.keys()
        return (theId, data)

    def getXY(self, x, y):
        theId, data = super(SoilFile, self)._getXY(x, y)
        vars = data.keys()
        return (theId, data)



    def validids(self, latlimit=None):
        lats = self.nc.variables['lat'][:]
        lons = self.nc.variables['lon'][:]
        subset_lats = []

        if latlimit != None:
            subset_lats = list(lats[ (lats >= latlimit[0]) & (lats < latlimit[1]) ])
        else:
            subset_lats = lats
        
        # use slm, also remove
        mask = self.nc.variables['slm'][:] * np.where(self.nc.variables['TopDep'][0,:] < 0.0, 0, 1)
        data = []

        # we go from north to south -> reverse 
        for j, lat in zip(range(len(lats)), reversed(list(lats))):
            for i, lon in enumerate(lons):
                if (lat in subset_lats) and (mask[ len(list(lats))-j-1,i] == 1):
                    data.append( (lat, lon, j, i) )

        return data

class RefFile( NCFile ):
    def __init__(self, fname, id_map=None):
        super(RefFile, self).__init__( fname, id_map )
        self.ref = True

    def get(self, lat, lon):
        theId, data = super(RefFile, self)._get(lat, lon)
        vars = data.keys()
        return (theId, data)

    def getXY(self, x, y):
        theId, data = super(RefFile, self)._getXY(x, y)
        vars = data.keys()
        return (theId, data)


    def validids(self, latlimit=None):
        lats = self.nc.variables['latitude'][:]
        lons = self.nc.variables['longitude'][:]
        subset_lats = []

        if latlimit != None:
            subset_lats = list(lats[ (lats >= latlimit[0]) & (lats < latlimit[1]) ])
        else:
            subset_lats = lats
        
        # use slm, also remove
        mask = self.nc.variables['mask'][:]
        data = []


        abort = 1000
        print( "WARNING: debug mode, we only use 1000 coords max !!!" )

        # we go from north to south -> reverse 
        cnt = 1
        for j, lat in zip(range(len(lats)), reversed(list(lats))):
            for i, lon in enumerate(lons):
                if (lat in subset_lats) and (mask[ len(list(lats))-j-1,i] == 1):
                    data.append( (lat, lon, j, i) )

                    cnt += 1
                    if cnt == abort:
                        return data

        return data


def main(sids, clids, coords, idxs, options, OUTPATH):

    # note: internal ids will follow the domain of the netcdf file !

    global soil_file
    global fire_file
    global ref_file

    global climate_dbpath


    # Date (start & end)
    ystart = [int(x) for x in options.years.split("-")][0]
    yend   = [int(x) for x in options.years.split("-")][1]

    sdate_dt = dt.datetime(ystart, 1, 1)
    edate_dt = dt.datetime(yend+1, 1, 1)
    sdate = "%d-%d-%d" % (sdate_dt.year, sdate_dt.month, sdate_dt.day)
    edate = "%d-%d-%d" % (edate_dt.year, edate_dt.month, edate_dt.day)
    days  = (edate_dt - sdate_dt).days

    mana  = options.mana


    # (A) create site,setup file
    # ============================

    print( "   building site/ setup file" )
    
    nc       = SoilFile( soil_file )
    nc_admin = RefFile( ref_file ) 
    id_map_other = nc.create_id_map( coords )
    

    sites  = []  # holds the xml chunks for site file
    setups = []  # holds the xml chunks for setup file

    # gdds
    speciesparas = []  # holds the xml chunks for speciesparameter file

    addFlag = {} #dict key:coord val: bool, indicates if a site has no valid soil info


    print( "   >>> now building r-r-v rotations !!!" )

    # iterate over all sites
    for sid in sids:

        coord = (D_info[sid]['lat'], D_info[sid]['lon'])

        # lid == local id (ignore?)
        lid, data  = nc.getXY( D_info[sid]['x'], D_info[sid]['y'] )

        _, data_admin = nc_admin.getXY( D_info[sid]['x'], D_info[sid]['y'] )


        # if we do not have any agriculture in this cell skip and do not produce datafiles
        # DISABLED FOR PHIL
        #if data_admin["agrifrac"] == 0.0 or data_admin['rice_fr'] == 0.0:
        #    if sid == 3995149:
        #        print 'SKIPPING'
        #
        #            continue

        # translate
        rrot = int(data_admin['rice_rot'])

        # ids:
        # 100 (r-r-r)
        # 200 (r-r)
        # 300 (r-r-v)
        # 400 (r)
        # 500 (r-v) 


        if 'riceonly' in options.mana:

            if rrot == 1:
                rrot_trans = 400    # single rot rice (rest fallow)
            elif rrot == 2:
                rrot_trans = 200    # double rot rice
            elif rrot == 3:
                rrot_trans = 100    # triple rice
            else:
                print("rotation ZERO !!!", sid, Dprovince[ int(data_admin['provinceid']) ])
                #rrot_trans = 200
                exit(-1)

        elif 'mixed' in options.mana:
            # second set
            if rrot == 1:
                rrot_trans = 500
            elif rrot == 2:
                rrot_trans = 300
            elif rrot == 3:
                rrot_trans = 100
            else:
                print( "rotation ZERO !!!", sid, Dprovince[ int(data_admin['provinceid']) ])
                #rrot_trans = 200
                exit(-1)

        else:
            print( "rotation selection unknown" )
            exit(1)

        manaid = str(int(data_admin['provinceid']) * 1000 + rrot_trans)

        D_info[sid]['manaid'] = manaid #str(int(data_admin['provinceid']))
        D_info[sid]['area']   = data_admin['agrifrac']*0.01 * D_info[sid]['area']
        D_info[sid]['region']   = Dregion[ int(data_admin['regionid']) ]
        D_info[sid]['province'] = Dprovince[ int(data_admin['provinceid']) ]

        data2 = translateDataFormat( data )

        # create site & setup xml file instances
        site     = SiteXML(  ystart, yend, lat=D_info[sid]['lat'], lon=D_info[sid]['lon'], id=sid )


        if options.tempadjust:
            cid = sid
        else:
            cid = D_info[sid]['climid']

        setup = SetupXML( ystart, yend, lat=D_info[sid]['lat'], lon=D_info[sid]['lon'], \
                id=sid, cid=cid, acid=0, \
                mana=D_info[sid]['manaid'], area=D_info[sid]['area'], \
                region=D_info[sid]['region'], province=D_info[sid]['province'])
        
        # gdds
        speciespara = SpeciesParaXML( ystart, yend, \
                lat=D_info[sid]['lat'], lon=D_info[sid]['lon'], id=sid)

        PARS1 = [('MAX_TDD',     D_info[sid]['gdd1']*0.9, None),
                 ('GDDBUDSTART', D_info[sid]['gdd1']*0.5, None)]
        PARS2 = [('MAX_TDD',     D_info[sid]['gdd2']*0.9, None),
                 ('GDDBUDSTART', D_info[sid]['gdd2']*0.5, None)]
        PARS3 = [('MAX_TDD',     D_info[sid]['gdd3']*0.9, None),
                 ('GDDBUDSTART', D_info[sid]['gdd3']*0.5, None)]

        # get max. value (surrogate for missing corn info)
        PARS_CORN = [max([x[0] for x in [PARS1, PARS2, PARS3]], key=operator.itemgetter(1)),
                     max([x[1] for x in [PARS1, PARS2, PARS3]], key=operator.itemgetter(1)),
                     ('GRAINCN', 45.0,  None),
                     ('GRAIN',    0.5,  None),
                     ('ROOT',     0.15, None),
                     ('SLAMAX',  25.0,  None),
                     ('SLAMIN',  25.0,  None),
                     ('STRAW',    0.35, None),
                     ('VCMAX25', 80.0,  None),
                     ('FALEAF',   0.65, None),
                     ('FYIELD',   0.2,  None),
                     ('DOC_RESP_RATIO', 0.5, None),
                     ('MAINTENANCE_TEMP_REF', 35.0, None),
                     ('MC_LEAF',  0.015,None),
                     ('MC_STEM',  0.005,None),
                     ('MC_ROOT',  0.001,None),
                     ('MC_STORAGE', 0.0005, None),
                     ('N_DEF_FACTOR', 1.0, None),
                     ('NCFRTOPT', 0.005, None),
                     ('NCSAPOPT', 0.007, None),
                     ('NCFOLMIN', 0.007, None),
                     ('NCFOLOPT', 0.01,  None)]





        # run 2016-08-01
        #PARSALL = [('N_DEF_FACTOR', 1.5, None),
        #           ('SLAMAX',      16.0, None),
        #           ('SLAMIN',      16.0, None),
        #           ('GRAIN',      0.425, None),
        #           ('STRAW',      0.425, None),
        #           ('ROOT',        0.15, None),
        #           ('VCMAX25',     60.0, None)]

        # run 2016-08-02
        PARSALL = [('N_DEF_FACTOR', 1.5, None),
                   ('SLAMAX',      17.0, None),
                   ('SLAMIN',      17.0, None),
                   ('GRAIN',      0.45, None),
                   ('STRAW',      0.4, None),
                   ('ROOT',        0.15, None),
                   ('VCMAX25',     50.0, None)]


        PARS1 += PARSALL
        PARS2 += PARSALL
        PARS3 += PARSALL

        speciespara.addSpecies('IR72', 'IR72', [], group='crop', master='create')
        speciespara.addSpecies('IR72-0', 'IR72', PARS1, group='crop', master='append')
        speciespara.addSpecies('IR72-1', 'IR72', PARS2, group='crop', master='append')
        speciespara.addSpecies('IR72-2', 'IR72', PARS3, group='crop', master='append')
        speciespara.addSpecies('IR72-3', 'IR72', PARS1, group='crop', master='append')
        speciespara.addSpecies('IR72-4', 'IR72', PARS2, group='crop', master='append')
        speciespara.addSpecies('IR72-5', 'IR72', PARS3, group='crop', master='append')
        speciespara.addSpecies('IR72-6', 'IR72', PARS1, group='crop', master='append')
        speciespara.addSpecies('IR72-7', 'IR72', PARS2, group='crop', master='append')
        speciespara.addSpecies('IR72-8', 'IR72', PARS3, group='crop', master='append')

        # add corn info
        speciespara.addSpecies('CORN', 'CORN', PARS_CORN, group='crop', master='none')


        # (1) site file

        # litter layer
        addFlag[coord] = False

        # changed check: added thickness check
        if data2[0]['topd'] >= 0.0 and ( (data2[0]['botd'] - data2[0]['topd'])*10 > 0):
            addFlag[coord] = True


        # only 1 layer !!!
        for l in range(3):
            if data2[l]['topd'] >= 0.0:
                data2[l]['depth'] = (data2[l]['botd'] - data2[l]['topd'])*10
                if l in [0,1]:
                    split = 10
                elif l in [2,3]:
                    split = 4
                else:
                    split = 2
                data2[l]['split'] = split

                # default iron percentage
                data2[l]['iron']  = 0.01

                data2[l].pop('topd')
                data2[l].pop('botd')
                
                site.addSoilLayer( data2[l], litter=False, accuracy=cmap )

        if addFlag[coord] == True:
            sites.append( site )
            setups.append( setup )

            # for gdds   
            speciesparas.append( speciespara )

    # merge site chunks into common site file
    xml = ET.Element( "ldndcsite" )
    for scnt, site in enumerate(sites):
        x = site.xml.find("description")
        if scnt == 0:
            xml.append( x )
        a = site.xml
        a.remove( x )
        xml.append( a )
    strOut = MD.parseString(ET.tostring(xml)).toprettyxml()
    open(OUTPATH + "/VN_arable_%s_site.xml" % options.sclass, 'w' ).write( strOut )

    # merge setup chunks into common setup file
    xml = ET.Element( "ldndcsetup" )
    for scnt, setup in enumerate(setups):
        if scnt == 0:
            xml.append(setup.xml.find("description"))
        xml.append(setup.xml.find("setup"))

    strOut = MD.parseString(ET.tostring(xml)).toprettyxml()
    # setup file now differ in their event ids, located in mana subfolder  !!!    
    open(OUTPATH + "/VN_arable_%s" % options.mana + "/VN_arable_%s_%s_setup.xml" % (options.mana, options.sclass), 'w' ).write( strOut )


    # only output if D1 (!!!)

    if options.sclass == 'D1':
        print( 'Bulding species parameters file' )
        # merge speciesparameter chunks into common file
        xml = ET.Element( "ldndcspeciesparameters" )
        for scnt, speciespara in enumerate(speciesparas):
            if scnt == 0:
                xml.append(speciespara.xml.find("description"))
            #print  MD.parseString( ET.tostring(speciespara.xml) ).toprettyxml()
            xml.append(speciespara.xml.find("speciesparameters"))

        strOut = MD.parseString(ET.tostring(xml)).toprettyxml()
        open(OUTPATH + "/VN_arable_speciesparams.xml", 'w' ).write( strOut )


    # (B) create project file
    # =======================
    print( "   building project file" )

    cscen_str = ""
    if options.cscen != None:
        cscen_str += "_" + options.cscen

    if options.tempadjust:
        cfname = "../VN_arable_climate_%s%s_tds.txt" % (options.climarchive, cscen_str)
    else:
        cfname = "../VN_arable_climate_%s%s.txt" % (options.climarchive, cscen_str)

    # setup file now differ in their event ids !!!

    data = {'setup': 'VN_arable_%s_%s_setup.xml' % (options.mana, options.sclass) ,
             'site': '../VN_arable_%s_site.xml' % options.sclass    ,
             'clim': cfname,
            'event': 'VN_arable_%s_mana.xml' % options.mana,
          'airchem': '../VN_arable_airchem.txt',
         #'sitepara': '../VN_arable_siteparams.xml',
      'speciespara': '../VN_arable_speciesparams.xml'}

    Dinfo = {'opath': OUTPATH, 'mana': "VN_arable_%s" % options.mana, 'soil': options.sclass}

    try:
        os.mkdir('%(opath)s/%(mana)s' % Dinfo)
    except:
        pass
    
    project = ProjectXML( ystart, yend, lat=-1, lon=-1, days=days, sinkprefixcounter="%03r-" )
    project.addFiles( data, sourceprefix='%(opath)s/%(mana)s/' % Dinfo, \
                            sinkprefix='%(opath)s/%(mana)s/%(mana)s_output/%(mana)s_%(soil)s_' % Dinfo )
    project.write(filename='%(opath)s/%(mana)s/%(mana)s_%(soil)s.xml' % Dinfo )


    # (C) create joined climate file
    print( "   building joined climate file" )
    climate_file_path = os.path.join(OUTPATH, "VN_arable_climate_%s%s.txt" % (options.climarchive, cscen_str))

    if options.tempadjust == True and options.noclimate == False:
        # individual climate for all sids
        # range constraint to 1980 - 2050 for the moment due to size
        YEARS = range(1995,2013)

        print( "   temperature adjustment requested, this will take a while" )
        print( "   climate data limited to %d-%d" % (YEARS[0], YEARS[-1]) )
        
        # modify file name, due to size write directly to zip file
        climate_file_path = climate_file_path[:-4] + "_tds.txt" 
        
        # adjust temp by adiabatic lapse rate and elevation model
        ADIABATICLR = 0.65    # 0.5 degrees by 100m

        bar = progressbar.ProgressBar(maxval=len(sids), term_width=80, \
                widgets=[progressbar.Bar('=', ' %s [' % "   modifying climate", ']'), ' ', progressbar.Percentage()]).start()

        for cnt, sid in enumerate(sids):
            topodiff   = D_info[sid]['topodiff']

            climate_id = D_info[sid]['climid']

            if cnt == 0:
                # open file
                fout = gzip.open( climate_file_path + '.gz', 'w')

            climate_dbfile = os.path.join(climate_dbpath, "climate_%08d.txt.gz" % climate_id)
            climate_lines = gzip.open( climate_dbfile, 'rb' ).readlines()
            
            

            # parse climate file and modify temp(s)
            if cnt == 0:
                # write global section (only once)
                globalHeader = """%%global
	time     = "%d-1-1/1"
	%%cuefile = "cuefile.txt"

""" % YEARS[0]
                # write all header lines
                for lcnt, line in enumerate(climate_lines):
                    if "climate\n" in line:
                        break

                fout.write( "".join( climate_lines[0:lcnt-1] ) )
                fout.write( '%%cuefile = "cuefile.txt\n\n')

            # remove those lines
            climate_lines = climate_lines[lcnt-1:]

            # write all lines until we see the %data flag
            for lcnt, line in enumerate(climate_lines):
                if "id" in line:
                    climate_lines[lcnt] = "        id = %d\n" % sid

                if "elevation" in line:
                    climate_lines[lcnt] = "        elevation = %d\n" % D_info[sid]['elevation']

                if "latitude" in line:
                    climate_lines[lcnt] = "        latitude  = %f\n" % D_info[sid]['lat']

                if "longitude" in line:
                    climate_lines[lcnt] = "        longitude = %f\n" % D_info[sid]['lon']

                if "data\n" in line:
                    break

            # fix elevation line
            #for i, l in enumerate(climate_lines[0:lcnt]):
            #    if "elevation" in l:
            #        climate_lines[i] = "        elevation = %d" % 

            COMMENT = """    # Data height-adjusted (TOPODIFF)
    #   original climid (0.5deg res): %d"
    #   new id == cid (native res):   %d"
    #
    #   elevation difference to orginal [m]: %d"
    #     adiabatic correction [deg / 100m]: %.2f"
    #   
    #   WARNING: tavg & tamp are not corrected (for now)
""" % (climate_id, sid, D_info[sid]['topodiff'], ADIABATICLR)
            
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
            df = df[(df.year >= YEARS[0]) & (df.year <= YEARS[-1])]

            # adjust temp columns
            if df['tavg'].mean(axis=0) > -99:
                df['tavg'] = df['tavg'] - (ADIABATICLR * (topodiff / 100.0))
                df['tavg'] = df.tavg.round(1)
            
            df['tmin'] = df['tmin'] - (ADIABATICLR * (topodiff / 100.0))
            df['tmax'] = df['tmax'] - (ADIABATICLR * (topodiff / 100.0))
            df['tmin'] = df.tmin.round(1)
            df['tmax'] = df.tmax.round(1)

            out = df.to_csv(None, sep="\t", index=False, header=False) + '\n'
            fout.write( out )

            bar.update(cnt)

        bar.finish()

        fout.close()

    elif options.noclimate == False:
        # default clim (0.5 deg)
        for cnt, climate_id in enumerate(clids):
            climate_id = int(climate_id)
            if cnt == 0:
                # open file
                fout = open( climate_file_path, 'w')
            
            climate_dbfile = os.path.join(climate_dbpath, "climate_%08d.txt.gz" % climate_id)
            climate_lines = gzip.open( climate_dbfile, 'rb' ).readlines()

            if cnt > 0:
                climate_lines = climate_lines[4:]   # ignore first four lines (global section, unsafe!!!)

            fout.write( "".join( climate_lines ) )
        
        fout.close()
    else:
        print( "   skipping climate files as requested" )


    #print "   merging management and parameter files from province database"
    #os.system("python mergeProvinceData.py")



class MyParser( OptionParser ):
    def format_epilog(self, formatter):
        return self.epilog


if __name__=="__main__":

    parser = MyParser( "usage: %prog [options]", epilog=
"""

Phillipines GHG Emission Inventory Preprocessor (v1)

Use this tool to create model driving data for Phillipines inventory simulations
___________________________________________
2018/08/15, christian.werner@senckenberg.de
""")

    parser.add_option("-a", "--archive", dest="climarchive", default="NoRESM1", \
            help="climate data archive: PGF2, NoRESM1, IPSL,...")
    parser.add_option("--management", dest="mana", default="pf-control", \
            help="management regime (i.e. pf-control)")
    parser.add_option("-y", "--years", dest="years", default="1999-2015", \
            help="years (i.e. 1999-2015)")
    parser.add_option("-s", "--scenario", dest="cscen", default=None, \
            help="climate scenario: rcp2p6, rcp4p5, rcp6p0, rcp8p5")
    parser.add_option("--soil", dest="sclass", default="D1", \
            help="soil class: subsoil group D1, D2, ...")
    parser.add_option("-t", dest="tempadjust", action="store_true", default=True, \
            help="produces per-pixel climate files with adjusted temperature")
    parser.add_option("--no-climate", dest="noclimate", action="store_true", default=False, \
            help="do not create climate file")


    (options, args) = parser.parse_args()

    y = options.years.split("-")


    # for transect AU_natt use high-res data

    soil_file  = "data/VN_WISESOIL_%s.nc" % options.sclass
    ref_file   = "data/VN_MISC4.nc"
    gdd_file   = "data/VN_gdds.nc"

    climate_dbpath = "clim_pgf2/db/"
    #climate_dbpath = "/Volumes/Drobo/data/global/climate/db/%s/%s" % (options.climarchive, options.cscen)
    #climate_dbpath = "/Volumes/Drobo/data/global/climate/db/%s/%s" % (options.climarchive, options.cscen)

    if os.path.isdir( climate_dbpath ) == False:
        # check if climate data for this archive and scenario
        # exists at this location
        print( "Could not find climate db folder:" )
        print( climate_dbpath )
        exit(1)

    # Lat, lon tuples to be processed
    coords = []

    # get simulation ids & coords
    with( Dataset(ref_file, 'r') ) as nc:
        lats       = nc.variables["lat"][:]
        lons       = nc.variables["lon"][:]
        simmask    = nc.variables["simmask"][:]
        cellids    = nc.variables["cid"][:]
        climids    = nc.variables["climid"][:]
        area       = nc.variables["area_ha"][:]
        regionid   = nc.variables["regionid"][:]
        provinceid = nc.variables["provinceid"][:]
        ricefrac   = nc.variables["rice_fr"][:]
        uplafrac   = nc.variables["upla_fr"][:]
        elevation  = nc.variables["elevation"][:]
        rot        = nc.variables["rice_rot"][:]

        if options.tempadjust:
            topodiff = nc.variables["topodiff"][:]

    with( Dataset(gdd_file, 'r') ) as nc:
        gdds1      = nc.variables["gdd1"][:]
        gdds2      = nc.variables["gdd2"][:]
        gdds3      = nc.variables["gdd3"][:]

    # populate info and other dictionaries
    D_info = {}

    allids = []
    coords = []
    idxs   = []
    sids   = []
    clids  = []
    
    Dprovince = {}
    Dregion   = {}

    Lprovince = []
    Lregion   = []
    Lpids     = []
    Lrids     = []

    for cnt, line in enumerate(open("../gridlist_vietnam64.txt", 'r')):
        if cnt > 0:
            s   = line[:-1].split()
            x   = s[0].split("-")
            pid = int(s[1])
            rid = int(s[2])
            if x[0] not in Lregion:
                Lregion.append(x[0])
                Lrids.append(rid)

            if x[1] not in Lprovince:
                Lprovince.append(x[1])
                Lpids.append(pid)

    Dregion   = dict(zip( Lrids, Lregion ))
    Dprovince = dict(zip( Lpids, Lprovince ))


    # some climate ids are not available: use neighbors instead
    DuseOtherClim = {}
    
    # disabled for the moment
    #DuseOtherClim[107137] = 107136
    #DuseOtherClim[104253] = 104252
    #DuseOtherClim[117213] = 116492


    allids = []

    A_rcid = np.array(range(720*360)).reshape((360,720))
    
    l_1ks = []
    for j in range(360):
        for i in range(720):
            l_1ks.append(j*1000+i)
    A_1kid = np.array(l_1ks).reshape((360,720))

    def convert_cid_system(running_cid):
        ''' take a consecutive climid and return the 1k equivalent '''
        pos = np.nonzero(A_rcid == running_cid)
        return int(A_1kid[pos])


    for i in range(len(simmask)):
        for j in range(len(simmask[0])):
            if simmask[i,j] == 1:
                cid = int(cellids[i,j])
                if climids[i,j] > 0 and climids[i,j] < 10000000:
                    clid = int(climids[i,j])


                    # we now use the general formula (ignore clid from VN_MISC4)
                    #lons = np.arange(-179.75, 180.0, 0.5)
                    #lats = np.arange(89.75,-90, -0.5)

                    #ix = np.where( lons == cd.lon.values )[0]
                    #jx = np.where( lats == cd.lat.values )[0]
                    
                    # global climid index

                    clid = convert_cid_system(clid)

                    # replace clid with valid clid
                    if clid in DuseOtherClim.keys():
                        clid = DuseOtherClim[clid]

                    clids.append( clid )
                    D_info[cid] = {'lat': lats[i], 'lon': lons[j], \
                            'climid': clid, 'x': j, 'y': i, \
                            'area': area[i,j], \
                            'region': Dregion[regionid[i,j]], \
                            'province': Dprovince[provinceid[i,j]], \
                            'elevation': elevation[i,j]
                            }

                    if options.tempadjust:
                        D_info[cid].update( {'topodiff': topodiff[i,j]} )

                    # gdds
                    D_info[cid].update( {'gdd1': gdds1[i,j]} )
                    D_info[cid].update( {'gdd2': gdds2[i,j]} )
                    D_info[cid].update( {'gdd3': gdds3[i,j]} )

                    allids.append(int(provinceid[i,j]))
                    sids.append(cid)
                    coords.append( (lats[i], lons[j]) )
                    idxs.append( (i,j) )

    clids = list(set(clids))

    # call the main function
    OUTPATH="../projects/regional/VN_arable"

    main(sids, clids, coords, idxs, options, OUTPATH)

    exit(0)
