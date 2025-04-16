#!/usr/bin/env python3
# cehck LTA and write the cycle for which we have observations in the grid file

import os, sys, re, pickle
import numpy as np
#from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
#from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
#from astropy.utils import iers
#iers.IERS_A_URL='https://maia.usno.navy.mil/ser7/finals2000A.all'
#iers_a = iers.IERS_A.open(iers.IERS_A_URL)
#iers.earth_orientation_table.set(iers_a)
from astropy.table import Table
from astropy.coordinates import EarthLocation, SkyCoord, match_coordinates_sky
from astropy.time import Time
import astropy.units as u
import zlib, base64

survey_projects = 'LT16_004,LT14_002,LC12_017,LC9_016,LC8_031,LC15_011,LC18_007,LC18_020,LC20_025,LC20_039,LC20_011' # list of projects related with the LBA survey
projects = survey_projects.split(',')

lofar_location = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m) # LOFAR

# obs ids that are know to have bad data
bad_obs_ids = pickle.load(open(os.path.dirname(os.path.realpath(__file__))+'/LoLSS-bkp/badobsids.pickle', 'rb'))
#bad_obs_ids += [2002720,2002893,2002900,2003137,2003474,2003537,2004173,2004180,2004187,2004194,2004201,2004208,2004215,2004222,2004229,2004244,2004259,2004266,2004273,2004288,2004323,2004330,2004411,2004552,2004559,2004869,2004925,2004932,2005019,2005236,2005243,2005547,2005771,2005778,2005792,2005829,2006414] # to be removed from LTA

# The class of data to query
cls = CorrelatedDataProduct
re_cal = re.compile('.*3[c|C](196|296|295|380).*') # 3c296 is a typo in LTA

if not os.path.exists("update_allsky-grid.pickle"):
    obs_all = [] # [field_id,project,obs_id,lst,antennaset,field_ra,field_dec,demix,ignoretarget,hrs]
    for project in projects:
            print('Checking project: %s' % project)
            query_observations = Observation.select_all().project_only(project)
            for observation in query_observations:
                obs_id = int(observation.observationId)
                antennaset = observation.as_dict()['Observation.antennaSet']
                if 'HBA' in antennaset: continue


                
                # TODO get demix info
                #dataproduct_query = cls.observations.contains(observation)
                #dataproduct_query &= cls.isValid == 1
                #parset = dataproduct_query[0].as_dict()['CorrelatedDataProduct.pipeline.parset.content']
                #parset = zlib.decompress(base64.b64decode(''.join(parset))).decode()
                #parset = dict(x.split("=",1) for x in re.split("ObsSW.", parset))
                #demix = parset['Observation.ObservationControl.PythonControl.PreProcessing.demix_always']
                #ignoretarget = parset['Observation.ObservationControl.PythonControl.DPPP.demixer.ignoretarget']
                demix = ''
                ignoretarget = ''

                # unfortunately there are too many cases with wrond metadata, set this by hand
                #hrs = np.rint(observation.as_dict()['Observation.duration']/3600)
                #if hrs != 1:
                #    print('%i: Strange hrs: %i' % (obs_id, hrs))
                if obs_id == 2021053: hrs = 8
                elif obs_id == 813284: hrs = 5
                elif obs_id in [2021470,2021477,2021519,2021526,2021533,2021540,2021554,2021561,2021568,2021575,2021582,2025030,2025622,2037219,2037226,2037233,2037268,2039336,2040684,2040784,833652,833662]: hrs = 4
                elif obs_id == 813294: hrs = 3
                elif obs_id in [2023765,2036150,2036157,2036164,2036171,2036178,2036192,2036199,2039326,2039361,2039368]: hrs = 2
                elif obs_id in [2006261,2007870]: 
                    hrs = 0
                    continue # do not add these obs
                else: hrs = 1

                print('Checking obs_id: %i (hrs:%i)' % (obs_id, hrs))
    
                for subarray in observation.subArrayPointings:
                    subarray_dict = subarray.as_dict()
                    field_id = subarray_dict['SubArrayPointing.targetName'].strip().replace(' ','_').split('_')[-1]
                    print(subarray_dict['SubArrayPointing.targetName'],subarray_dict['SubArrayPointing.targetName'].strip().replace(' ','_'),field_id)
                    field_ra = subarray_dict['SubArrayPointing.pointing.rightAscension']
                    field_dec = subarray_dict['SubArrayPointing.pointing.declination']

                    # ignore calibrators
                    if re_cal.match(field_id):
                        continue

                    # ignore period with correlator problems
                    time = subarray_dict['SubArrayPointing.startTime']
                    # find HA
                    astrotime = Time(time, format='datetime', scale='utc')
                    lst = astrotime.sidereal_time('mean', lofar_location.lon)

                    if time.year == 2021 and ( (time.month==2 and time.day>=8) or (time.month>2 and time.month<8) or ( time.month==8 and time.day<=3) ):
                        print('Add BUG obs to the list: %s' % (field_id))
                        obs_all.append([field_id,'bug',obs_id,0,antennaset,field_ra,field_dec,demix,ignoretarget,hrs])
                    elif obs_id in bad_obs_ids:
                        print('Add BAD obs to the list: %s' % (field_id))
                        obs_all.append([field_id,'bad',obs_id,0,antennaset,field_ra,field_dec,demix,ignoretarget,hrs])
                    else: 
                        print('Add obs to the list: %s (LST: %f)' % (field_id, lst.hour))
                        obs_all.append([field_id,project,obs_id,lst,antennaset,field_ra,field_dec,demix,ignoretarget,hrs])

    # add manual things:

    # Commissioning obs on SPARSE 12/2020
    commissioning_obs = {800556: [['P184+37', 'P187+40', 'P188+35'],
  [184.0686541666666, 187.9613333333333, 188.88933333333333],
  [37.22016916666666, 39.75032611111111, 34.73559944444445]],
 800546: [['P182+35', 'P187+37', 'P191+40'],
  [182.59013333333328, 187.32009166666663, 191.3290916666666],
  [34.70261472222222, 37.23719777777778, 39.768018611111124]],
 800536: [['P184+40', 'P185+35', 'P190+37'],
  [184.5944541666666, 185.73941249999996, 190.5722375],
  [39.73269083333333, 34.7190913888889, 37.25426666666666]],
 800586: [['P182+35', 'P187+37', 'P191+40'],
  [182.59013333333328, 187.32009166666663, 191.3290916666666],
  [34.70261472222222, 37.23719777777778, 39.768018611111124]],
 800576: [['P184+40', 'P185+35', 'P190+37'],
  [184.5944541666666, 185.73941249999996, 190.5722375],
  [39.73269083333333, 34.7190913888889, 37.25426666666666]],
 800566: [['P184+37', 'P187+40', 'P188+35'],
  [184.0686541666666, 187.9613333333333, 188.88933333333333],
  [37.22016916666666, 39.75032611111111, 34.73559944444445]],
 800616: [['P184+40', 'P185+35', 'P190+37'],
  [184.5944541666666, 185.73941249999996, 190.5722375],
  [39.73269083333333, 34.7190913888889, 37.25426666666666]],
 800606: [['P184+37', 'P187+40', 'P188+35'],
  [184.0686541666666, 187.9613333333333, 188.88933333333333],
  [37.22016916666666, 39.75032611111111, 34.73559944444445]],
 800596: [['P182+35', 'P187+37', 'P191+40'],
  [182.59013333333328, 187.32009166666663, 191.3290916666666],
  [34.70261472222222, 37.23719777777778, 39.768018611111124]]}

    for obs_id in commissioning_obs.keys():
        for obs,ra,dec in zip(commissioning_obs[obs_id][0],commissioning_obs[obs_id][1],commissioning_obs[obs_id][2]):
            obs_all.append([obs,'comm2020',obs_id,-1,'LBA Sparse Even',ra,dec,'[Cas,Cyg]','false',1])

    pickle.dump(obs_all, open( "update_allsky-grid.pickle", "wb" ))

else:
    print ("WARNING: use old pickle file.")
    obs_all = pickle.load( open( "update_allsky-grid.pickle", "rb" ) )

# now update the grid file
grid = Table.read('allsky-grid.fits')
grid['hrs'] = 0
grid['cycle'] = ''
grid['obsid'] = 0
grid['LST'] = None
grid['antset'] = ''
all_Skycoords = SkyCoord(grid['ra'],grid['dec'])
for obs in obs_all:
    obs[0] = obs[0].strip()
    # fix some typos
    if obs[0] == 'PP033+66': obs[0] = 'P033+66'
    if obs[0] == 'PP219+37': obs[0] = 'P219+37'
    if obs[0] == '094+59': obs[0] = 'P094+59'
      
    # fix some wrong names
    #if obs[1] == 'LC12_017' and obs[0] == 'P174+57': obs[0] = 'P176+60'
    if (obs[1] == 'LT14_002' or obs[1] == 'LC18_020') and obs[0] == 'P142+49': obs[0] = 'P142+42'
    if obs[1] == 'LT16_004' and obs[0] == 'P351+28': obs[0] = 'P321+28'

    try:
        idx = np.where(grid['name'] == obs[0].upper())[0][0]
    except:
        if obs[0] != 'Coma' and obs[0] != 'M31' and obs[0] != 'Cring' and not 'PSO' in obs[0] \
            and not 'Test' in obs[0]  and not 'test' in obs[0] and not 'A'in obs[0] and obs[0] != '3c48' \
            and obs[0] != '3C48' and obs[0] != 'ref' and not 'Zwcl'in obs[0] and not 'lba'in obs[0] and not 'hba'in obs[0]:
            print('WARNING: missing %s in the grid' % obs)
        continue
    
    dist = SkyCoord(grid['ra'][idx]*u.deg,grid['dec'][idx]*u.deg).separation(SkyCoord(obs[5]*u.deg,obs[6]*u.deg))
    if dist > 1*u.arcmin:
        match = match_coordinates_sky(SkyCoord(obs[5]*u.deg,obs[6]*u.deg), all_Skycoords)
        if match[1] < 1*u.arcmin:
            idx_correct = match[0]
            name_correct = grid['name'][idx_correct]
            print('WARNING: wrong coord for %s (id: %i) - %s (should be: %f %f - it is: %f %f) -> assigning to: %s (%f %f)' % (obs[0],obs[2],obs[1],grid['ra'][idx],grid['dec'][idx],obs[5],obs[6],name_correct,grid['ra'][idx_correct],grid['dec'][idx_correct]))
            idx=int(idx_correct)
        else:
            print('WARNING: wrong coord for %s (id: %i) - %s (should be: %f %f - it is: %f %f) -> no other pointing found!' % (obs[0],obs[2],obs[1],grid['ra'][idx],grid['dec'][idx],obs[5],obs[6]))
            continue

    # this fixes obs with multiple entries of the same field (e.g. 849992)
    if obs[2] in grid['obsid'][idx]: 
        print('WARNING: obsid %i already in the db for field %s. Ignoring.' % (obs[2],obs[0]) )
        continue

    # remove known problematic obsids (no data):
    if obs[2] in [790836,801468,801478,801492,805972,818060,828830,833618,2002010,2002272,2019066,2019093,2021187,2021463,2021484,2021678,2023187,2035650,2035741,2035762,2035769,2035776,2035888,2036185,2039382]: 
        print('WARNING: obsid %i has missing data. Ignoring.' % (obs[2]) )
        continue
    # remove known problematic obsids (flagged data):
    if obs[2] in [2007856, 2014278, 2014285, 2015837, 2015844, 2006261]:
        print('WARNING: obsid %i has flagged data. Ignoring.' % (obs[2]) )
        continue

    # all ok, add
    if not obs[1] == 'bug' and not obs[1] == 'bad': grid['hrs'][idx] += obs[9] # most of the time is 1
    idxcell = list(grid['obsid'][idx]).index(0)
    grid['cycle'][idx][idxcell] = obs[1]
    grid['obsid'][idx][idxcell] = obs[2]
    grid['antset'][idx][idxcell] = obs[4]
    grid['demix'][idx][idxcell] = obs[7]
    grid['ignoretarget'][idx][idxcell] = obs[8]
    try:
        grid['LST'][idx][idxcell] = obs[3].hour
    except:
        grid['LST'][idx][idxcell] = 0


grid.write('allsky-grid.fits', overwrite=True)
