#!/usr/bin/env python3
# cehck LTA and write the cycle for which we have observations in the grid file

import os, sys, re, pickle
import numpy as np
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
from astropy.utils import iers
iers.IERS_A_URL='https://maia.usno.navy.mil/ser7/finals2000A.all'
iers_a = iers.IERS_A.open(iers.IERS_A_URL)
iers.earth_orientation_table.set(iers_a)
from astropy.table import Table
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u

survey_projects = 'LT16_004,LT14_002,LC12_017,LC9_016,LC8_031,LC18_020' # list of projects related with the LBA survey
projects = survey_projects.split(',')

lofar_location = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m) # LOFAR

# obs ids that are know to have bad data
bad_obs_ids = [834468,2002342,827912,2002321,828950,828670,828870,854892,844722,844956,2001946,843486,844592,847236,846084,846784,846154,844532,842870,843246,844946,856844,845380,842910,844712,842298,844886,844866,845410,842930,846774,843476,846074,844936,841418,842980,843346,842214,854322,844612,844652,844702,844876,841568,842900,843426,843466,841478,844926,844996,2002410,840096,2002459,2002211,854882,844642,844692,857598,843456,843554,836831,828860,841086,839796,841140,838890,836821,831968,833148,831958,831668,835435,828086,855786,856954,844732,844582,846794,844622] # LC16_004
bad_obs_ids += [845086,846714,845732,847996,845686,846674,845056,845798,845186,791676,791726,786033,789424,801616,791966,792096,791796,792066,792746,792036,792226,792106,844552,801706,801746,845066,847438,802360,801736,846624,847478,845096,855268,845714,846644,845742,845650,848016,801716,801606,846634,801756,847448,801766,845106,853972,844542,845156,845610,803744,803894,845778,854372,845166,845116,847986,845046,845076,847458,846704,845788,845676,845126,845176,845696,846654,790736,854052,787165,855178,786983,790666,789214,786525,790896,787637,786615,787647,787737,784095] # LC14_002
bad_obs_ids += [2002720,2002893,2002900,2003137,2003474,2003537,2004173,2004180,2004187,2004194,2004201,2004208,2004215,2004222,2004229,2004244,2004259,2004266,2004273,2004288,2004323,2004330,2004411,2004552,2004559,2004869,2004925,2004932,2005019,2005236,2005243,2005547,2005771,2005778,2005792,2005829,2006414] # to be removed from LTA

# The class of data to query
cls = CorrelatedDataProduct
re_cal = re.compile('.*3[c|C](196|296|295|380).*') # 3c296 is a typo in LTA

if not os.path.exists("update_allsky-grid.pickle"):
    obs_all = []
    for project in projects:
            print('Checking project: %s' % project)
            query_observations = Observation.select_all().project_only(project)
            for observation in query_observations:
                obs_id = int(observation.observationId)
                antennaset = observation.as_dict()['Observation.antennaSet']
    
                print('Checking obs_id: %i' % obs_id)
                for subarray in observation.subArrayPointings:
                    subarray_dict = subarray.as_dict()
                    field_id = subarray_dict['SubArrayPointing.targetName'].split('_')[-1]
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
                        obs_all.append([field_id,'bug',obs_id,0,antennaset,field_ra,field_dec])
                    elif obs_id in bad_obs_ids:
                        print('Add BAD obs to the list: %s' % (field_id))
                        obs_all.append([field_id,'bad',obs_id,0,antennaset,field_ra,field_dec])
                    else: 
                        print('Add obs to the list: %s (LST: %f)' % (field_id, lst.hour))
                        obs_all.append([field_id,project,obs_id,lst,antennaset,field_ra,field_dec])

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
            obs_all.append([obs,'comm2020',obs_id,-1,'LBA Sparse Even',ra,dec])

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
for obs in obs_all:
    obs[0] = obs[0].strip()
    # fix some name errors
    if obs[0] == 'PP033+66': obs[0] = 'P033+66' # has a wrong name in LTA
    if obs[0] == 'PP219+37': obs[0] = 'P219+37' # has a wrong name in LTA
    if obs[0] == '094+59': obs[0] = 'P094+59' # has a wrong name in LTA
    if obs[0] == 'P142+49': obs[0] = 'P142+42' # has a wrong name in LTA

    try:
        idx = np.where(grid['name'] == obs[0].upper())[0][0]
    except:
        print('WARNING: missing %s in the grid' % obs)
        continue

    dist = SkyCoord(grid['ra'][idx]*u.deg,grid['dec'][idx]*u.deg).separation(SkyCoord(obs[5]*u.deg,obs[6]*u.deg))
    if dist > 1*u.arcmin:
        print('WARNING: wrong coord for %s (should be: %f %f - it is: %f %f)' % (obs[0],grid['ra'][idx],grid['dec'][idx],obs[5],obs[6]))
        continue

    if not obs[1] == 'bug' and not obs[1] == 'bad': grid['hrs'][idx] += 1
    idxcell = list(grid['obsid'][idx]).index(0)
    #print(obs,idxcell)
    grid['cycle'][idx][idxcell] = obs[1]
    grid['obsid'][idx][idxcell] = obs[2]
    grid['antset'][idx][idxcell] = obs[4]
    try:
        grid['LST'][idx][idxcell] = obs[3].hour
    except:
        grid['LST'][idx][idxcell] = 0


grid.write('allsky-grid.fits', overwrite=True)
