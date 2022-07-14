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
from astropy.coordinates import EarthLocation
from astropy.time import Time
import astropy.units as u

survey_projects = 'LT16_004,LT14_002,LC12_017,LC9_016,LC8_031' # list of projects related with the LBA survey
projects = survey_projects.split(',')

lofar_location = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m) # LOFAR

# The class of data to query
cls = CorrelatedDataProduct
re_cal = re.compile('.*3[c|C](196|295|380).*')

if not os.path.exists("update_allsky-grid.pickle"):
    obs_all = []
    for project in projects:
            print('Checking project: %s' % project)
            query_observations = Observation.select_all().project_only(project)
            for observation in query_observations:
                obs_id = int(observation.observationId)
    
                print('Checking obs_id: %i' % obs_id)
                for subarray in observation.subArrayPointings:
                    subarray_dict = subarray.as_dict()
                    field_id = subarray_dict['SubArrayPointing.targetName'].split('_')[-1]

                    # ignore calibrators
                    if re_cal.match(field_id):
                        continue

                    # ignore period with correlator problems
                    time = subarray_dict['SubArrayPointing.startTime']
                    # find HA
                    astrotime = Time(time, format='datetime', scale='utc')
                    lst = astrotime.sidereal_time('mean', lofar_location.lon)
                    if time.year == 2021 and ( (time.month==2 and time.day>=8) or (time.month>2 and time.month<8) or ( time.month==8 and time.day<=3) ):
                        print('Add BAD obs to the list: %s' % (field_id))
                        obs_all.append([field_id,'bad',obs_id,0])
                    else: 
                        print('Add obs to the list: %s (LST: %f)' % (field_id, lst.hour))
                        obs_all.append([field_id,project,obs_id,lst])

    # add manual things:

    # Commissioning obs on SPARSE 12/2020
    for obs in ['P182+35','P184+37','P184+40','P185+35','P187+37','P187+40','P188+35','P190+37','P191+40']:
        obs_all.append([obs,'comm2020',-1,0])
        obs_all.append([obs,'comm2020',-1,0])
        obs_all.append([obs,'comm2020',-1,0])

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
for obs in obs_all:
    obs[0] = obs[0].strip()
    try:
        idx = np.where(grid['name'] == obs[0].upper())[0][0]
    except:
        print('WARNING: missing %s in the grid' % obs)
        continue
    if not obs[1] == 'bad': grid['hrs'][idx] += 1
    idxcell = list(grid['obsid'][idx]).index(0)
    #print(obs,idxcell)
    grid['cycle'][idx][idxcell] = obs[1]
    grid['obsid'][idx][idxcell] = obs[2]
    try:
        grid['LST'][idx][idxcell] = obs[3].hour
    except:
        grid['LST'][idx][idxcell] = 0

grid.write('allsky-grid.fits', overwrite=True)
