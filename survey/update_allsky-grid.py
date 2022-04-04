#!/usr/bin/env python3
# cehck LTA and write the cycle for which we have observations in the grid file

import os, sys, re, pickle
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
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
                        obs_all.append([field_id,'bad',obs_id,'-'])
                    else: 
                        print('Add obs to the list: %s' % (field_id))
                        obs_all.append([field_id,project,obs_id,lst])

    # add manual things:

    # Commissioning obs on SPARSE 12/2020
    for obs in ['P182+35','P184+37','P184+40','P185+35','P187+37','P187+40','P188+35','P190+37','P191+40']:
        obs_all.append([obs,'comm2020'])
        obs_all.append([obs,'comm2020'])
        obs_all.append([obs,'comm2020'])

    pickle.dump(obs_all, open( "update_allsky-grid.pickle", "wb" ))

else:
    print ("WARNING: use old pickle file.")
    obs_all = pickle.load( open( "update_allsky-grid.pickle", "rb" ) )

# now update the grid file
grid = Table.read('allsky-grid.fits')
grid['hrs'] = 0
grid['cycle'] = ''
grid['obsid'] = ''
grid['LST'] = ''
for obs in obs_all:
    obs[0] = obs[0].strip()
    idx = (grid['name'] == obs[0].upper())
    if sum(idx) == 0: 
        print('WARNING: missing %s in the grid' % obs)
        continue
    if not obs[1] == 'bad': grid['hrs'][idx] += 1
    if grid['cycle'][idx].tolist()[0] == '':
        grid['cycle'][idx] = obs[1]
        grid['obsid'][idx] = obs[2]
        grid['LST'][idx] = obs[3]
    else:
        grid['cycle'][idx] = grid['cycle'][idx].tolist()[0]+','+obs[1]
        grid['obsid'][idx] = grid['obsid'][idx].tolist()[0]+','+obs[2]
        grid['LST'][idx] = grid['LST'][idx].tolist()[0]+','+obs[3]

grid.write('allsky-grid.fits', overwrite=True)
