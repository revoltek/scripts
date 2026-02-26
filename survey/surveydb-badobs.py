#!/usr/bin/python3

import os, sys, glob, pickle

dirs_to_check = ['amp_AF','amp_BF','amp_BP','amp_RES','iono_iono','pa_amp1d','pa_phase','iono_ph']

def move_to_bad(bad_obs_ids, dir_to_check):
    for file in glob.glob(dir_to_check+'/*png'):
        bad_id = int(file.split('/')[-1].split('_-_')[0][2:])
        if bad_id in bad_obs_ids:
            if not os.path.exists('bad/%i' % bad_id):
                os.system('mkdir %s' % ('bad/%i' % bad_id))
            os.system('mv %s bad/%i' % (file, bad_id))
            print('move: %s -> bad/%i' % (file, bad_id))


bad_obs_ids = [int(d.split('/')[1]) for d in glob.glob('bad/[0-9]*')]

#bad_obs_ids += [834468,2002342,827912,2002321,828950,828670,828870,854892,844722,844956,2001946,843486,844592,847236,846084,846784,846154,844532,842870,843246,844946,856844,845380,842910,844712,842298,844886,844866,845410,842930,846774,843476,846074,844936,841418,842980,843346,842214,854322,844612,844652,844702,844876,841568,842900,843426,843466,841478,844926,844996,2002410,840096,2002459,2002211,854882,844642,844692,857598,843456,843554,836831,828860,841086,839796,841140,838890,836821,831968,833148,831958,831668,835435,828086,855786,856954,844732,844582,846794,844622] # LC16_004
#bad_obs_ids += [845086,846714,845732,847996,845686,846674,845056,845798,845186,791676,791726,786033,789424,801616,791966,792096,791796,792066,792746,792036,792226,792106,844552,801706,801746,845066,847438,802360,801736,846624,847478,845096,855268,845714,846644,845742,845650,848016,801716,801606,846634,801756,847448,801766,845106,853972,844542,845156,845610,803744,803894,845778,854372,845166,845116,847986,845046,845076,847458,846704,845788,845676,845126,845176,845696,846654,790736,854052,787165,855178,786983,790666,789214,786525,790896,787637,786615,787647,787737,784095] # LC14_002

files_bad = glob.glob('bad/*png')
bad_obs_ids += [int(file_bad.split('/')[-1].split('_-_')[0][2:]) for file_bad in files_bad]
bad_obs_ids = sorted(list(set(bad_obs_ids)))
print('Bad IDs:', bad_obs_ids)

for dir_to_check in dirs_to_check+['bad']:
    print('Check dir: %s' % dir_to_check)
    move_to_bad(bad_obs_ids, dir_to_check)

with open('/home/fdg/scripts/survey/LoLSS-bkp/badobsids.pickle', 'wb') as f:
    pickle.dump(bad_obs_ids, f)
