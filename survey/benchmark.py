#!/usr/bin/python3

import os, sys, re, glob
from statistics import mean

fields = None
outerorsparse = '' # "o" OR "s"
fields = ['P022+84o']

if fields is None:
    fields = glob.glob('P*'+outerorsparse)
else:
    fields = [field for field in fields if outerorsparse in field]

def duration_to_hours(s):
    """
    Parse a string like '1d 5h 54m 3s', '1 day, 5h 54m', or '0h 05m' into total hours (float).
    Examples:
    print(duration_to_hours("dd-parallel_P127+39s # 1d 5h 54m 3s"))  # includes seconds
    print(duration_to_hours("quality_P127+39s # 0h 05m"))          # 0.0833...
    print(duration_to_hours("foo # 15m"))                         # 0.25
    """
    # take only part after '#'
    if '#' in s:
        s = s.split('#', 1)[1].strip()

    # regex with optional groups for days, hours, minutes, seconds
    # supports "1d", "1 day", "1 day,", "5h", "54m", "3s" (case-insensitive)
    pattern = r'(?:(\d+)\s*(?:d|day|days))?[, ]*\s*(?:(\d+)\s*h)?\s*(?:(\d+)\s*m)?\s*(?:(\d+)\s*s)?'
    m = re.search(pattern, s, flags=re.IGNORECASE)
    if not m:
        return 0.0  # no match

    days = int(m.group(1) or 0)
    hours = int(m.group(2) or 0)
    minutes = int(m.group(3) or 0)
    seconds = int(m.group(4) or 0)

    return days * 24 + hours + minutes / 60.0 + seconds / 3600.0

def find_all_sections(sections_dict, s):
    """
    Return a list of all keys whose patterns appear as substrings in s.
    """
    matches = []
    for key, patterns in sections_dict.items():
        for pat in patterns:
            if pat.lower() in s.lower():   # case-insensitive
                matches.append(key)
                break                      # avoid duplicate key if multiple patterns match
    return matches

durations = {'copy':[],'timesplit':[],'dd-serial':[],'dd-parallel':[],'quality':[], 'saveproducts':[]}
sections_ddp = {
    'imaging':['imaging','image'], 
    'solving':['solve'] ,
    'predict':['init_model','add_patches','intreg_predict'],
    'apply':['corrupt_iono', 'corrFR', 'corr_sf_sols', 'final_fr_corr','subfield_corr_tec', 'correct-sidelobe'],
    'extregion':['xtreg'],
    'subfield':['subfield'],
    'faraday':['fr','FR'],
    'amp':['amp'],
    'flag':['flag_residuals_lr'],
    'data_manipulation':['addcol', 'subtract']}
fract_ddp = {s:[] for s in sections_ddp.keys()}
sections_dds = {
    'imaging':['image','imaging','lres','lressub'], 
    'solving':['calibrate'] ,
    'predict':['predict'],
    'output':['leakage','output'],
    'flag':['flag'],
    'data_manipulation':['add_columns', 'subtract', 'fullsub', 'shift', 'beamcorr', 'merge_h5', 'interpsol', 'remove_col'],
    'final imaging':['imaging'],
    'final imaging lres':['lres']}
fract_dds = {s:[] for s in sections_dds.keys()}

for field in fields:

    if not os.path.exists(field+'/logs/PiLL.walker'):
        print(f'ERROR: missing {field}/logs/PiLL.walker')
        continue
   # open and process each line
    with open(field+'/logs/PiLL.walker') as f:
        print(f"Workin on {field}")
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            name = line.split('#', 1)[0].split('_')[0].strip()
            duration_hrs = duration_to_hours(line)
            print(f"{name}: {duration_hrs:.2f} h") 
            if duration_hrs == 0: continue # skip missing Saveproducts
            durations[name].append(duration_hrs)

    # dd-parallel
    field_durations_ddp = {s:[] for s in sections_ddp.keys()}
    if not os.path.exists(field+'/logs/pipeline-ddparallel.walker'):
        print(f'ERROR: missing {field}/logs/pipeline-ddparallel.walker')
        continue
    with open(field+'/logs/pipeline-ddparallel.walker') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            name = line.split('#', 1)[0].strip()
            duration_hrs = duration_to_hours(line)
            #print(f"{name}: {duration_hrs:.2f} h") 
            if duration_hrs == 0: continue # skip non contributions

            matches = find_all_sections(sections_ddp, name)
            for match in matches:
                field_durations_ddp[match].append(duration_hrs)
            if len(matches) == 0: print(f'ERROR: no match for {name}!')
        
        # sum up all times per section
        for section in sections_ddp.keys():
            sumsection = sum(field_durations_ddp[section])
            fract_ddp[section].append(sumsection/durations['dd-parallel'][-1])
            #print(f"{section}: {sumsection} h [out of {durations['dd-parallel'][-1]} h]")

    # dd-serial
    field_durations_dds = {s:[] for s in sections_dds.keys()}
    if not os.path.exists(field+'/logs/pipeline-ddserial.walker'):
        print(f'ERROR: missing {field}/logs/pipeline-ddserial.walker')
        continue
    with open(field+'/logs/pipeline-ddserial.walker') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            name = line.split('#', 1)[0].strip()
            duration_hrs = duration_to_hours(line)
            #print(f"{name}: {duration_hrs:.2f} h") 
            if duration_hrs == 0: continue # skip non contributions

            matches = find_all_sections(sections_dds, name)
            for match in matches:
                field_durations_dds[match].append(duration_hrs)
            if len(matches) == 0: print(f'ERROR: no match for {name}!')
        
        # sum up all times per section
        for section in sections_dds.keys():
            sumsection = sum(field_durations_dds[section])
            fract_dds[section].append(sumsection/durations['dd-serial'][-1])
            #print(f"{section}: {sumsection} h [out of {durations['dd-serial'][-1]} h]")

means = {k: (mean(v) if v else 0.0) for k, v in durations.items()}

print("Means (hours):")
for k, m in means.items():
    print(f"{k}: {m:.2f} h")

print("Benchmark DDP (%):")
for k, m in fract_ddp.items():
    percent = mean(m)*100
    print(f"{k}: {percent:.2f}%")

print("Benchmark DDS (%):")
for k, m in fract_dds.items():
    percent = mean(m)*100
    print(f"{k}: {percent:.2f}%")
