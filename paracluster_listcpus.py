#!/usr/bin/python

# list_cpus.py
#
# list by node the jobs running on the paracluster and show number of free/used processors per node
# written by: Wendy Williams wwilliams@strw.leidenuniv.nl
# last modified: 25 feb 2015


import os
import re
import sys
import numpy as np

os.system('pbsnodes > tempnodelist')
os.system('showq > tempqlist')

with open('tempnodelist') as f:
    lines = f.readlines() #grep "jobs = " templist
    
with open('tempqlist') as f:
    qlines = f.readlines() #grep "jobs = " templist
    
#nodelist = [item.strip() for item in lines if re.match("para..\.strw\.leidenuniv\.nl", item)]

#jobslist = [item.replace("     jobs = ","") for item in lines if re.match("     jobs = *", item)]

## take care if there is an empty node
jobslist = []
nodelist = []
jobsline = ""
for line in lines:
    
    if re.match("para..\.strw\.leidenuniv\.nl", line):
        node = line.strip()
    if re.match("     jobs = *", line):
        jobsline = line.replace("     jobs = ","").strip()
    if line.strip() == "":
        jobslist.append(jobsline)
        nodelist.append(node)
        node = ""
        jobsline = ""
        
       

   
qjobnamelist = [item.strip().split()[0] for item in qlines if re.search("Running", item)]
qjobnamelist_BatchHold = [item.strip().split()[0] for item in qlines if re.search("BatchHold", item)]
qjobuserlist = [item.strip().split()[1] for item in qlines if re.search("Running", item)]
qjobuserlist_BatchHold = [item.strip().split()[1] for item in qlines if re.search("BatchHold", item)]
  

if len(jobslist) != len(nodelist):
    print "Error: jobslist doesn't match nodelist"
    sys.exit()

for i in range(len(nodelist)):
    jobs = jobslist[i]
    if len(jobs) > 0:
        jobs = jobs.split(",")
    joblist = []
    uniquejoblist = []
    for j in jobs:
        thisjob = j.strip().replace(".para33.strw.leidenuniv.nl","").split("/")[1]
        joblist.append(thisjob)
        if thisjob not in uniquejoblist:
            uniquejoblist.append(thisjob)
            
            
    tally = [joblist.count(job)  for job in uniquejoblist]
    
    print
    print "NODE = {node} - {nused} used ({nfree} free)".format(node=nodelist[i],nused=len(joblist),nfree=64-len(joblist))
    #print "JOBS = "
    for j in range(len(uniquejoblist)):
        if uniquejoblist[j] in qjobnamelist:
            user = qjobuserlist[qjobnamelist.index(uniquejoblist[j])]
        elif uniquejoblist[j] in qjobnamelist_BatchHold:
            user = qjobuserlist_BatchHold[qjobnamelist_BatchHold.index(uniquejoblist[j])] + " - BatchHold"
        else:
            user = "?"
        print " {job:10s} - {count:3d} ({user})".format(job=uniquejoblist[j], count=tally[j], user=user)
        
        
os.system("rm -rf tempnodelist")
os.system("rm -rf tempqlist")
