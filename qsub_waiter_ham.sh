#!/bin/bash
# run the argument in qsub and return when the process is finished
# check if the called program might use more CPUs or not

sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
shopt -s expand_aliases

# get number of processors
proc=$1
shift

# ugly workaround to catch qsub errors and resubmit
until [[ $id =~ ^[0-9]+$ ]]; do

# For DEBUG output add:
#PBS -o output-\$PBS_JOBID
# instead of /dev/null

    cmd="#!/bin/bash
#PBS -N LBApipe
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=$proc
#PBS -j oe
#PBS -o /dev/null
echo \"\${PBS_JOBID} (proc:${proc}, node: \`/bin/hostname\`) - ${@}\" >> commands.log
source /home/lofar/init-lofar.sh
source /home/lofar/init-lofar-test.sh
PATH=\"/home/stsf309/scripts:${PATH}\"
${@}"

    # call the command and capture the stdout
    id=`qsub /dev/stdin << EOF | perl -pe 's:^\D+(\d+).*$:$1:'
$cmd
EOF`

    # grep the process id
    id=`echo $id | cut -f1 -d"."`
    sleep 1
done

# status starts as Q (queue) and becomes E (executing)
status=`qstat $id | grep $id | awk '{print $5}'`

while [ "$status" != 'C' ] # while $status is not "C" (complete)
    do
        sleep $sleep_time
        status=`qstat $id | grep $id | awk '{print $5}'`
    done
