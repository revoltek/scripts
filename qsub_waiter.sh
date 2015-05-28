#!/bin/bash
# run the argument in qsub and return when the process is finished
# check if the called program might use more CPUs or not

sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
shopt -s expand_aliases

# choose number of processors to reserve
case $@ in
    *calibrate-stand-alone*) proc=1 ;;
    *NDPPP*) proc=1 ;;
    *wsclean*) proc=5 ;;
    *casa*) proc=5 ;;
    *awimager*) proc=5 ;;
    *) proc=1
esac

# ugly workaround to catch qsub errors and resubmit
until [[ $id =~ ^-?[0-9]+$ ]]; do

    cmd="#!/bin/bash
#PBS -N fdg
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=$proc
#PBS -j oe
#PBS -o output-\$PBS_JOBID
source /home/lofar/init-lofar.sh
source /home/lofar/init-lofar-test.sh
cd \$PBS_O_WORKDIR
echo \"\$PBS_JOBID - $@\" >> commands.log
$@"

    # call the command and capture the stdout
    id=`qsub << EOF
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
