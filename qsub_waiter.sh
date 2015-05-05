#!/bin/bash
# run the argument in qsub and return when the process is finished
# check if the called program might use more CPUs or not

sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
shopt -s expand_aliases

# choose number of processors to reserve
case $@ in
    *calibrate-stand-alone*) proc=1 ;;
    *NDPPP*) proc=6 ;;
    *wsclean*) proc=6 ;;
    *casa*) proc=6 ;;
    *awimager*) proc=6 ;;
    *) proc=1
esac

# call the command and capture the stdout
id=`qsub << EOF
#!/bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=$proc
#PBS -o ./out
#PBS -e ./err
source /home/lofar/init-lofar.sh
cd $PBS_O_WORKDIR
$@
EOF`

sleep 2

# grep the process id
id=`echo $id | cut -f1 -d"."`

# status starts as Q (queue) and becomes E (executing)
status=`qstat $id | grep $id | awk '{print $5}'`

while [ "$status" != 'C' ] # while $status is not "C" (complete)
    do
        sleep $sleep_time
        status=`qstat $id | grep $id | awk '{print $5}'`
    done
