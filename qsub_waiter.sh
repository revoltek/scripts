#!/bin/bash
# run the argument in qsub and return when the process is finished
# check if the called program might use more CPUs or not

shopt -s expand_aliases
sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls

# choose number of processors to reserve
case $@ in
    *calibrate-stand-alone*) proc=1 ;;
    *NDPPP*) proc=6 ;;
    *wsclean*) proc=6 ;;
    *casa*) proc=6 ;;
    *awimager*) proc=6 ;;
    *) proc=1
esac

echo "Number of processors: $proc"

# call the command and capture the stdout
id=`qsub << EOF
#!/bin/bash
#PBS -N fdgjob
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=$proc
#PBS -o ./out
#PBS -e ./err
$@
EOF`

# grep the process id
id=`echo $id | cut -f1 -d"."`
echo $@

# status starts as Q (queue) and becomes E (executing)
status=`qstat $id | grep $id | awk '{print $5}'`

while [ "$status" != 'C' ] # while $status is not "C" (complete)
    do
        sleep $sleep_time
        status=`qstat $id | grep $id | awk '{print $5}'`
    done
