#!/bin/bash
# run the argument in qsub and return when the process is finished
# check if the called program might use more CPUs or not

sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
shopt -s expand_aliases

# get number of processors
proc=$1
shift

# escaped command
escapedcmd=`echo $@ | LC_ALL=C sed -e 's/[^a-zA-Z0-9,._+@%/-]/\\\&/g;'`

# ugly workaround to catch qsub errors and resubmit
until [[ $id =~ ^[0-9]+$ ]]; do

# For DEBUG output add:
#PBS -o output-\$PBS_JOBID
# instead of /dev/null
# To limit to a node, e.g. para11, use:
#PBS -l nodes=1:ppn=$proc:para11

    cmd="""#!/bin/bash
#PBS -M astro@voo.it
#PBS -m a
#PBS -N LBApipe
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=$proc
#PBS -j oe
#PBS -o /dev/null
cd \$PBS_O_WORKDIR
echo \"\${PBS_JOBID} (proc:${proc}, node: \`/bin/hostname\`) - \'${escapedcmd}\'\" >> commands.log
source /net/para34/data1/oonk/tjd_upd/lofim.sh
export PATH=\"${PATH}\"
export PYTHONPATH=\"${PYTHONPATH}\"
${@}"""

    # call the command and capture the stdout
    id=`qsub -q lofarq /dev/stdin << EOF | perl -pe 's:^\D+(\d+).*$:$1:'
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
