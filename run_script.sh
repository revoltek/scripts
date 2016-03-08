#!/bin/bash 

# setup PYTHONPATH so it can find the pydal module

if [ -z "$PYTHONPATH" ]
then
  export PYTHONPATH=/app/usg/release/lib/python
else
  export PYTHONPATH=/app/usg/release/lib/python:${PYTHONPATH}

fi

echo ""
echo PYTHONPATH = ${PYTHONPATH}
echo ""

command='source /app/scripts/doPyrap'
# run the example
#command='python msinfo.py /lifs001/L2009_13427/SB0.MS'
#command='python closure.py SB0.MS'
command='python baseline.py /home/fdg/data/calibLBA/outdp3-L2009_13581_SB10.MS'
#command='python autoflagger.py'
#command='python uvcoverage.py SB0.MS'
echo ""
echo Running: ${command}
echo ""

$command
