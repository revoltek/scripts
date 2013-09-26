#!/bin/bash                                                                             
                                                                                        
# ./make_manipimage.sh dir                                                            
                                         
dir=${1}                                 

do_convolve=false
res=50
do_regrid=false
do_shift=true
rainc='-4.84813e-6' #+1 arcsec in rad
decinc='3.2036e-5' #-6.5 arcsec in rad
do_rescale=true
do_average=true

prefix=''
suffix=''

if  $do_convolve ; then
  for img in `ls -d $dir/${prefix}*.image${suffix}`; do
    echo CONVOLVING: $img - resolution: $res
    casapy --nogui --nologger -c casa-convolve.py $img $res
  done
  suffix=${suffix}-convolved${res}
fi

if  $do_regrid ; then
  for img in `ls -d $dir/${prefix}*.image${suffix}`; do                                                          
    echo REGRIDDING: $img
    casapy --nogui --nologger -c casa-regrid.py $img                     
  done
  suffix=${suffix}-regridded        
fi

if  $do_shift ; then
  for img in `ls -d $dir/${prefix}*.image${suffix}`; do                                                          
    echo SHIFT: $img - ra: $rainc - dec: $decinc
    casapy --nogui --nologger -c casa-shift.py $img $rainc $decinc
  done
  suffix=${suffix}-shifted
fi

if  $do_rescale ; then
  for img in `ls -d $dir/${prefix}*.image${suffix}`; do                                                          
    echo RESCALING: $img
    casapy --nogui --nologger -c casa-rescale.py $img                     
  done
  suffix=${suffix}-rescaled
fi

i=1                                                                                     
if  $do_average ; then
  for regexp in {00[0-9],01[0-9],02[0-9],03[0-9],04[0-9],05[0-9],06[0-9],07[0-9],08[0-9],09[0-9],10[0-9],11[0-9],12[0-9],13[0-9],14[0-9],15[0-9],16[0-9],17[0-9],18[0-9],19[0-9],20[0-9],21[0-9],22[0-9],23[0-9],24[0-9]}; do
    echo "AVERAGING: " `ls -d $dir/${prefix}*${regexp}.image${suffix}`
    image_average.py $dir/${prefix}*${regexp}.image${suffix}
    mv averaged.img $dir/avg${i}${prefix}.image${suffix}-averaged
    i=$(( i+1 ))
  done
  suffix=${suffix}-averaged
fi

rm casapy*
