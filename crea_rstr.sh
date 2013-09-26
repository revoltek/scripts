#!/bin/sh

# arg 1 is the model rstr file with XXX.image as the imagename to be substituted

for file in `ls -d *image`; do
	cat $1 | sed s/XXX.image/${file}/ > ${file}.rstr
	casapy --nogui --nologger -c write_jpeg.py ${file}.rstr
done
