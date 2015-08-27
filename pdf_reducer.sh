#!/bin/bash

# Usage: pdf_reducer.sh file.pdf

echo "Filesize old: "$(stat -c%s "$1")
pdftops $1
ps2pdf  -dPDFSETTINGS=/ebook ${1/pdf/ps}
rm ${1/pdf/ps}
echo "Filesize new: "$(stat -c%s "$1")
echo "Done."
