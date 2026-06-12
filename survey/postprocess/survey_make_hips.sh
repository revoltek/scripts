#!/bin/sh

# module load aladin
# module load jdk

java -Xmx160000m -jar /iranet/soft/aladin/aladin-12.060/Aladin.jar -hipsgen maxthread=64 in=/homes/fdg/storage/surveytgts/mosaics/healpix out=/homes/fdg/storage/surveytgts/mosaics/hips id="ivo://astron.nl/P/LoLSS_DR2" creator=FrancescoDeGasperin dataRange=-0.0005,0.003 target="202.47 47.195" verbose=3 INDEX TILES DETAILS