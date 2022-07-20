#!/bin/sh

java -Xmx160000m -jar ~/opt/src/Aladin/Aladin.jar -hipsgen maxthread=64 in=/home/fdg/data/LBAsurvey/mosaic-dde/mosaics out=/home/fdg/data/LBAsurvey/mosaic-dde/hips id=LoLSS_DR1 creator=FrancescoDeGasperin pixelRange="-5E-4 0.003" target="202.47 47.195" INDEX TILES DETAILS
