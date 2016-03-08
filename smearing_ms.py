#!/usr/bin/env python

"""
smearing.py

This script is intended to calculate the maximum allowable channel size
(to avoid bandwidth smearing) and timeslot (to avoid time smearing).

It uses the following equations to estimate these parameters:

dt = 1200 D/B
dv = (0.1 D/B) * v

where D is the desired field size, B is the maximum baseline length,
t is time and v is frequency.

Since the field-of-view is F=1.22*l/D where l is the wavelength,
and the baseline can be expressed in lambdas as L=B/l, we can write
D/B = 1.22/(F*L) in the above

dt = 1200 * 1.22 / (F*L)
dv = (0.1 * 1.22 / (F*L)) * v

These are relaxed for the imaging case:

dt = 4127 * 1.22 / (F*L)
dv = (0.3 * 1.22 / (F*L)) * v

Both options are taken care of in the script (see the help).

When not using verbose mode, there are just two outputs (on the same line):
averaging_time [seconds] averaging_frequency [Hz]

Written by George Heald
v1.0 completed 4/11/2010
v1.1 refactored, John Swinbank, 4/11/2010
"""

import pyrap.tables as pt
import logging

ATEAM = (("19h59m28.3", "40d44m02"),
    ("23h23m24", "58d48m54"),
    ("05h34m31.95", "22d00m52.1"),
    ("12h30m49.4", "12d23m28"))

STATION_DIAMETER = {
    "LBA": 87.0,
    "HBA": 35.0,
    "LBA_INNER": 30.0}


def read_ms(logger, msname, ateam, diameter=None):
    def get_station_diameter(table):
        histable = pt.table(table.getkeyword('HISTORY'), ack=False)
        for line in histable.getcell('APP_PARAMS', 0):
            try:
                key, value = line.split("=")
            except:
                pass
            if key == "Observation.antennaSet":
                antenna_set = value
                break
        if antenna_set == "LBA_INNER":
            logger.debug("LBA_INNER mode")
            return STATION_DIAMETER["LBA_INNER"]
        elif antenna_set[:3] == "LBA":
            logger.debug("LBA_(OUTER,SPARSE,X,Y) mode")
            return STATION_DIAMETER["LBA"]
        elif antenna_set[:3] == "HBA":
            logger.debug("HBA mode")
            return STATION_DIAMETER["HBA"]
        else:
            logger.error("Failed to identify antenna set")

    def field_size_ateam(table):
        logging.debug('Computing field size for A-team')
        fieldtable = table.getkeyword('FIELD').split()[1]
        taqloutput = pt.taql("calc from %s calc max(angdist (DELAY_DIR[0,], [%s]))" % (fieldtable, ", ".join(",".join(src) for src in ATEAM))  )
        return taqloutput[0]

    def field_size_nominal(table, wavelength, diameter):
        if not diameter:
            diameter = get_station_diameter(table)
        logger.debug("Station diameter %f m" % diameter)
        return 1.22*wavelength/diameter

    t = pt.table(msname, readonly=True, ack=False)
    interval = t.getcell('INTERVAL', 0)
    swtable = t.getkeyword('SPECTRAL_WINDOW')
    tsw = pt.table(swtable, readonly=True, ack=False)
    freq = tsw.getcell('REF_FREQUENCY', 0)
    wavelength = 299792458./freq
    maxbl = pt.taql("calc sqrt(max([select sumsqr(UVW[0:1]) from %s]))" % msname)[0] / wavelength
    chwidth = tsw.getcell('CHAN_WIDTH', 0)[0]

    if ateam:
        fieldsize = field_size_ateam(t)
    else:
        fieldsize = field_size_nominal(t, wavelength, diameter)

    logger.debug('Frequency is %f MHz'%(freq/1.e6))
    logger.debug('Wavelength is %f m'%(wavelength))
    logger.debug('Maximum baseline length is %f m = %f lambdas'%(maxbl*wavelength,maxbl))
    logger.debug('Integration time is %f sec'%(interval))
    logger.debug('Channel width is %f Hz'%(chwidth))
    logger.debug('Field size is %f degrees'%(fieldsize*180./3.14159))

    return fieldsize, maxbl, freq, interval, chwidth


def calculate_sizes(logger, imaging, fieldsize, maxbl, freq, interval, chwidth):
    def smearing_imaging(fieldsize, maxbl, freq):
        logger.debug("Imaging mode")
        #fieldsize=0.5*3.14115/180.
        return 4127.*1.22/(fieldsize*maxbl), 0.3*1.22/(fieldsize*maxbl)*freq

    def smearing_noimaging(fieldsize, maxbl, freq):
        logger.debug("Calibration mode")
        return 1200.*1.22/(fieldsize*maxbl), 0.1*1.22/(fieldsize*maxbl)*freq

    if imaging:
        dt, dv = smearing_imaging(fieldsize, maxbl, freq)
    else:
        dt, dv = smearing_noimaging(fieldsize, maxbl, freq)

    if dt < interval:
        logger.warn("Warning! Averaging time (%f) less than interval."%(dt))
        dt = interval
    if dv < chwidth:
        logger.warn("Warning! Averaging frequency (%f) less than channel width."%(dv))
        dv = chwidth

    return dt, dv


def main(msname, options):
    # Run in stadalone mode
    ateam = options.ateam
    imaging = options.imaging
    if options.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")
    else:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    logger = logging.getLogger()

    dt, dv = calculate_sizes(logger, imaging, *read_ms(logger, msname, ateam))
    logger.debug("Output values are averaging_time [s] and averaging_freq [Hz]")
    logger.info("%f, %f" % (dt, dv))


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options] MSname"
    opt = optparse.OptionParser(usage)
    opt.add_option('-a','--ateam',help='Calculate field size based on A-team positions rather than nominal field size for this frequency and array? [default False]',action='store_true',default=False)
    opt.add_option('-i','--imaging',help='Calculate based on imaging case rather than calibration case? [default False]',action='store_true',default=False)
    opt.add_option('-v','--verbose',help='Give helpful output?',action='store_true',default=False)
    (options, args) = opt.parse_args()
    if len(args) != 1:
            opt.error('Incorrect number of arguments, see help.')
    main(args[0], options)

