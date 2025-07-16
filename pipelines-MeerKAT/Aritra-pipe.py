import casatasks as casa
import numpy as np
from casacore import tables as t
import os
import sys
import logging

###########################################################
## THIS PART COULD BE CONVERTED AS INPUT THROUGH PARFILE ##
###########################################################

# User defined variables


msFile = '/path/to/ms/file.ms'
#msFile = '/home/abasu/1234567890_sdp_l0.ms'

FluxCal = 'FluxCal'
PolCal = 'PolCal'
PhaseTargetDic = {'PhaseCal1':'Target1', 'PhaseCal2':'Target2', 'PhaseCal3':'Target3'} # PhaseCal <--> Target pairs

inpath = '/path/to/cal.tables/' # use ./ for $PWD

myuvrange = '>150m' # uv range for calibration
ref_ant = 'm008'    # reference antenna [should be the same as for 1GC Stokes I calibration]

do_split = False

tab = {'K_tab': 'delay.cal',
       'B_tab': 'bandpass.cal',
       'G_tab': 'gain.cal',
       'T_tab': 'Tgain.cal',
       'F_tab': 'fluxscale.cal',
       # Pol cal tables #
       'Kcross_tab': 'kcross.cal',
       'Xf_tab': 'Xf.cal',
       'Dfdc_tab': 'Dfdc.cal',
       'Df_tab': 'Df.cal'}


###########################################################
###########################################################

# Setting actual paths
for name in tab:
    tab[name] = os.path.join(inpath, tab[name])

# Setting up logfiles
log_file = os.path.join(msFile + '.log')
casa_log = os.path.join(msFile + '_casa.log')

log_level = logging.DEBUG
logging.basicConfig(filename=log_file,
        format='%(asctime)s %(name)s:%(funcName)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M', filemode='w',
        level=log_level)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logger = logging.getLogger(__name__)

logger.info('Simple polarization analysis script')
logger.info('Written by: Aritra Basu [TLS, Tautenburg]; on: 11 July 2025')
logger.info('       (a cloudy, cool noon, Berlin)')
logger.info('** Prepared for polarization workshop in Bologna (14-18 July 2024) **')


def rename_casalog(casa_log):
    """
    By default, casa creates a log file casataks is imported
    Rename this file to casa_log and delete the old file
    """
    # Get the name of the current casa log
    old_log_filename = casa.casalog.logfile()
    # Point the casa logsink to the new file
    casa.casalog.setlogfile(filename=casa_log)
    # Delete the old file
    os.remove(old_log_filename)


def split_into_spw(ms_name):
    """
    Split the input MS into 16 SPW. To save space, the output MS replaces
    the input MS [DEACTIVATED].
    """
    temp_ms = ms_name.replace(ms_name.split('/')[-1], ms_name.split('/')[-1]+'16spw')
    casa.mstransform(vis=ms_name, outputvis=temp_ms, datacolumn='data', \
                     regridms=True, nspw=16, spw='0:210~3841')
    #delete_item(ms_name)
    #os.rename(temp_ms, ms_name)
    return ms_name.split('/')[-1]+'16spw'


def delete_item(item):
    """
    Delete the specified file or directory
    """
    if os.path.isfile(item):
        os.remove(item)
    else:
        shutil.rmtree(item)


rename_casalog(casa_log)

############## Here 1GC calibration goes in ##############

if do_split:
    msFile = split_into_spw(msFile)

casa.gaincal(vis=msFile, field=FluxCal, gaintype=’K’, caltable=tab['K_tab'])

casa.gaincal(vis=msFile, field=FluxCal, gaintype=’G’, calmode=’p’, caltable=tab['G_tab'], gaintable=[tab['K_tab']])

casa.bandpass(vis=msFile, field=FluxCal, bandtype=’B’, caltable=tab['B_tab'], gaintable=[tab['K_tab'], tab['G_tab']])

casa.gaincal(vis=msFile, field=CalField, gaintype=’T’, calmode=’ap’, caltable=tab['T_tab'], gaintable=[tab['B_tab'], tab['K_tab'], tab['G_tab']])

casa.fluxscale(vis=msFile, reference=FluxCal, fluxtable=F_tab, caltable=tab['T_tab'], transfer=PolCal,GainCal)

#casa.setjy(vis=msFile, field=PolCal, standard=YourFavModel, fluxdensity=[poly], spix=[poly], polindex=[poly], polangle=[poly])

############## FLAGS ##############
# Flagging using tricolor



###################################################
####     Polarization calibration starts      #####
###################################################

# Solve for KCROSS
logger.info('')
logger.info('Solving Kcross with %s' %PolCal)
casa.gaincal(vis=msFile,
    caltable=tab['Kcross_tab'],
    field=PolCal,
    uvrange=myuvrange,
    refant=ref_ant,
    solint='inf',
    parang=True,
    combine='',
    gaintype='KCROSS',
    gaintable=[tab['B_tab'], tab['G_tab'], tab['K_tab'], tab['T_tab']],
    gainfield=['', '', '', PolCal],
    interp=['linear', 'linear', 'nearest', 'nearest'])

# Solve for polarization angle
# Xf that is constant within a scan
logger.info('')
logger.info('Solving cross-hand calibration with %s' %PolCal)
casa.polcal(vis=msFile,
   caltable=tab['Xf_tab'],
   field=PolCal,
   uvrange=myuvrange,
   solint='inf',
   poltype='Xf',
   refant=ref_ant,
   combine='scan',
   preavg=-1.,
   gaintable=[tab['B_tab'], tab['G_tab'], tab['K_tab'], tab['Kcross_tab'], tab['T_tab']],
   gainfield=['', '', '', '', PolCal],
   interp=['linear', 'linear', 'nearest', '', 'nearest'])

# Solve for leakage
# First solve for a DC term that is constant across all scans
# Then solve for a freq dependent term that is constant across a scan
logger.info('')
logger.info('Solving DC leakage across scans with %s' %FluxCal)
casa.polcal(vis=msFile,
   caltable=tab['Dfdc_tab'],
   field=FluxCal,
   uvrange=myuvrange,
   combine='scan',
   solint='inf',
   poltype='D',
   refant='',
   gaintable=[tab['B_tab'], tab['G_tab'], tab['K_tab'], tab['Kcross_tab'], tab['Xf_tab'], tab['T_tab']],
   gainfield=[FluxCal, FluxCal, FluxCal, '', '', FluxCal],
   interp=['nearest', 'nearest', 'nearest', '', '','nearest'])
logger.info('')
logger.info('Solving leakage with %s' %FluxCal)
casa.polcal(vis=msFile,
   caltable=tab['Df_tab'],
   field=FluxCal,
   uvrange=myuvrange,
   combine='',
   solint='inf',
   poltype='Df',
   refant='',
   gaintable=[tab['B_tab'], tab['G_tab'], tab['K_tab'], tab['Dfdc_tab'], tab['Kcross_tab'], tab['Xf_tab'], tab['T_tab']],
   gainfield=[FluxCal, FluxCal, FluxCal, '', '', '', FluxCal],
   interp=['nearest', 'nearest', 'nearest', '', '', '', 'nearest'])


###################################################
####     Polarization calibration apply       #####
###################################################

# Apply to FluxCal 
logger.info('')
logger.info('Applying solutions to FluxCal: %s' %FluxCal)
casa.applycal(vis=msFile,
  field=FluxCal,
  gaintable=[tab['K_tab'], tab['B_tab'], tab['F_tab'], tab['Dfdc_tab'], tab['Df_tab'], tab['Kcross_tab'], tab['Xf_tab'], tab['G_tab']],
  interp=['nearest', 'linear', 'nearest', '', '', '', '', 'nearest'],
  gainfield=['', '', FluxCal, '', '', '', '', ''],
  parang=True, calwt=False, flagbackup=False)


# Apply to PolCal 
logger.info('')
logger.info('Applying solutions to PolCal: %s' %PolCal)
casa.applycal(vis=msFile,
  field=PolCal,
  gaintable=[tab['K_tab'], tab['B_tab'], tab['F_tab'], tab['Dfdc_tab'], tab['Df_tab'], tab['Kcross_tab'], tab['Xf_tab'], tab['G_tab']],
  interp=['nearest', 'linear', 'nearest', '', '', '', '', 'nearest'],
  gainfield=['', '', PolCal, '', '', '', '', ''],
  parang=True, calwt=False, flagbackup=False)


# Apply to PhaseCal & Targets 
for source in PhaseTargetDic:
    logger.info('')
    logger.info('Applying solutions to PhaseCal: %s' %source)
    casa.applycal(vis=msFile,
        field=source,
        gaintable=[tab['K_tab'], tab['B_tab'], tab['F_tab'], tab['Dfdc_tab'], tab['Df_tab'], 
                   tab['Kcross_tab'], tab['Xf_tab'], tab['G_tab']],
        interp=['nearest', 'linear', 'linear', '', '', '', '', 'nearest'],
        gainfield=['', '', source, '', '', '', '', ''],
        parang=True, calwt=False, flagbackup=False)

    logger.info('')
    logger.info('Applying solutions to Target: %s' %PhaseTargetDic[source])
    casa.applycal(vis=msFile,
        field=PhaseTargetDic[source],
        gaintable=[tab['K_tab'], tab['B_tab'], tab['F_tab'], tab['Dfdc_tab'], tab['Df_tab'],
                   tab['Kcross_tab'], tab['X_tab'], tab['G_tab']],
        interp=['nearest', 'linear', 'linear', '', '', '', '', 'nearest'],
        gainfield=['', '', source, '', '', '', '', ''],
        parang=True, calwt=False, flagbackup=False)
