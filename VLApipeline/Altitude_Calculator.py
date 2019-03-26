import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun

npoints=1000

Obj_name='J2355+4950'

#Obj=SkyCoord.from_name('ICRF J235509.4+495008')

#Obj=SkyCoord.from_name('ICRF J172727.6+453039')

#Obj=SkyCoord.from_name('PSZ1 G108.18-11.53')

#Obj=SkyCoord.from_name('MACSJ1752.0+4440')

Obj=SkyCoord.from_name('3c286')

#Obj=SkyCoord.from_name('3c48')

#Obj=SkyCoord.from_name('3c147')

#Obj=SkyCoord.from_name('3c138')

#Obj=SkyCoord.from_name('J071338.1+434917')

#print Phase_cal

#Phase_cal_ra='23:55:09.458169'
#Phase_cal_dec='+49:50:08.34001'

#Sun_ra='22:13:30.0'
#Sun_dec='-10:59:42.0'

VLA_longitude=-107.61777778*u.deg

VLA=EarthLocation(lat=34.07874917*u.deg,lon=VLA_longitude,height=2124*u.m)

#print VLA

#Test_longitude=0*u.deg

#Time at VLA = UTC - 7

utcoffset=-7*u.hour ######Check UTC and sidereal time difference

date='2019-02-20'

midnight=Time(date+' 00:00:00', scale='utc')

#midnight=0*u.hour+0*u.minute+0*u.second

GMT_at_midnight=midnight-utcoffset #GMT at VLA midnight i.e. 07:00:00 GMT

#print midnight

time_range=np.linspace(0,24,npoints)*u.hour

VLA_time=midnight+time_range 

GMT_time=GMT_at_midnight+time_range

#print GMT_time

LST=GMT_time.sidereal_time('apparent',longitude=VLA_longitude) #Local sidereal time from GMT
#print LST[0]
#print LST

Obj_altaz=Obj.transform_to(AltAz(obstime=GMT_time,location=VLA)) #Phase_cal AltAz coordinates using GMT

#print Phase_cal_altaz.alt[0]

Obj_alt=Obj_altaz.alt

Sun=get_sun(GMT_time)

#print Sun[0]

Sun_altaz=Sun.transform_to(AltAz(obstime=GMT_time,location=VLA))

Sun_alt=Sun_altaz.alt
'''
sep=Obj.separation(Sun)

#print sep

diff_alt=Obj_alt[0]-Sun_alt[0]
diff_az=Obj_altaz.az[0]-Sun_altaz.az[0]
dist=np.sqrt(diff_alt**2+diff_az**2)
print dist
'''

'''

max_sun=max(Sun_alt)

ind=-1
for i in Sun_alt:
	ind+=1
	if i==max_sun:
		max_ind=ind
	else:
		continue

print VLA_time[max_ind]
'''

'''
max_alt=max(Obj_alt)

ind=-1
for i in Phase_cal_alt:
	ind+=1
	if i==max_alt:
		max_ind=ind
	else:
		continue
'''
#Obs_start=(16*u.hour+13*u.minute)
#Obs_end=(6*u.hour+32*u.minute)
Obs_start=(10*u.hour+59*u.minute)
Obs_end=(0*u.hour+45*u.minute)
Obs_start_array=Obs_start*np.linspace(1,1,npoints)
Obs_end_array=Obs_end*np.linspace(1,1,npoints)

min_elevation=(15*u.deg)*np.linspace(1,1,npoints)

vertical_line=np.linspace(0,90,npoints)


xmarks=np.array(range(25))
plt.plot(LST,Obj_alt,linestyle='none', marker='.',markersize=1,color='green',label=Obj_name)
#plt.plot(GMT_time,Obj_alt,linestyle='none', marker='.',markersize=1,color='green',label=Obj_name)
#plt.plot(LST,Sun_alt,linestyle='none',marker='.',markersize='0.5',color='orange',label='Sun')
plt.plot(Obs_start_array,vertical_line,color='black')
plt.plot(Obs_end_array,vertical_line,color='black')
plt.plot(LST,min_elevation,color='red')
plt.axvspan(Obs_end/u.hour,Obs_start/u.hour,facecolor='black',alpha=0.5)
plt.xticks(xmarks)
plt.xlim(0,24)
plt.xlabel('LST')
plt.ylim(0,90)
plt.ylabel('Altitude [deg]')
plt.legend(loc='upper left')
plt.show()




'''
def convert_ra_dec_to_degrees(ra,dec):
	ra_col_count=0
	dec_col_count=0
	ra_hr,ra_min,ra_sec='','',''
	dec_deg,dec_min,dec_sec='','',''
	tot_hr=24
	tot_min=24*60
	tot_sec=24*60*60
	for a in ra:
		if ra_col_count<1:
			if a==':':
				ra_col_count+=1
			else:
				ra_hr+=a
		elif ra_col_count==1:
			if a==':':
				ra_col_count+=1
			else:
				ra_min+=a			
		else:
			ra_sec+=a
	ra_hr_deg=float(ra_hr)*360/tot_hr
	ra_min_deg=float(ra_min)*360/tot_min
	ra_sec_deg=float(ra_sec)*360/tot_sec
	ra_deg=ra_hr_deg+ra_min_deg+ra_sec_deg
	print ra_deg

	for b in dec:
		if dec_col_count<1:
			if b==':':
				dec_col_count+=1
			else:
				dec_deg+=b
		elif dec_col_count==1:
			if b==':':
				dec_col_count+=1
			else:
				dec_min+=b			
		else:
			dec_sec+=b
	#dec_deg=dec_deg[1:-1]
	#print dec_deg,dec_min,dec_sec
	#if '+' in dec_deg:
			
	dec_deg_deg=float(dec_deg)
	dec_min_deg=float(dec_min)/60
	dec_sec_deg=float(dec_sec)/3600
	dec_full_deg=dec_deg_deg+dec_min_deg+dec_sec_deg	
	print dec_full_deg

convert_ra_dec_to_degrees(Phase_cal_ra,Phase_cal_dec)
convert_ra_dec_to_degrees(Sun_ra,Sun_dec)'''
