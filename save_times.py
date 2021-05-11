#################################################################################
# A set of functions to create and manipulate date arrays
# typically in the format [d,m,y]
#################################################################################
def save_time(start_year, end_year, ctype='standard'):

	#creates day, mon, and year for the given date range
	#ctype should be either 'standard' (default), 'noleap', or '360_day'

	import numpy as np
	#import numexpr as ne

	if ctype in ['365-day', '365_day']:
   		ctype = 'noleap'

	elif ctype in ['proleptic_gregorian', 'gregorian']:
    		ctype = 'standard'


	num_years = end_year - start_year + 1

	if ctype == '360_day':
    		ndy = 360
	else:
    		ndy = 366


	length_time = num_years*ndy

	#create array for year
	years = np.arange(start_year,end_year+1)
	year = np.empty((ndy,num_years))

	for i in range(num_years):
    		year[:,i] = np.tile(years[i],(ndy))

	year = np.reshape(year, (length_time), order='F')

	mon_nday = [31,29,31,30,31,30,31,31,30,31,30,31]	
	#create array for mon
	if ctype == '360_day':
    		q = 1
    		mon_1 = np.empty((360,))
    		for k in np.arange(12):
        		mon_1[q:q+30-1,1] = np.tile(k+1,(30))
        		q = q+30
    
	else:
		
			yr_days1 = [np.tile(h+1,mon_nday[h]) for h in range(12)]
			mon_1 = np.hstack(yr_days1)
#    		jan = np.tile(1,(31))
#    		feb = np.tile(2,(29))
#    		mar = np.tile(3,(31))
#    		apr = np.tile(4,(30))
#    		may = np.tile(5,(31))
#    		jun = np.tile(6,(30))
#    		jul = np.tile(7,(31))
#    		aug = np.tile(8,(31))
#    		sep = np.tile(9,(30))
#    		oct = np.tile(10,(31))
#    		nov = np.tile(11,(30))
#    		dec = np.tile(12,(31))
#    
#    		mon_1 = np.concatenate((jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec))


	mon = np.tile(mon_1, (num_years))

	#create array for day

	if ctype == '360_day':
    		day_2 = np.tile(np.arange(1,30),(12,1))
    		day_1 = np.reshape(day_2.T, (360,), order='F')
	else:

			day_0 = [np.arange(mon_nday[h])+1 for h in range(12)]
			day_1 = np.hstack(day_0)
#    		jan_d = np.arange(31)+1
#    		feb_d = np.arange(29)+1
#    		mar_d = np.arange(31)+1
#    		apr_d = np.arange(30)+1
#    		may_d = np.arange(31)+1
#    		jun_d = np.arange(30)+1
#    		jul_d = np.arange(31)+1
#    		aug_d = np.arange(31)+1
#    		sep_d = np.arange(30)+1
#    		oct_d = np.arange(31)+1
#    		nov_d = np.arange(30)+1
#    		dec_d = np.arange(31)+1
#    
#    		day_1 = np.concatenate((jan_d.T, feb_d.T, mar_d.T, apr_d.T, may_d.T, jun_d.T, jul_d.T, aug_d.T, sep_d.T, oct_d.T, nov_d.T, dec_d.T))


	day = np.tile(day_1, (num_years))

	#remove the 29-Feb except for leap years
	if ctype == 'standard':
		xind = np.ones((len(day)))
		xind[((day == 29) & (mon == 2)) & (((year % 4) != 0) | (((year % 100)==0) & ((year % 400)!=0)))] = 0
		xind = xind.astype('bool')
		year = year[xind]
		mon = mon[xind]
		day = day[xind]
    
	elif ctype == 'noleap':
		xind = np.ones((len(day)))
		xind[((day == 29) & (mon == 2))] = 0
		xind = xind.astype('bool')
		year = year[xind]
		mon = mon[xind]
		day = day[xind]


	dmy = np.array((day, mon, year)).T
	dmy = dmy.astype(int)
	
	return dmy


###########################################################################################
def save_time_m(start_year, end_year, start_mon=1, end_mon=12):

	#creates mon and year for the given date range
	#optional input for starting and ending months, if full year isn't needed

	import numpy as np

	num_years = end_year - start_year + 1

	length_time = num_years*12

	#create array for year
	years = np.arange(start_year,end_year+1)
	year = np.empty((12,num_years))

	for i in range(num_years):
    		year[:,i] = np.tile(years[i],(12))

	year = np.reshape(year, (length_time), order='F')

	#create array for mon
	mons = np.arange(1,13)
	mon = np.tile(mons, num_years)


	my = np.array((mon, year)).T
	my = my.astype(int)

	if start_mon != 1:
		ind = np.ones((len(my[:,0])))
		ind[(my[:,1]==years[0])&(my[:,0]<start_mon)] = 0.
		ind = ind.astype('bool')
		my = my[ind,:]
	if end_mon != 12:
		ind = np.ones((len(my[:,0])))
		ind[(my[:,1]==years[-1])&(my[:,0]>end_mon)] = 0.
		ind = ind.astype('bool')
		my = my[ind,:]

	return my

##############################################################################
def time_labels(dmy):

	#creates string labels for each year of daily data, to be used 
	# as axis labels for plotting
	
	import numpy as np
	
	#determine the index value of each new year
	d = np.diff(dmy[:,-1])
	ind = np.array(np.nonzero(d==1))+1
	ind = np.squeeze(ind.T)
	ind = np.insert(ind, 0, 0)

	yr_st = dmy[ind,-1]
	yr_st = yr_st.astype('|S4')

	return (ind, yr_st)

###############################################################################
def save_nc_format(dmy):
	
	# takes a dmy array and outputs an array of double in the absolute time
	#  format for netcdf files ("day as %Y%m%d.%f" which is 20020425.0)
	#  this is accepted as a valid input for the CDOs	
	
	import numpy as np

	nd,n = np.shape(dmy)
	
	dmys = dmy.astype(str)
	
	date1 = np.zeros((nd,),dtype='S12')
	
	for i in range(nd):
		if dmy[i,0]<10:
			if dmy[i,1]<10:
				date1[i] = dmys[i,2]+'0'+dmys[i,1]+'0'+dmys[i,0]
			else:
				date1[i] = dmys[i,2]+dmys[i,1]+'0'+dmys[i,0]
		elif dmy[i,1]<10:
			date1[i] = dmys[i,2]+'0'+dmys[i,1]+dmys[i,0]
		else:
			date1[i] = dmys[i,2]+dmys[i,1]+dmys[i,0]
		
		if n==4:
			date1[i] += str(dmy[i,3]/24.)[1:]
	
	date2 = date1.astype(float)
	
	return date2
	
###############################################################################
def save_time_h(start_date, end_date, step=1, start_hour=0, ctype='standard'):

	# returns a dmyh array (h = hour)
	# start_date and end_date should be in the form [d m y]
	# accepts input for the hour step (e.g., 1, 3, 12, etc.)
	# will only work properly if start_hour is an even step from 0

	import numpy as np
	
	dmy = save_time(start_date[2], end_date[2], ctype=ctype)
	
	ind1 = int(np.nonzero((dmy[:,0]==start_date[0])&(dmy[:,1]==start_date[1])&(dmy[:,2]==start_date[2]))[0])
	ind2 = int(np.nonzero((dmy[:,0]==end_date[0])&(dmy[:,1]==end_date[1])&(dmy[:,2]==end_date[2]))[0])
	
	dmy1 = dmy[ind1:ind2+1,:]	

	nd = np.size(dmy1[:,0])	
	
	if np.mod(24,step)!=0:
		raise ValueError('step must be a factor of 24')
	nh = 24/step
	
	dmy2 = np.repeat(dmy1,nh,axis=0)

	hrs = np.arange(0,24,step)
	hrs1 = np.tile(hrs,nd)
	rm = start_hour/step
		
	dmyh1 = np.hstack((dmy2,hrs1[:,None]))
	if rm==0:
		dmyh = dmyh1
	else:
		dmyh = dmyh1[rm:-(nh-rm),:]
	
	return dmyh

###############################################################################

def dmy_to_julianday(dmy, use_leap=1):
	
	# converts input dmy array into list of Julian Day values
	# if the use_leap flag is set to 1, will determine if leap days and adjust
	#    accordingly; if not, assumes 365-day years
	# currently not appropriate for 360-day calendars
	
	import numpy as np
	import datetime

	dmy = np.copy(dmy)	
	
	nd = np.size(dmy[:,0])
	
	jds = np.zeros((nd,2))
	jds[:,1] = dmy[:,2]
	
	#if use_leap is false, replace all years with non-leap year 1997
	if use_leap==0:
		dmy[:,2] = 1997
	
	for i in range(nd):
		dt = datetime.datetime(dmy[i,2],dmy[i,1],dmy[i,0])
		jds[i,0] = dt.timetuple().tm_yday
	
	jds = jds.astype('int')
	
	return jds
	
###############################################################################
	
def remove_leap(x,dmy):
	
	#removes leap days from x, assuming time is the first dimension
	
	import numpy as np
	
	if np.size(x,axis=0) != np.size(dmy,axis=0):
		raise ValueError('size of dmy must match size of time dimension of x')
	
	ind = np.ones((np.size(x,axis=0)))
	ind[(dmy[:,1]==2) & (dmy[:,0]==29)] = 0
	ind = ind.astype('bool')
	
	x1 = np.copy(x)
	x1 = x1[ind,...]
	dmy1 = np.copy(dmy)[ind,:]
	
	return x1,dmy1
	
###############################################################################

def convert_nctime(df1, calendar='standard', flag='dmy'):
	
	#reads the time from a netcdf file and converts to a dmy array
	#df1 is the file read with netCDF4.Dataset
	#flag indicates output array; default is dmy, but my and y are also options

	import numpy as np
	from netCDF4 import num2date

	time = df1.variables['time'][:]
	units = df1.variables['time'].units
	
	dates = num2date(time,units=units,calendar=calendar)
	
	n = len(dates)
	
	dmy = np.zeros((n,3))
	for i in range(n):
		dmy[i,0] = dates[i].day
		dmy[i,1] = dates[i].month
		dmy[i,2] = dates[i].year
		
	dmy = dmy.astype('int')
	
	if flag=='y':
		return dmy[:,2]
	elif flag=='my':
		return dmy[:,1:]
	else:
		return dmy

###############################################################################
	
	
	
	
	
	








