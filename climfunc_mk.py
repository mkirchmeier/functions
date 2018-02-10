#################################################################################
# A set of functions commonly needed for climate data
# See save_times.py for the functions to create the appropriate time inputs
#################################################################################
def remove_seasonal(x, my, startyear=0, endyear=0):

	# returns values of x with monthly means (of the period specified) removed
	# x should have time as its first dimension
	# if startyear or endyear are unspecified, use bounds of my
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(my)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	if np.shape(my)[1]==3:
		my = my[:,1:]	#if dmy, remove day column

	if startyear==0:
		startyear = my[0,1]	#if no startyear is given, use first entry
	if endyear==0:
		endyear = my[-1,1]	#if no endyear is given, use last
		
	y = np.zeros(n)
	my2 = my[(my[:,1]>=startyear) & (my[:,1]<=endyear),:]
	x2 = x[(my[:,1]>=startyear) & (my[:,1]<=endyear),...]

	for k in xrange(12):
		y[my[:,0]==k+1,...] = x[my[:,0]==k+1,...] - np.nanmean(x2[my2[:,0]==k+1,...],axis=0)

	return y

#################################################################################
def remove_seasonal2(x, my, seasons=[[12,1,2],[3,4,5],[6,7,8],[9,10,11]], startyear=0, endyear=0, cflag=1, tflag=0):

	# returns values of x with seasonal means (of the period specified) removed
	# x should have time as its first dimension
	# if startyear or endyear are unspecified, use bounds of my
	# seasons should be a list of month numbers grouped into the desired seasons
	#	default is DJF,MAM,JJA,SON
	# 	0s are returned for months not included in seasons
	# if cflag is set to 1, calculate seasonal means first
	# if tflag is set to 1, also return seasonal means
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(my)[0]:
		raise ValueError('time should be the first dimension of the input data')

	if isinstance(seasons[0],list)==0:	#if only one season provided, nest list
		seasons = [seasons]
	
	if np.shape(my)[1]==3:
		my = my[:,1:]	#if dmy, remove day column

	if startyear==0:
		startyear = my[0,1]	#if no startyear is given, use first entry
	if endyear==0:
		endyear = my[-1,1]	#if no endyear is given, use last


	if cflag==1:
		x1,ys = calc_seasonal(x,my,seasons=seasons,tflag=1)

		
		y = np.zeros(np.shape(x1))
		ys2 = ys[(ys[:,0]>=startyear) & (ys[:,0]<=endyear),:]
		x2 = x1[(ys[:,0]>=startyear) & (ys[:,0]<=endyear),...]
		sm = np.zeros(np.shape(x1)); sm = sm[:len(seasons),...]	

		for k in xrange(len(seasons)):
			y[ys[:,1]==k+1,...] = x1[ys[:,1]==k+1,...] - np.nanmean(x2[ys2[:,1]==k+1,...],axis=0)
			sm[k,...] = np.nanmean(x2[ys2[:,1]==k+1,...],axis=0)

			
	else:
		y = np.zeros(n)
		my2 = my[(my[:,1]>=startyear) & (my[:,1]<=endyear),:]
		x2 = x[(my[:,1]>=startyear) & (my[:,1]<=endyear),...]
		sm = np.zeros(n); sm = sm[:len(seasons),...]
	
		for k in xrange(len(seasons)):
			y[np.in1d(my[:,0],seasons[k]),...] = x[np.in1d(my[:,0],seasons[k]),...] - np.nanmean(x2[np.in1d(my2[:,0],seasons[k]),...],axis=0)
			sm[k,...] = np.nanmean(x2[np.in1d(my2[:,0],seasons[k]),...],axis=0)


	if tflag==1:		
		return y, sm
	else:
		return y
	
#################################################################################
def remove_seasonal_dec(x, my, yds, flag=0):

	# returns the values of x with the monthly means removed by decade
	# x should have time as its first dimension
	# yds is an array of the bounds of each decade
	# flag to return array of (time,...) if 0 and (years, months, ...) if 1

	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(my)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	if np.size(n)>1:
		x1 = np.reshape(x, (-1,12,n[1]), order='C')	#reshape to (years, months, ...)
	else:
		x1 = np.reshape(x, (-1,12), order='C')
	
	y = np.empzerosty(np.shape(x1))
	y[:] = np.nan
	nd = np.size(yds)-1
	yrs2 = np.unique(my[:,1]) #list of years

	for m in np.arange(nd):
		inds = np.squeeze(np.nonzero((yrs2>=yds[m])&(yrs2<yds[m+1])))
		y[inds,...] = np.nanmean(x1[inds,...],axis=0)
	
	if flag==1:
		z = x1 - y
	else:
		y1 = np.reshape(y, n, order='C')
		z = x - y1 
	return z

#################################################################################
def remove_clim_daily(x, dmy, startyear=0, endyear=0, c_flag=0):

	# returns values of x with daily means (of the period specified) removed
	# x should have time as its first dimension
	# assumes any leap days have been removed
	# if startyear and/or endyear is 0, use first (or last)
	# if c_flag = 1, also return grid of climatologies
	
	import numpy as np
	import save_times as st

	n = np.shape(x)
	if n[0]!=np.shape(dmy)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	nd = 365
	dates = st.save_time(2000,2001,'365_day')
	
	y = np.zeros(n)

	if startyear==0:
		startyear = dmy[0,2]	#if no startyear is given, use first entry
	if endyear==0:
		endyear = dmy[-1,2]	#if no endyear is given, use last
		
	if np.size(n)==1:
		c = np.zeros((365))
		n2 = [1]
	else:
		n1 = list(n)
		n1[0] = nd
		c = np.zeros(n1)
		n2 = np.ones(int(np.size(n1)))
		n2 = list(n2)
	
	dmy2 = dmy[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),:]
	x2 = x[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),...]
	
	for i in np.arange(nd):
		
		ind1 = np.zeros(np.size(dmy[:,0]))
		ind1[(dmy[:,0]==dates[i,0])&(dmy[:,1]==dates[i,1])] = 1
		ind1 = ind1.astype('bool')
		
		ind2 = np.zeros(np.size(dmy2[:,0]))
		ind2[(dmy2[:,0]==dates[i,0])&(dmy2[:,1]==dates[i,1])] = 1
		ind2 = ind2.astype('bool')
		
		c[i,...] = np.nanmean(x2[ind2,...],axis=0)
		n2[0] = np.sum(ind1)
		y[ind1,...] = x[ind1,...] - np.tile(c[i,...],n2)

	if c_flag==1:
		return (y,c)
	else: 
		return y

#################################################################################
def remove_clim_mon(x, dmy, startyear, endyear, c_flag=0):

	# returns values of x with monthly means (of the period specified) removed
	# x should have time as its first dimension
	# assumes any leap days have been removed
	# if c_flag = 1, also return grid of climatologies
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(dmy)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	mons = np.arange(12)+1
	
	y = np.zeros(n)
	
	if np.size(n)==1:
		c = np.zeros((12))
		n2 = [1]
	else:
		n1 = list(n)
		n1[0] = 12
		c = np.zeros(n1)
		n2 = np.ones(int(np.size(n1)))
		n2 = list(n2)
	
	dmy2 = dmy[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),:]
	x2 = x[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),...]
	
	for k in np.arange(12):
		
		ind1 = np.zeros(np.size(dmy[:,0]))
		ind1[dmy[:,1]==mons[k]] = 1
		ind1 = ind1.astype('bool')
		
		ind2 = np.zeros(np.size(dmy2[:,0]))
		ind2[dmy2[:,1]==mons[k]] = 1
		ind2 = ind2.astype('bool')
		
		c[k,...] = np.nanmean(x2[ind2,...],axis=0)
		n2[0] = np.sum(ind1)
		y[ind1,...] = x[ind1,...] - np.tile(c[k,...],n2)

	if c_flag==1:
		return (y,c)
	else: 
		return y

#################################################################################
def ratio_clim_mon(x, dmy, startyear, endyear, c_flag=0):

	# returns values of x as a ratio of the monthly means (of the period specified)
	# x should have time as its first dimension
	# assumes any leap days have been removed
	# if c_flag = 1, also return grid of climatologies
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(dmy)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	mons = np.arange(12)+1
	
	y = np.zeros(n)

	if np.size(n)==1:
		c = np.zeros((12))
		n2 = [1]
	else:
		n1 = list(n)
		n1[0] = 12
		c = np.zeros(n1)
		n2 = np.ones(int(np.size(n1)))
		n2 = list(n2)
	
	dmy2 = dmy[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),:]
	x2 = x[(dmy[:,2]>=startyear) & (dmy[:,2]<=endyear),...]
	
	for k in np.arange(12):
		
		ind1 = np.zeros(np.size(dmy[:,0]))
		ind1[dmy[:,1]==mons[k]] = 1
		ind1 = ind1.astype('bool')
		
		ind2 = np.zeros(np.size(dmy2[:,0]))
		ind2[dmy2[:,1]==mons[k]] = 1
		ind2 = ind2.astype('bool')
		
		c[k,...] = np.nanmean(x2[ind2,...],axis=0)
		n2[0] = np.sum(ind1)
		y[ind1,...] = x[ind1,...] / np.tile(c[k,...],n2)

	if c_flag==1:
		return (y,c)
	else: 
		return y

#################################################################################
def seasonal_means(x, my, startyear=0, endyear=0):

	# returns the monthly means (of the period specified) of x
	# x should have time as its first dimension
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(my)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	if np.shape(my)[1]==3:
		my = my[:,1:]	#if dmy, remove day column
		
	if startyear==0:
		startyear = my[0,1]	#if no startyear is given, use first entry
	if endyear==0:
		endyear = my[-1,1]	#if no endyear is given, use last
	
	
	my2 = my[(my[:,1]>=startyear) & (my[:,1]<=endyear),:]
	
	y = np.empty(n); y = y[:12,...]
	x2 = x[(my[:,1]>=startyear) & (my[:,1]<=endyear),...]
	
	for k in xrange(12):
		y[k,...] = np.nanmean(x2[my2[:,0]==k+1,...],axis=0)

	return y

###############################################################################
def calc_monthly(x,dmy,method='mean'):
	
	#returns an array of monthly means calculated from daily data
	#x should have time as its first dimension
	#method can be 'mean' for monthly averages or 'sum' for monthly totals

	import numpy as np
	
	n = np.shape(x)
	if n[0] != np.shape(dmy)[0]:
		raise ValueError('time should be the first dimension of the input data')
	
	yrs = np.unique(dmy[:,2])
	ny = len(yrs)
	
	y1 = np.zeros((12,ny)+n[1:])
	
	for k in xrange(12):
		for j in xrange(ny):
			if method is 'mean':
				y1[k,j,...] = np.nanmean(x[(dmy[:,1]==k+1)&(dmy[:,2]==yrs[j]),...],axis=0)
			elif method is 'sum':
				y1[k,j,...] = np.nansum(x[(dmy[:,1]==k+1)&(dmy[:,2]==yrs[j]),...],axis=0)
			else:
				raise ValueError('invalid method')
				
	y = np.reshape(y1,(12*ny,)+n[1:],order='F')
	
	return y

#################################################################################
def calc_seasonal(x, my, seasons=[[12,1,2],[3,4,5],[6,7,8],[9,10,11]], tflag=0):

	# returns an array of seasonal means of x
	# x should have time as its first dimension
	# seasons should be a list containing month numbers grouped into seasons
	# if tflag is 1, return year, season array
	
	import numpy as np

	n = np.shape(x)
	if n[0]!=np.shape(my)[0]:
		raise ValueError('time should be the first dimension of the input data')

	if isinstance(seasons[0],list)==0:	#if only one season provided, nest list
		seasons = [seasons]
	
	ns = len(seasons)
	
	if np.shape(my)[1]==3:
		my = my[:,1:]	#if dmy, remove day column
		
	yrs = np.unique(my[:,1])
	ny = len(yrs)
	
	y1 = np.zeros((ns,ny)+n[1:])

	
	for k in xrange(ns):
		for j in xrange(ny):
			y1[k,j,...] = np.nanmean(x[(my[:,1]==yrs[j])&(np.in1d(my[:,0],seasons[k])),...],axis=0)

	y = np.reshape(y1,(ns*ny,)+n[1:],order='F')

	if tflag==1:
		my2 = np.zeros((ns*ny,2))
		my2[:,0] = np.repeat(yrs,ns)
		my2[:,1] = np.tile(np.arange(ns)+1,ny)
		my2 = my2.astype('int')
		
		return y,my2
	else:
		return y

#################################################################################
def runmean(x, d, ts='None', axis=0):
	# returns the running means of the data in x with window d
	# if x has more than one dimension, compute along dimension axis
	# adapted from https://gordoncluster.wordpress.com/2014/02/13/python-numpy-how-to-generate-moving-averages-efficiently-part-2/
	# if the time step array is provided, the new accompanying ts array will be output
	
	import numpy as np
	
	n = np.shape(x)
	nd = np.size(n)
	
	if nd > 3:
		raise ValueError('Currently cannot handle inputs larger than 3D. Sorry.')
	
	w = np.repeat(1.0, d)/d
	

	if nd>1:
		x1 = np.swapaxes(x,0,axis) #make axis of interest first dimension
		xa1 = np.empty(np.shape(x1))
		xa1 = xa1[:-(d-1),...] #set up array less the tail entries
		n1 = np.shape(xa1)
		
		if nd == 2:
			for i in np.arange(n1[1]):
				xa1[:,i] = np.convolve(x1[:,i],w,'valid')
		elif nd==3:
			for i in np.arange(n1[1]):
				for j in np.arange(n1[2]):
					xa1[:,i,j] = np.convolve(x1[:,i,j],w,'valid')
			
		xa = np.swapaxes(xa1,0,axis) #put the axes back in the original order
	else:
		xa = np.convolve(x, w, 'valid')
		
	if ts is 'None':	
		return xa
	else:
		q = np.floor(d/2.)
		if d%2==0:
			ts2 = ts[(q-1):-q]
			ts2 = ts2 - (ts[1]-ts[0])/2.
		else:
			ts2 = ts[q:-q]
			
		return (xa,ts2)

#################################################################################
def nonrunmean(x, yrs, d, axis=0, remove='last'):
	
	# returns the means of the data in x with non-verlapping window d
	# example, compute non-overlapping 5-year means
	# if x has more than one dimension, compute along dimension ax
	# input years of data and output will be center year of bins
	# the remove flag can take inputs 'first' or 'last' to determine whether 
	#    extra values are removed from the beginning or end of the array

	import numpy as np
	
	n = np.shape(x)
	nd = np.size(n)
	
	if nd > 3:
		raise ValueError('Currently cannot handle inputs larger than 3D. Sorry.')
	
	yrs2 = np.squeeze(yrs)
	if np.size(np.shape(yrs2))>1:
		raise ValueError('Please make sure the input yrs is a vector')
	if np.size(yrs2)!=n[axis]:
		raise ValueError('Size of yrs does not match size of x')
	if remove=='last':
		yrs1 = yrs2[:-(np.size(yrs2)%d)]
	else:
		yrs1 = yrs2[(np.size(yrs2)%d):]
	yrs3 = np.reshape(yrs1, (-1,d))
	yrs5 = np.mean(yrs3,axis=1)


	if nd>1:
		x1 = np.swapaxes(x,0,axis) #make axis of interest first dimension
		n1 = np.shape(x1)
		
		if nd == 2:
			if remove=='last':
				x2 = x1[:-(n[0]%d),:]
			else:
				x2 = x1[(n[0]%d):,:]
			y = np.reshape(x2,(-1,d,n1[1]))  #separate into periods of length d
			xa1 = np.nanmean(y,axis=1)
		elif nd==3:
			if remove=='last':
				x2 = x1[:-(n[0]%d),:,:]
			else:
				x2 = x1[(n[0]%d):,:,:]
			y = np.reshape(x2,(-1,d,n1[1],n1[2]))  #separate into periods of length d
			xa1 = np.nanmean(y,axis=1)
			
		xa = np.swapaxes(xa1,0,axis) #put the axes back in the original order
	else:
		if remove=='last':
			x1 = x[:-(n[0]%d)]
		else:
			x1 = x[(n[0]%d):]
		y = np.reshape(x1, (-1,d))  #separate into periods of length d
		xa = np.nanmean(y,axis=1)
		
	return (xa,yrs5)
							
#################################################################################
def shift_latlon(lat1, lon1, axisa=0, axiso=1, add_ext=1):
	
	# returns lats and lons corresponding to the SW corner of the gridbox, with input at centers
	# if input lat/lon are the same size as the data, supply axis on which to shift
	#
	# us add_ext flag to add an extra lat and lon so that plotting will work correctly
	
	import numpy as np
	
	lat = np.copy(lat1); lon = np.copy(lon1)
					
	#only accept 1- and 2-D input for lat and lon
	nla = np.shape(lat)
	nlo = np.shape(lon)
	if np.size(nla)>2 or np.size(nlo)>2:
		raise ValueError('Please supply 1- or 2-D lat/lon arrays.')
	if np.size(nla)!=np.size(nlo):
		raise ValueError('Curious that lat and lon do not have the same number of dimensions.')
	
	
	#make sure lon is 0-359
	if np.any(lon<0):
		lon[lon<0] += 360	
	
	
	if np.size(nla)==1 or nla[0]==1 or nla[1]==1:	#if 1-D lat/lon input
		difo = np.diff(lon)/2
		difo = np.append(difo,difo[-1]) #just use same value for last
		difa = np.diff(lat)/2
		difa = np.append(difa,difa[-1])
	
		lat2 = lat - difa
		lon2 = lon - difo
		
		if add_ext==1:
			lat2 = np.append(lat2,lat2[-1]+2*difa[-1])
			lon2 = np.append(lon2,lon2[-1]+2*difo[-1])
			
	else:     #if 2-D lat/lon input
		difo = np.diff(lon,axis=axiso)/2
		difo1 = np.repeat(difo,[1]*(nlo[axiso]-2)+[2],axis=axiso)

		difa = np.diff(lat,axis=axisa)/2
		difa1 = np.repeat(difa,[1]*(nla[axisa]-2)+[2],axis=axisa)
		
		lat2 = lat - difa1
		lon2 = lon - difo1
		
	
	#make sure lon is 0-359
	if np.any(lon2<0):
		lon2[lon2<0] += 360	
	
	return lat2,lon2

################################################################################
def gridboxarea(lat1,lon1):
	
	# computes a nlat x nlon array of areas (in sq. km) for each grid box
	# the lat/lon arrays should represent the center of the grid box
	#
	# follows the methodology of http://ferret.pmel.noaa.gov/Ferret/faq/averages-integrals-on-the-sphere
	#


	import numpy as np	

	lat = np.copy(lat1); lon = np.copy(lon1)
	
	R = 	6371	#radius of the earth in km
	
	if np.size(np.shape(lat))>1 or np.size(np.shape(lon))>1:
		raise ValueError('inputs must be vectors (i.e., single column or row)')
	
	#make sure all lons are between 0 and 360	
	if np.any(lon<0):
		lon[lon<0] += 360
		ind = np.argsort(lon)
		sflag = 1
		lon = np.sort(lon)
	else:
		sflag=0
	
	#if decending lats, flip for analysis
	if (lat[1]-lat[0])<0:
		dflag=1
		lat = np.flipud(lat)
	else:
		dflag=0
	
	nla = np.size(lat)
	nlo = np.size(lon)
	
	gba = np.zeros((nla,nlo))
	
	#convert lat and lon from degrees to radians
	latr = lat*np.pi/180	
	lonr = lon*np.pi/180
	
	#calculate width of grid boxes
	difo = np.diff(lonr)
	difo = np.append(difo,difo[-1]) #just use same value for last
	difa = np.diff(latr)
	difa = np.append(difa,difa[-1])	
		
	for i in np.arange(nla):
		for j in np.arange(nlo):	
			gba[i,j] = (R**2)*np.cos(latr[i])*difo[j]*2*np.sin(difa[i]/2)
	
	if dflag==1:
		gba = np.flipud(gba)
	
	if sflag==1:
		gba = gba[:,ind]
		
	return gba

###############################################################################
def dist_globe(lat1,lon1,lat2,lon2):
	
	# calculates the distance in m between two lat,lon points on the globe
	# expects inputs in degrees
	
	import numpy as np
	
	#make sure the originals in the calling script are safe
	lon1 = np.copy(lon1)
	lon2 = np.copy(lon2)
	
	#lons on [0, 360]
	if lon1<0:
		lon1 += 360
	if lon2<0:
		lon2 += 360
	
	#convert to radians
	latr1 = lat1*np.pi/180
	latr2 = lat2*np.pi/180
	lonr1 = lon1*np.pi/180
	lonr2 = lon2*np.pi/180

	r = 	6371000	#radius of the earth in m

	delo = np.abs(lonr2-lonr1)
	
	num = np.sqrt((np.cos(latr2)*np.sin(delo))**2+(np.cos(latr1)*np.sin(latr2)-np.sin(latr1)*np.cos(latr2)*np.cos(delo))**2)
	den = np.sin(latr1)*np.sin(latr2)+np.cos(latr1)*np.cos(latr2)*np.cos(delo)
	
	d = r*np.arctan2(num,den)
	
	return d
	
################################################################################
def quantmap(x,y,n=400,pmap='None',o_flag=0):
	
	# adjusts the data in x by matching the quantiles of y
	# x and y should be vectors
	# n indicates the number of quantiles to use
	# will need to apply check for unphysical values on returned x1
	# if o_flag==1, will return pmap list with info for mapping percentiles
	# option to input pmap from previous mapping to be used this time
	
	import numpy as np
	from scipy import stats, linalg
	
	nx = np.shape(x)
	if np.size(nx) == 2:
		if nx[0] or nx[1] != 1:
			raise ValueError('x must be a vector')
	elif np.size(nx) > 2:
		raise ValueError('x must be a vector')
	ny = np.shape(y)
	if np.size(ny) == 2:
		if ny[0] or ny[1] != 1:
			raise ValueError('y must be a vector')
	elif np.size(ny) > 2:
		raise ValueError('y must be a vector')	

	ndx = np.size(x)
	x1 = np.zeros((ndx,))	

	#rank data
	r = stats.mstats.rankdata(np.ma.masked_invalid(x))	
		
	if pmap is 'None':
		#remove nans
		xz = x[~np.isnan(x)]
		yz = y[~np.isnan(y)]	

		#determine percentile values
		p_step = 100./n
		p = np.arange(0,100,p_step)+p_step
		
		psx = np.percentile(xz,p)
		psy = np.percentile(yz,p)
		
		offset = psy - psx
		
		#linear regression for top percentiles
		x2 = p[p>=95]
		o2 = offset[p>=95]
		x3 = np.hstack((np.ones((np.size(x2),1)), x2[:,None]))
		at = np.dot(linalg.inv(np.dot(x3.T,x3)),np.dot(x3.T,o2[:,None]))	
	
		#linear regression for bottom percentiles
		x4 = p[p<=5]
		o3 = offset[p<=5]
		x5 = np.hstack((np.ones((np.size(x4),1)), x4[:,None]))
		ab = np.dot(linalg.inv(np.dot(x5.T,x5)),np.dot(x5.T,o3[:,None]))
	else:
		offset = pmap[0]
		p = pmap[1]
		ab = pmap[2]
		at = pmap[3]
		
	
	for t in np.arange(ndx):
		if r[t]==0:
			x1[t] = np.nan
		else:
			pn = (100./ndx) * (r[t] - 0.5)
			
			if pn <= 0.5:
				os = ab[1]*pn + ab[0]
			elif pn > 99:
				os = at[1]*pn + at[0]
			else:
				ind = np.argmin(np.abs(p - pn))
				os = offset[ind]
			
			x1[t] = x[t] + os
	
	if o_flag==0:
		return x1
	else:
		pmap = [offset, p, ab, at]
		return (x1,pmap)

################################################################################
def quantmap_del(xf,xh,pmap):
	
	# applies quantile mapping to the time series of xf (future) using the fit
	# pmap from function quantmap above and maintaining the delta with xh (hist)
	
	import numpy as np
	from scipy import stats
	
	nx = np.shape(xf)
	if np.size(nx) == 2:
		if nx[0] or nx[1] != 1:
			raise ValueError('x must be a vector')
	elif np.size(nx) > 2:
		raise ValueError('x must be a vector')


	ndx = np.size(xf)
	xf1 = np.zeros((ndx,))	

	#rank data
	rf = stats.mstats.rankdata(np.ma.masked_invalid(xf))	

	offset = pmap[0]
	p = pmap[1]
	ab = pmap[2]
	at = pmap[3]
	
	psf = np.percentile(xf,p)
	psh = np.percentile(xh,p)
	delt = psf/psh
	
	for t in np.arange(ndx):
		if rf[t]==0:
			xf1[t] = np.nan
		else:
			pn = (100./ndx) * (rf[t] - 0.5)
			
			if pn <= 0.5:
				os = ab[1]*pn + ab[0]
			elif pn > 99:
				os = at[1]*pn + at[0]
			else:
				ind = np.argmin(np.abs(p - pn))
				os = offset[ind]
			
			xf1[t] = (xf[t] + os)*delt[p==pn]
	
	return xf1


###############################################################################
def detrend(x,time,period='None',axis=0):
	
	#returns x with a linear trend removed
	#the time array should match the dimension of x along which the trend is 
	#    calculated (axis)
	#to calculate the trend over a subset of the time series, input an array
	#    with the starting and ending points (period)
	#    the default is to use the entire time period
	
	import numpy as np
	import statfunc_mk as sf
	
	n = np.shape(x)	
	nd = np.size(n)
	
	if n[axis] != np.size(time):
		raise ValueError('time should match the axis of interest in x')
	if nd > 3:
		raise ValueError('Sorry, this function can only handle up to 3 dimensional input right now.')
	if period is not 'None':
		if np.size(period)!=2:
			raise ValueError('period must have 2 entries')
		
	time = np.squeeze(time)	#remove any extra dimensions
	
	
	if nd>1:
		x1 = np.swapaxes(x,0,axis) #make axis of interest first dimension
		n1 = np.shape(x1)
		
		if nd == 2:
			xa1 = np.zeros(n1)
			for i in np.arange(n1[1]):
				x2 = x1[:,i]
				if period is 'None':
					beta1 = sf.regr(time,x2,add_int=1)
				else:
					beta1 = sf.regr(time[(time>=period[0])&(time<=period[1])],
							x2[(time>=period[0])&(time<=period[1]),...],add_int=1)
				xa1[:,i] = x2 - (beta1[0]+beta1[1]*time)

		elif nd==3:
			xa1 = np.zeros(n1)
			for i in np.arange(n1[1]):
				for j in np.arange(n1[2]):
					x2 = x1[:,i,j]
					if period is 'None':
						beta1 = sf.regr(time,x2,add_int=1)
					else:
						beta1 = sf.regr(time[(time>=period[0])&(time<=period[1])],
								x2[(time>=period[0])&(time<=period[1]),...],add_int=1)
					xa1[:,i,j] = x2 - (beta1[0]+beta1[1]*time)
				
		xa = np.swapaxes(xa1,0,axis) #put the axes back in the original order
	else:
		if period is 'None':
			beta1 = sf.regr(time,x,add_int=1)
		else:
			beta1 = sf.regr(time[(time>=period[0])&(time<=period[1])],
					x[(time>=period[0])&(time<=period[1]),...],add_int=1)
		xa = x - (beta1[0]+beta1[1]*time)
	
	return xa
	

###############################################################################
def growdegday(ta,tn,th,end_dates='None',oflag=1):
	
	# calculates the growing degree days for the given threshold th
	# ta and tn are time series of daily mean temperature and daily minimum
	#    temperature for a single year
	# end_dates is an optional array of Julian dates at which to report the GDD
	#    default is for the whole year growing season
	# will return Julian dates of the start and end of accumulation of GDD
	#    unless oflag is set to 0

	import numpy as np
	
	if np.shape(ta)!=np.shape(tn):
		raise ValueError('inputs ta and tn must be the same shape')
	if np.size(np.shape(np.squeeze(ta)))!=1:
		raise ValueError('inputs ta and tn must be single time series')
	if np.size(ta)!=365:
		if np.size(ta)!=366:
			raise ValueError('function only accepts one year at a time')
	
	nd = np.arange(np.size(ta))	
	
	#start 10 days after first over threshold				
	ind1 = np.extract(ta>th,nd)
	if np.size(ind1)==0:	#if no days over threshold, return NaNs
		gdd = np.nan; sd = np.nan; ed = np.nan;
		return (gdd,sd,ed)
	else:
		ta[:ind1[0]+10] = 0
		sd = ind1[0]+10+1
	
	#end Oct 31 (jd=304)
	ta[(304-1):] = 0
	#or after first fall frost (first day with tmin < 0 after Aug. 1)
	ind2 = np.extract(tn[213:]<0.0,nd[213:])
	if np.size(ind2)!=0:
		ta[ind2[0]:] = 0
		ed = ind2[0]+1
		if ed > 304:
			ed = 304
	else:
		ed = 304			
	
	#determine contribution from each day
	ta = ta - th
	#any less than 0 don't count (temps below threshold)
	ta[ta<0] = 0

	if end_dates is 'None':
		gdd = np.sum(ta)
	else:
		nd = np.size(end_dates)
		gdd = np.zeros((nd,))
		for i in np.arange(nd):
			gdd[i] = np.sum(ta[:end_dates[i]])
	
	return (gdd,sd,ed)

###############################################################################
def calc_streaks(x,th,direc='geq',out_time=1):
	
	#calculates the number and length of streaks in the data x that exceed the
	# provided threshold th
	# direction defaults to 'geq' for greater than or equal to th but also
	# accepts 'g', 'l', or 'leq'
	#x must be a single time series
	#set out_time to 0 if the timing of the streaks is not desired
	#assumes NaNs do not meet the criteria
	
	import numpy as np
	
	nx = np.shape(x)
	if np.size(nx) == 2:
		if nx[0] or nx[1] != 1:
			raise ValueError('x must be a vector')
	elif np.size(nx) > 2:
		raise ValueError('x must be a vector')
	
	x = np.squeeze(x)
	
	y = np.zeros(len(x))
		
	if direc=='geq':
		y[x>=th] = 1
	elif direc=='g':
		y[x>th] = 1
	elif direc=='leq':
		y[x<=th] = 1
	elif direc=='l':
		y[x<th] = 1
	else:
		raise ValueError('please provide a valid direction')
		

	y1 = np.diff(y)
	nd1 = np.arange(len(y1),dtype='int')

	ind1 = np.extract(y1==1,nd1)	#start yes
	ind2 = np.extract(y1==-1,nd1)	#start no
	
	if y[~np.isnan(y)][0]==1:	#if starts yes
		ind1 = np.hstack((np.argwhere(~np.isnan(y))[0],ind1))
	if y[~np.isnan(y)][-1]==1:	#if ends yes
		ind2 = np.hstack((ind2,np.argwhere(~np.isnan(y))[-1]))
	if len(ind1)>len(ind2) and ind1[-1]>ind2[-1]:
		ind1 = ind1[:-1]
	if np.any(ind2-ind1)<0:
		raise ValueError('oops...negative durations')

	streaks = ind2 - ind1
	if y[~np.isnan(y)][0]==1:
		streaks[0]+=1
	if y[~np.isnan(y)][-1]==1:
		streaks[-1]+=1
	
	if out_time==1:
		return streaks,ind1+1
	else:
		return streaks


###############################################################################