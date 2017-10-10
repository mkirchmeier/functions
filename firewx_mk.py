###############################################################################
# A set of functions for calculating the FWI System indices from the Canadian
# Forest Fire Danger Rating System
# See http://cwfis.cfs.nrcan.gc.ca/background/summary/fwi
###############################################################################
def FFMC(tair,relh,wspd,pr,restart='None',fwi_flag=0):

	# returns the Fine Fuel Moisture Code (DMC) based on van Wagner and 
	#	Pickett (1985)
	
	import numpy as np
	
	#check dimensions
	if np.any([np.shape(x)!=np.shape(tair) for x in [relh,pr]]):
		raise ValueError('all major input variables must be the same shape')
		
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if restart is 'None':
		restart = np.zeros(np.shape(tair))
		restart[0] = 1	#assume all from same season
	elif np.size(restart)!=np.size(tair):
		raise ValueError('restart array must be same size as mon array')
		
	wspd = wspd*3.6	#from m/s to km/h
	
	#if relh > 100 set to 100 --> added to avoid nan results
	relh[relh>100] = 100.	
	
	n = np.size(tair)
	ffmc = np.zeros((n,))
	ms = np.zeros((n,))
	
	for i in np.arange(n):
		if restart[i]==1:
			f0 = 85
		else:
			f0 = ffmc[i-1]
			k = 1
			while np.isnan(f0):
				f0 = ffmc[i-1-k]	#if previous value is NaN, take last available
				k += 1
				if k>5:		#if going back more than 5 days, reset to default
					f0 = 85
			
		if np.isnan(tair[i]) or np.isnan(pr[i]) or np.isnan(relh[i]) or np.isnan(wspd[i]):
			ffmc[i] = np.nan 		#if missing input, return nan
			ms[i] = np.nan			
			continue
		
		m0 = 147.2 * (101 - f0) / (59.5 + f0)
		
		if pr[i] > 0.5:
			rf = pr[i] - 0.5
			if m0 <= 150:
				mr = m0 + 42.5 * rf * np.exp(-100/(251-m0))*(1-np.exp(-6.93/rf))
			else:
				mr = m0 + 42.5 * rf * np.exp(-100/(251-m0))*(1-np.exp(-6.93/rf)) + 0.0015 * (m0 - 150)**2 * rf**0.5
			if mr > 250:
				mr = 250
			m0 = mr
		
		ed = 0.942*relh[i]**0.679 + 11*np.exp((relh[i]-100)/10) + 0.18*(21.1 - tair[i])*(1-np.exp(-0.115*relh[i]))
		
		if m0 > ed:
			k0 = 0.424 * (1 - (relh[i]/100)**1.7) + 0.0694 * wspd[i]**0.5 * (1 - (relh[i]/100)**8)
			kd = k0 * 0.581 * np.exp(0.0365*tair[i])
			m = ed + (m0 - ed)*10**-kd
		else:
			ew = 0.618*relh[i]**0.753 + 10*np.exp((relh[i]-100)/10) + 0.18*(21.1 - tair[i])*(1-np.exp(-0.115*relh[i]))
			if m0 < ew:
				kl = 0.424*(1 - ((100-relh[i])/100)**1.7) + 0.694 * wspd[i]**0.5 * (1 - ((100-relh[i])/100)**8)
				kw = kl * 0.581*np.exp(0.0365*tair[i])
				m = ew - (ew - m0)*10**-kw
			else:
				m = m0
		
		f = 59.5 * (250 - m)/(147.2 + m)
		
		ffmc[i] = f
		ms[i] = m
		
	if fwi_flag==1:
		return ms
	else:
		return ffmc
		


###############################################################################
def DMC(tair,relh,prcp,mon,dmc0='None',restart='None'):
	
	# returns the Duff Moisture Code (DMC) based on van Wagner and Pickett (1985)
	
	import numpy as np
	
	#check dimensions
	if np.any([np.shape(x)!=np.shape(tair) for x in [relh,prcp]]):
		raise ValueError('all major input variables must be the same shape')
		
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(mon)!=np.size(tair):
		raise ValueError('size of mon array must match size of wx variables')
		
	if restart is 'None':
		restart = np.zeros(np.shape(tair))
		restart[0] = 1	#assume all from same season
	elif np.size(restart)!=np.size(mon):
		raise ValueError('restart array must be same size as mon array')
	
	if dmc0 is 'None':
		dmc0 = 6
	else:
		dmc0 = dmc0[~np.isnan(dmc0)]
	if np.size(dmc0)!=np.sum(restart):
		if np.size(dmc0)==1:
			dmc0 = np.repeat(dmc0,np.sum(restart))
		else:
			raise ValueError('dmc0 must either be a single value or one-per year')
	if np.sum(restart)==1:
		dmc0 = np.array([dmc0])
			
	les = np.array([6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0])	
	
	mon = mon.astype('int')
	le = les[mon-1]	
	
	n = np.size(tair)
	dmc = np.zeros((n,))
	
	j = 0
	for i in np.arange(n):

		if restart[i]==1:
			p0 = dmc0[j]
			j = j+1
		else:
			p0 = dmc[i-1]
			q = 1
			while np.isnan(p0):
				p0 = dmc[i-1-q]	#if previous value is NaN, take last available
				q += 1
				if q>15:		#if going back more than 15 days, reset to default
					p0 = dmc0[j-1]
			
		if np.isnan(tair[i]) or np.isnan(prcp[i]) or np.isnan(relh[i]):
			dmc[i] = np.nan 		#if missing input, return nan
			continue			

		if prcp[i] > 1.5:		
			re = 0.92*prcp[i] - 1.27
		
			mo = 20 + np.exp(5.6348 - p0/43.43)
			
			if p0 <= 33:
				b = 100 / (0.5 + 0.3*p0)
			elif p0 <= 65:
				b = 14 - 1.3*np.log(p0)
			else:
				b = 6.2*np.log(p0) - 17.2
			
			mr = mo + 1000*re / (48.77 + b*re)
			
			pr = 244.72 - 43.43*np.log(mr-20)
			if pr < 0:
				pr = 0
			
		else:
			pr = p0
		
		if tair[i] < -1.1:
			t = -1.1
		else:
			t = tair[i]
			
		k = 1.894 * (t + 1.1)*(100-relh[i])*le[i]*10**-6
		
		p = pr + 100*k
		
		dmc[i] = p
		
	return dmc
	
###############################################################################
def DC2(tair,pr,mon,dc0='None',restart='None'):
	
	# returns the Drought code (DC) based on van Wagner and Pickett (1985)
	
	import numpy as np
	
	#check dimensions
	if np.shape(tair) != np.shape(pr):
		raise ValueError('all major input variables must be the same shape')
		
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
		
	if np.size(mon)!=np.size(tair):
		raise ValueError('size of mon array must match size of wx variables')
	
	if restart is 'None':
		restart = np.zeros(np.shape(tair))
		restart[0] = 1	#assume all from same season
	elif np.size(restart)!=np.size(mon):
		raise ValueError('restart array must be same size as mon array')
	
	if dc0 is 'None':
		dc0 = 15
	else:
		dc0 = dc0[~np.isnan(dc0)]
	if np.size(dc0)!=np.sum(restart):
		if np.size(dc0)==1:
			dc0 = np.repeat(dc0,np.sum(restart))
		else:
			raise ValueError('dc0 must either be a single value or one per year')
	if np.sum(restart)==1:
		dc0 = np.array([dc0])	
		
	lfs = np.array([-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6])	
	
	mon = mon.astype('int')
	lf = lfs[mon-1]	
	
	n = np.size(tair)
	dc = np.zeros((n,))
	
	j = 0
	for i in np.arange(n):
		if restart[i]==1:
			d0 = dc0[j]
			j = j+1
		else:
			d0 = dc[i-1]
			k = 1
			while np.isnan(d0):
				d0 = dc[i-1-k]	#if previous value is NaN, take last available
				k += 1
				if k>60:		#if going back more than 60 days, reset to default
					d0 = dc0[j-1]

		if np.isnan(tair[i]) or np.isnan(pr[i]):
			dc[i] = np.nan 		#if missing input, return nan
			continue

		if pr[i]>2.8:
			rd = 0.83*pr[i] - 1.27
			
			q0 = 800 * np.exp(-d0/400)
			
			qr = q0 + 3.937*rd
			
			dr = 400 * np.log(800/qr)
			if dr < 0:
				dr = 0
		
		else:
			dr = d0
		
		if tair[i] < -2.8:
			t = -2.8
		else:
			t = tair[i]
			
		v = 0.36*(t+2.8) + lf[i]
		if v < 0:
			v = 0
		
		d = dr + 0.5*v		
		
		dc[i] = d
		
	return dc
	
###############################################################################
def FWI(tair,relh,wspd,pr,mon,dc0='None',dmc0='None',restart='None',return_flag='FWI'):
	
	# returns the Fire Weather Index based on van Wagner and Pickett (1985)
	# for daily data
	#
	# inputs:
	#	tair in C
	#	relh in %
	#	wspd in m/s
	#	pr in mm/day
	#
	# only accepts input variables that are single values or time series
	#	if data contains more dimensions, please call FWI within loop
	# the restart input is a variable of the same length as the weather inputs
	#	that contains 1s to indicate a new fire season and 0s everywhere
	#	else; when restart==1, previous values of various indices are not
	#	carried over
	# relies on variables calculated with other fire indices

	import numpy as np
	
	#check dimensions
	if np.any([np.shape(x)!=np.shape(tair) for x in [wspd,relh,pr]]):
		raise ValueError('all major input variables must be the same shape')
		
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(mon)!=np.size(tair):
		raise ValueError('size of mon array must match size of wx variables')
	
	#convert units
	wspd = wspd*3.6	#from m/s to km/h
	
	if restart is 'None':
		restart = np.zeros(np.shape(tair))
		restart[0] = 1	#assume all from same season
	elif np.size(restart)!=np.size(mon):
		raise ValueError('restart array must be same size as mon array')
	
	#variables needed from other indices
	m = FFMC(tair,relh,wspd/3.6,pr,restart,fwi_flag=1)
	p = DMC(tair,relh,pr,mon,dmc0=dmc0,restart=restart)
	d = DC2(tair,pr,mon,dc0=dc0,restart=restart)
	
	##FWI
	
	fw = np.exp(0.05039*wspd)
	ff = 91.9*np.exp(-0.1386*m)*(1.0 + (m**5.31)/(4.93*10**7))
	
	r = 0.208 * fw * ff
	
	if return_flag is 'isi' or return_flag is 'ISI':
		return r
	
	u = np.zeros(np.shape(p))
	ind = np.zeros(np.shape(p))	
	ind[p <= 0.4*d] = 1
	ind = ind.astype('bool')
	u[ind] = 0.8 * p[ind]*d[ind] / (p[ind] + 0.4*d[ind])
	u[ind==False] = p[ind==False] - (1 - 0.8*d[ind==False]/(p[ind==False] + 
				0.4*d[ind==False]))*(0.92 + (0.0114*p[ind==False])**1.7)
	
	if return_flag is 'bui' or return_flag is 'BUI':
		return u
	
	if return_flag is 'ib':
		return [r,u]
	
	fd = np.zeros(np.shape(u))
	fd[u<=80.] = 0.626*(u[u<=80.]**0.809) + 2.0
	fd[u>80.] = 1000. / (25. + 108.64*np.exp(-0.023*u[u>80.]))
	
	b = 0.1 * r * fd
	
	s = np.zeros(np.shape(b))
	s[b<=1.] = b[b<=1.]
	s[b>1.] = np.exp(2.72*(0.434*np.log(b[b>1.]))**0.647)

	fwi = s

	if return_flag is 'dict' or return_flag is 'd':
		rd = {}
		rd['fwi'] = fwi
		rd['bui'] = u
		rd['isi'] = r
		rd['dc'] = d
		rd['dmc'] = p
		rd['ffmc'] = 59.5 * (250 - m)/(147.2 + m)
		
		return rd

	return fwi

###############################################################################
def DSR(fwi):

	dsr = 0.0272 * (fwi**1.77)

	return dsr

###############################################################################	
def fire_season_length(tair,dmy,ths=[12,5],end_after='None'):
	
	# calculates the length of the fire season for each year based on the 
	#  definition in Wotton and Flannigan (1993)
	#  default thresholds use 12 C for the start date and 5 C for the end date,
	#  but alternate thresholds can be supplied
	#  the optional input end_after is a Julian date after which to start looking
	#     for the fire-season-end criteria; default is the start date
	# tair should be a single time series of full-year daily values in C
	
	import numpy as np
	
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(dmy[:,0])!=np.size(tair):
		raise ValueError('length of dmy array must match size of tair')	
	
	yrs = np.unique(dmy[:,2])
	ny = np.size(yrs)
	
	sd = np.zeros((ny,))
	ed = np.zeros((ny,))
	
	for i in np.arange(ny):
		
		tairy = tair[dmy[:,2]==yrs[i]]
		
		nd = np.arange(np.size(tairy))+1
		xl1 = np.nan
		xl2 = np.nan
		
		#start date
		ts1 = np.zeros(np.shape(tairy))
		ts1[tairy>=ths[0]] = 1
		ts2 = np.diff(ts1,axis=0)

		ind1 = np.extract(ts2==1,nd)
		ind2 = np.extract(ts2==-1,nd)
		if np.size(ind1)==0:
			sd[i] = np.nan
			ed[i] = np.nan
			continue
		ind1 = ind1[ind1!=364]	#if last day starts run, remove
		ind2 = ind2[ind2!=1]		#if first day ends run, remove
		if ind2[0]<ind1[0]:
			ind1 = np.hstack((1,ind1))	#if year starts with temps above, add start date at 01 Jan
		if np.size(ind1)>np.size(ind2):
			xl1 = ind1[-1]
			ind1 = ind1[:-1]	#if year ends with temps above, remove final
		ni = ind2 - ind1
		if np.any(ni<0):
			raise ValueError('year starts with fire season')
		ns = nd[ind1[ni>=3]]
		if np.size(ns)==0:
			sd[i] = xl1+3
		else:
			sd[i] = ns[0]+2

		#end date
		if np.isnan(sd[i]):
			ed[i] = np.nan     #if start date does not exist, neither will end date
			continue
		
		te1 = np.zeros(np.shape(tairy))
		if end_after is 'None':
			ea = sd[i]
		else:
			ea = np.max((end_after,sd[i]))
		te1[(tairy<=ths[1])&(nd>ea)] = 1
		te2 = np.diff(te1,axis=0)

		ind3 = np.extract(te2==1,nd)
		ind4 = np.extract(te2==-1,nd)
		if np.size(ind3)>np.size(ind4):
			xl2 = ind3[-1]
			ind3 = ind3[:-1]	#remove final value; year ends with streak of temps below threshold
		ni2 = ind4 - ind3
		ne = nd[ind3[ni2>=3]]
		if np.size(ne)==0:
			if np.isnan(xl2) | (xl2+3 > nd[-1]):
				ed[i] = nd[-1]
			else:
				ed[i] = xl2+3
		else:
			ed[i] = ne[0]+2
	
	fsl = ed - sd + 1
	
	return (fsl,sd,ed)
	
###############################################################################
def start_stop(tair,dmy,ths=[6,6],r_flag=0):
	
	#  determines the dates to start and stop the fire weather calculations
	#  for each year; returns indices of data to use
	#  if r_flag==1, returns the corresponding restart array
	
	import numpy as np
	
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(dmy[:,0])!=np.size(tair):
		raise ValueError('length of dmy array must match size of tair')	

	(fsl,sd,ed) = fire_season_length(tair,dmy,ths=ths,end_after=182)

	yrs = np.unique(dmy[:,2])
	ny = np.size(yrs)
	
	n = np.size(tair)	
	
	ind = np.zeros((n,))
	restart = np.zeros((n,))
	for j in np.arange(ny):
		if np.isnan(sd[j]):
			continue
		i1 = np.extract(dmy[:,2]==yrs[j],np.arange(n))[0]
		ind[i1+sd[j]-1:i1+ed[j]] = 1
		restart[i1+sd[j]-1] = 1
	ind = ind.astype('bool')
	
	restart=restart[ind]	
	
	if r_flag==1:
		return (ind,restart)
	else:
		return ind

###############################################################################
def start_stop2(tair,pr,snow,dmy,ths=[6,6],jdsnow=60,r_flag=0):
	
	#  determines the dates to start and stop the fire weather calculations
	#  for each year; returns indices of data to use and inital values for
	#  DC and DMC
	#  ths includes the temperatures for the start/stop calculation if snow
	#  is not present; jdsnow is the Julian day after which to choose a snow
	#  start date
	#  if r_flag==1, returns the corresponding restart array
	#  see http://cwfis.cfs.nrcan.gc.ca/background/dsm/fwi for method
	#  ignoring overwintering for now
	#  expects snow in cm, consider snow free if < 1 cm
	
	import numpy as np
	
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(dmy[:,0])!=np.size(tair):
		raise ValueError('length of dmy array must match size of tair')	

	(fsl,sd,ed) = fire_season_length(tair,dmy,ths=ths,end_after=182)

	yrs = np.unique(dmy[:,2])
	ny = np.size(yrs)
	
	n = np.size(tair)	
	
	ind = np.zeros((n,))
	dmc0 = np.zeros((ny,)); dmc0[:] = np.nan
	dc0 = np.zeros((ny,)); dc0[:] = np.nan
	restart = np.zeros((n,))
	
	for j in np.arange(ny):
		#if too cold, then skip
		if np.isnan(sd[j]):
			continue			
		#determine if sig. snow cover
		snowyr1 = snow[(dmy[:,2]==yrs[j])&(dmy[:,1]<=2)]
		if np.sum(np.isnan(snowyr1))/float(np.size(snowyr1))>=.10:		#skip if missing too much data
			continue
		mdepth = np.mean(snowyr1)
		scover = np.size(snowyr1[snowyr1>0.1])/float(np.size(snowyr1))
		if mdepth>10 and scover>=.75:
			#determine first instance of snow free for 3 days
			snowyr = snow[(dmy[:,2]==yrs[j])]
			if np.sum(np.isnan(snowyr))/float(np.size(snowyr))>=.10:		#skip if missing too much data
				continue
			snowyr2 = snowyr[jdsnow-1:]	#start after 01 March or supplied date
			if np.all(snowyr2>=1):	#if snow cover never/barely leaves, continue to next year
				continue
			elif np.size(snowyr[snowyr<1])/float(np.size(snowyr))<=.10:
				continue
			
			ns = np.size(snowyr2)
			nd = np.arange(ns)+jdsnow
			nd1 = np.arange(ns)
			ts1 = np.zeros((ns))
			ts1[snowyr2<1] = 1
			ts2 = np.diff(ts1)
	
			ind1 = np.extract(ts2==1,nd1)	#start snow-free periods
			ind2 = np.extract(ts2==-1,nd1)	#start snow-cover periods
			
			if ts1[0]==1:	#if starts snow-free
				ind1 = np.hstack((np.array([0]),ind1))
			if ts1[-1]==1:	#if ends snow-free
				ind2 = np.hstack((ind2,np.array([nd1[-1]])))
			if np.size(ind1)>np.size(ind2) and ind1[-1]>ind2[-1]:
				ind1 = ind1[:-1]
			if np.any(ind2-ind1)<0:
				raise ValueError('oops...negative durations')
			
			ni = ind2-ind1
			sd1 = nd[ind1[ni>=3][0]]+4	   #start after 3 days of no snow (after 01 March or supplied date)		
			
			dmc0[j] = 6
			dc0[j] = 15
			
			#determine return of snow cover (use as end date if earlier than temp-based date)
			sfm = nd[ind1[np.argmax(ni)]]	#determine start of max snow-free period
			edp = nd[ind2[ind2+jdsnow>np.max((sfm,sd1))][0]]  #determine when snow returns after start date and max snow-free period
#			if sd1 > ed[j]:
#				ed1 = edp
#			else:
#				ed1 = np.min((edp,ed[j]))
			ed1 = edp
			
		else:
			if np.isnan(sd[j]):
				continue
			sd1 = sd[j]
			ed1 = ed[j]

			#determine # of days since precip
			pryr = pr[(dmy[:,2]==yrs[j])]
			pd = 0; w = 0
			while pd==0:
				pd = pryr[int(sd1)-1-w]
				w = w+1
			dmc0[j] = 2*(w-1)
			dc0[j] = 5*(w-1)

		i1 = np.extract(dmy[:,2]==yrs[j],np.arange(n))[0]
		ind[i1+sd1-1:i1+ed1] = 1
		restart[i1+sd1-1] = 1
		
	ind = ind.astype('bool')
	
	restart=restart[ind]
	
	dc0 = dc0[~np.isnan(dc0)]
	dmc0 = dmc0[~np.isnan(dmc0)]
	
	if r_flag==1:
		return (ind,dmc0,dc0,restart)
	else:
		return ind,dmc0,dc0	
	
###############################################################################
def start_stop3(tair,pr,snow,dmy,ths=[6,6],jdsnow=60,r_flag=0):
	
	#  same as start_stop2 but snow is a binary variable
	#  determines the dates to start and stop the fire weather calculations
	#  for each year; returns indices of data to use and inital values for
	#  DC and DMC
	#  ths includes the temperatures for the start/stop calculation if snow
	#  is not present; jdsnow is the Julian day after which to choose a snow
	#  start date
	#  if r_flag==1, returns the corresponding restart array
	#  see http://cwfis.cfs.nrcan.gc.ca/background/dsm/fwi for method
	#  ignoring overwintering for now
	
	import numpy as np
	
	if int(np.squeeze(np.shape(np.squeeze(tair)))) != np.size(tair):
		raise ValueError('only input for a single time series is supported')
	
	if np.size(dmy[:,0])!=np.size(tair):
		raise ValueError('length of dmy array must match size of tair')	

	(fsl,sd,ed) = fire_season_length(tair,dmy,ths=ths,end_after=182)

	yrs = np.unique(dmy[:,2])
	ny = np.size(yrs)
	
	n = np.size(tair)	
	
	ind = np.zeros((n,))
	dmc0 = np.zeros((ny,)); dmc0[:] = np.nan
	dc0 = np.zeros((ny,)); dc0[:] = np.nan
	restart = np.zeros((n,))
	
	for j in np.arange(ny):
		#if too cold, then skip
		if np.isnan(sd[j]):
			continue			
		#determine if sig. snow cover
		snowyr1 = snow[(dmy[:,2]==yrs[j])&(dmy[:,1]<=2)]
		if np.sum(np.isnan(snowyr1))/float(np.size(snowyr1))>=.10:		#skip if missing too much data
			continue
		scover = np.size(snowyr1[snowyr1>0])/float(np.size(snowyr1))
		if scover>=.75:
			#determine first instance of snow free for 3 days
			snowyr = snow[(dmy[:,2]==yrs[j])]
			if np.sum(np.isnan(snowyr))/float(np.size(snowyr))>=.10:		#skip if missing too much data
				continue
			snowyr2 = snowyr[jdsnow-1:]	#start after 01 March or supplied date
			if np.all(snowyr2>0):	#if snow cover never/barely leaves, continue to next year
				continue
			elif np.size(snowyr[snowyr<0.1])/float(np.size(snowyr))<=.10:
				continue
			
			ns = np.size(snowyr2)
			nd = np.arange(ns)+jdsnow
			nd1 = np.arange(ns)
			ts1 = np.zeros((ns))
			ts1[snowyr2==0] = 1
			ts2 = np.diff(ts1)
	
			ind1 = np.extract(ts2==1,nd1)	#start snow-free periods
			ind2 = np.extract(ts2==-1,nd1)	#start snow-cover periods
			
			if ts1[0]==1:	#if starts snow-free
				ind1 = np.hstack((np.array([0]),ind1))
			if ts1[-1]==1:	#if ends snow-free
				ind2 = np.hstack((ind2,np.array([nd1[-1]])))
			if np.size(ind1)>np.size(ind2) and ind1[-1]>ind2[-1]:
				ind1 = ind1[:-1]
			if np.any(ind2-ind1)<0:
				raise ValueError('oops...negative durations')
			
			ni = ind2-ind1
			if np.all(ni<3):
				continue
			sd1 = nd[ind1[ni>=3][0]]+4	   #start after 3 days of no snow (after 01 March or supplied date)
			if sd1 > ed[j]:		#if snow cover leaves after temperatre end date, then no fire season
				continue
			
			dmc0[j] = 6
			dc0[j] = 15
			
			#determine return of snow cover (use as end date if earlier than temp-based date)
			sfm = nd[ind1[np.argmax(ni)]]	#determine start of max snow-free period
			edp = nd[ind2[ind2+jdsnow>np.max((sfm,sd1))][0]]+1  #determine when snow returns after start date and max snow-free period
			ed1 = np.min((edp,ed[j]))
			
		else:
			if np.isnan(sd[j]):
				continue
			sd1 = sd[j]
			ed1 = ed[j]

			#determine # of days since precip
			pryr = pr[(dmy[:,2]==yrs[j])]
			pd = 0; w = 0
			while pd==0:
				pd = pryr[int(sd1)-1-w]
				w = w+1
			dmc0[j] = 2*(w-1)
			dc0[j] = 5*(w-1)

		i1 = np.extract(dmy[:,2]==yrs[j],np.arange(n))[0]
		ind[i1+sd1-1:i1+ed1] = 1
		restart[i1+sd1-1] = 1
		
	ind = ind.astype('bool')
	
	restart=restart[ind]
	
	dc0 = dc0[~np.isnan(dc0)]
	dmc0 = dmc0[~np.isnan(dmc0)]
	
	if r_flag==1:
		return (ind,dmc0,dc0,restart)
	else:
		return ind,dmc0,dc0	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
