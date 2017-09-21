#################################################################################
# A set of functions to convert between meteorological variables
#################################################################################
def slp_pres(slp,tair,elev):

	#slp in hPa, kPa, Pa, etc.
	#tair in C
	#elev in m
	
	import numpy as np
	
	if np.shape(tair)!=np.shape(elev):
		elev = np.squeeze(np.tile(elev,(np.shape(tair)[0],1,1)))	#making assumptions about shapes here
	
	pres = slp * np.exp(-elev / (29.263*(tair+273.15)))
	
	return pres

###############################################################################
def shum_relh(shum,tair,pres,punit='Pa'):
	
	#shum in kg/kg
	#tair in C
	#pres in supplied units
	
	if punit is 'hPa':
		e0 = 6.11
	elif punit is 'kPa':
		e0 = 0.611
	elif punit is 'Pa':
		e0 = 611
	else:
		raise ValueError('please supply valid pressure units')
	
	tair = tair + 273.15	
	
	import numpy as np	
	
	#constants
	ep = 0.622
	T0 = 273.15
	Rv = 461.5
	L  = 2.5*10**6	#assume liquid water
	
	es = e0 * np.exp((L/Rv)*((1/T0)-(1/tair)))	
	
	e = (pres*shum) / (ep + shum*(1-ep))
	
	relh = 100*e/es
	
	return relh
	
###############################################################################

def t_es(tair,punit='Pa'):
	
	#tair in C
	
	import numpy as np
	
	tair = tair + 273.15
	
	if punit is 'hPa':
		e0 = 6.11
	elif punit is 'kPa':
		e0 = 0.611
	elif punit is 'Pa':
		e0 = 611
	else:
		raise ValueError('please supply valid pressure units')
	
	#constants
	T0 = 273.15
	Rv = 461.5
	L  = 2.5*10**6	#assume liquid water
	
	es = e0 * np.exp((L/Rv)*((1/T0)-(1/tair)))	
	
	return es

###############################################################################

def vpd(moist,tair,mvar,pres=None,punit='Pa'):
	
	#tair in C
	#specify moisture variable provided using mvar
	#   current options: relh (relative humidity), huss (specific humidity), mr (mixing ratio)
	#   may need to supply pres
	#relh in % or huss/mr in kg/kg
	
	tair = tair + 273.15	
	
	import numpy as np	
	
	#constants
	T0 = 273.15
	Rv = 461.5
	L  = 2.5*10**6	#assume liquid water
	
	if punit is 'hPa':
		e0 = 6.11
	elif punit is 'kPa':
		e0 = 0.611
	elif punit is 'Pa':
		e0 = 611
	else:
		raise ValueError('please supply valid pressure units')
		
	es = e0 * np.exp((L/Rv)*((1/T0)-(1/tair)))	
	
	if mvar=='relh':
		e = es*moist/100
	elif mvar=='huss':
		ep = 0.622
		if pres==None:
			raise ValueError('pressure input needed for this moisture variable')
		e = (moist*pres) / (moist+ep*(1-moist))
	elif mvar=='mr':
		ep = 0.622
		if pres==None:
			raise ValueError('pressure input needed for this moisture variable')
		e = (moist*pres) / (moist+ep)
	else:
		raise ValueError('please provide approriate moisture variable')
	
	vpd = es - e
	
	return vpd
	
