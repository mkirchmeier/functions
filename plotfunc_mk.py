def plot_borders(map,ax,lw=1):
	from matplotlib.patches import Polygon
	
	file4 = '/HOME/mky/Research/Data/Canada_shape/Coastline/NOAA/GSHHS_shp/l/GSHHS_l_L1'
	file5 = '/HOME/mky/Research/Data/Canada_shape/Coastline/NOAA/WDBII_shp/l/WDBII_border_l_L1'
	file5a = '/HOME/mky/Research/Data/Canada_shape/Coastline/NOAA/WDBII_shp/l/WDBII_border_l_L2'
	file6 = '/HOME/mky/Research/Data/Canada_shape/Coastline/NOAA/GSHHS_shp/l/GSHHS_l_L2'

	map.readshapefile(file4, name='coast', drawbounds=True)
	map.readshapefile(file5, name='countries', drawbounds=True)
	map.readshapefile(file5a, name='states', drawbounds=True)
	map.readshapefile(file6, name='lake',drawbounds=False)

	for i in xrange(50):
		seg = map.lake[i]
		poly = Polygon(seg,facecolor='None',edgecolor='k',linewidth=lw)
		ax.add_patch(poly)
	
###############################################################################
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):

	#from http://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
	
	import matplotlib as mpl
	import numpy as np
	
	new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
	return new_cmap
