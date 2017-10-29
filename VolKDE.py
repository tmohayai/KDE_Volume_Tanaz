"""
Author: Tanaz A. Mohayai
Copyright (C) 2016-present by Tanaz A. Mohayai. All rights reserved.

This module uses the Kernel Density Estimation Technique from Python's Scipy module, "kde" [1], as defined in its default 
configuration and computes the density and volume associated with a distribution of muons, particles or data points in a 
4-dimensional position-momentum phase space, for a user-defined contour. 

This script, along with others in this directory serve the Muon Ionization Cooling Experiment, MICE and are part of the 
author's PhD thesis work on MICE [2, 3, 4, 5, 6]. 
There is continous progress on the analysis and modules presented here and the author remains the sole contributor. 
However, suggestions for this work are encouraged. This agreement may change without notice. 

The function parameters defined in this module along with their descriptions are as following: 

data: input data array, with each row representing the particles, muons or data points and each column representing the 
coordinates, i.e. x, px, y, py.

sample_size: size of the sample or number of muons - important for finding the user-defined contour. 

percent: percentage of the particles that reside inside the user-defined contour. 

i: the location/s along the MICE channel where the distrbutions are evaluated - this could be the virtual detector number 
in G4-beamline [7]. 

x_colm: column number in the data file corresponding to x position coordinates.

px_colm: column number in the data file corresponding to px momentum coordinates.

y_colm: column number in the data file corresponding to y position coordinates.

py_colm: column number in the data file corresponding to py momentum coordinates.

region_colm: column number in the data file corresponding to region number, representing the z positions or location/s of 
the particles or muons along the MICE channel.

N_mc: number of monte carlo points with default being 1e6. 

Note: the default values are relevant when an ICOOL's "for009.dat" [8] data array format is used.

Note: parameter "percent" has to be in form of a quotient (0.09 for 9% of the sample size). 

Please note that this module, along with others which are currently in author's posession are planned to be merged into 
the official MICE Analysis User Software, MAUS package [9].

[1] Scipy's "gaussian_kde()" module by R. Kern, http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html.
[2] T. A. Mohayai et al., "Novel Application of Non-parametric Density Estimation Technique in Muon Ionization Cooling Experiment," APS-DPF'17 Proceedings.
[3] T. A. Mohayai, "Novel Application of Kernel Density Estimation in MICE," MICE-Note-506.
[4] T. A. Mohayai et al., "Novel Implementation of Non-parametric Density Estimation in MICE," IPAC'17 Proceedings.
[5] T. A. Mohayai et al., "Simulated Measurements of Beam Cooling in Muon Ionization Cooling Experiment," NA-PAC'16 Proceedings.
[6] T. A. Mohayai et al., "Simulated Measurements of Beam Cooling in Muon Ionization Cooling Experiment," NA-PAC'16 Proceedings.
[7] T. Roberts, "G4beamline User's Guide", Muons, Inc (2013), http://www.muonsinc.com. 
[8] http://www.cap.bnl.gov/ICOOL/fernow/readme.html
[9] C. D. Tunnell, C. T. Rogers, "MAUS: MICE Analysis User Software", IPAC (2011).
"""

from pylab import *
import numpy as np
from math import *
import scipy.stats as st


def volkde(data, sample_size, i, percent, x_colm, px_colm, y_colm, py_colm, region_colm, N_mc):

	    # Loading each muon position and momentum
		x = data[data[:,region_colm]==i, x_colm] 
		px = data[data[:,region_colm]==i, px_colm] 
		y = data[data[:,region_colm]==i, y_colm] 
		py = data[data[:,region_colm]==i, py_colm] 
		
		data_array = np.vstack([x,px,y,py])
		
		# Computing the distribution of the four-dimensional muon position and momentum 
		data_kde_object = st.gaussian_kde(data_array)
		data_kde = data_kde_object.evaluate(data_array)
		

		# Binary search below added by Pavel Snopok, psnopok@iit.edu 
		numb_points=round(percent*sample_size)
		density_r=max(data_kde)
		density_l=0
		density=density_l
		count=(data_kde>density).sum()
		while count<>numb_points:
			density=(density_l+density_r)/2
			count=(data_kde>density).sum()
			if count>numb_points:
				density_l=density
			else:
				density_r=density
		# Density of the contour under study
		dens=density
		flag = (data_kde>density) 
		xmin=x[flag].min()
		xmax=x[flag].max()
		ymin=y[flag].min()
		ymax=y[flag].max()
		
		pxmin=px[flag].min()
		pxmax=px[flag].max()
		pymin=py[flag].min()
		pymax=py[flag].max()
		
		hyper_cube = (xmax-xmin)*(ymax-ymin)*(pxmax-pxmin)*(pymax-pymin)
		
		mcx = random(N_mc)*float(xmax-xmin)+xmin
		mcy = random(N_mc)*float(ymax-ymin)+ymin
		mcpx = random(N_mc)*float(pxmax-pxmin)+pxmin
		mcpy = random(N_mc)*float(pymax-pymin)+pymin
		mc_kde = data_kde_object.evaluate(np.vstack([mcx,mcpx,mcy,mcpy]))
		# Volume of the contour under study
		volume = (mc_kde>=dens).sum()*(hyper_cube/N_mc)

		output = np.column_stack((volume, dens, i))
		return np.savetxt('volume_density_region'+str(i)+'.dat', output)
		
if __name__=="__main__":
	"""
	Author: Tanaz A. Mohayai
	Copyright (C) 2016-present by Tanaz A. Mohayai, Illinois Institute of Technology. All rights reserved.
	
	This example demonstrates the use of the VolKDE module. The for009_dummy is a dummy input array 
	and is used for illustrative purposes, only. Please note that this example is also separately created in this directory.
	
	"""

	
	for i in range(1, 4):
		data = loadtxt('for009_dummy.dat', skiprows=2)
		
		x_colm, px_colm, y_colm, py_colm = 6, 9, 7, 10 
		
		percent = 0.09
		
		region_colm = 4
		
		N_mc = 1e3
		
		sample_size, columns= data[data[:,region_colm]==i,:].shape
		
		
	
	
