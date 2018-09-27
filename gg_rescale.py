# ggg_rescale.py
import numpy as np 
import os
import matplotlib.pyplot as plt

def gg_rescale(model_dir, model_file, rs, rp, gp, tp):

	"""
	gg_rescale
	This function will take a single Generic ATMO Grid file and scale it 
	with an input planetary system radii, gravity and temperature. 


	INPUTS:
	dir = string for the directory containing the file you want to scale
		e.g. '/Users/name/Documents/Generic_Goyal_Grid/Rainout_condensation/'
	file = string of the model name you want to run the input file should be selected such that it is closest to the gravity and temperature for the planet
		e.g. trans-iso-generic_1500_05_+0.0_0.56_0001_0.00_model.txt'
	rs = radius of the star in solar radii
	rp = radius of the planet in jupiter radii
	gp = gravity of the planet in cgs units
	tp = temperature of the planet in K

	OUTPUTS:
	wav = wavlength array for the model (microns)
	depth = (Rp/R*)^2 model
	"""


	# SET UP THE CONSTANTS
	kb = 1.380658E-16 # gm*cm^2/s^2 * Kelvin
	mu = 1.6726E-24 * 2.3 #g  cgs  Hydrogen + Helium Atmosphere
	tau = 0.56 # optical depth
	rsun = 69580000000 # cm
	rjup = 6991100000 # cm

	# Read in the file from the directory specified by the user
	# define the columns of data from the file read in 
	model_filename = os.path.join(model_dir,model_file)
	print(model_filename)

	model_wav, model_rp = np.loadtxt(model_filename, unpack=True)

	model_rprs = np.sqrt(model_rp) * (rp/rs)

	rstar = rs * rsun
	rpl = rp * rjup
	rp2 = rpl

	# Extract the starting temperature and gravity from the input filename
	t1 = float(model_file.split('_')[1])
	g1 = float(model_file.split('_')[2])
	g1 = g1*1e2 #rescaled to cm

	# Calculate the transmission spectrum based on the file parameters
	#	 to get the baseline radius
	h1 = (kb * t1) / (mu * g1)
	rp1 = np.sqrt(model_rp) * rsun #cm
	z1 = rp1 - (np.sqrt(model_rp[2000])*rsun) #cm
	epsig1 = tau * np.sqrt((kb * t1 * mu * g1) / (2. * np.pi * rp1)) * np.exp(z1 / h1)

	# Rescale based on the input parameters

	h2 = (kb * tp) / (mu * gp)
	z2 = h2 * np.log(epsig1 / tau * np.sqrt((2. * np.pi * rp2)/(kb * tp * mu * gp)))
	r2 = z2 + rp2

# Resort the data so it goes from lowest wavelength to highest
	srt = np.argsort(model_wav)
	wav = model_wav[srt]
	r2 = r2[srt]

	depth = (r2 / rstar)**2

# These are the OUTPUTS
	return wav, depth














if __name__ == '__main__':
	model_dir = './Data'
	model_file = 'trans-iso-generic_0400_05_+0.0_0.35_0001_0.06_model.txt'
	# model_file = 'trans-iso-generic_1500_05_+0.0_0.56_0001_0.00_model.txt'
	rs = 1.100 # Rstar
	rp = 1.200 # RJ
	gp = 610 # cgs
	tp = 450 # K

	# result = gg_rescale(model_dir, model_file, rs, rp, gp, tp)
	# This will give you a tuple with wavelength and transit depth.
	# Or do the following to put them into your own arrays
	wav, rprs = gg_rescale(model_dir, model_file, rs, rp, gp, tp)


	plt.plot(wav, rprs)
	plt.xlim(0.3,15)
	plt.xscale('log')
	plt.tight_layout()
	plt.show()





