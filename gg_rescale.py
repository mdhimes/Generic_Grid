# ggg_rescale.py
import numpy as np 
import os
import matplotlib as mpl
# mpl.use("TkAgg") # uncomment if working woth MacXCode in Mac 10.14 or higher
import matplotlib.pyplot as plt

def gg_rescale(model_dir, model_file, radius_star, radius_planet, gravity_planet, temperature_planet, mean_molecular_weight=2.3):

	"""
	gg_rescale
	This function will take a single Generic ATMO Grid file and scale it 
	with an input planetary system radii, gravity and temperature. 


	INPUTS:
	dir = string for the directory containing the file you want to scale
		e.g. '/Users/name/Documents/Generic_Goyal_Grid/Rainout_condensation/'
	file = string of the model name you want to run the input file should be selected such that it is closest to the gravity and temperature for the planet
		e.g. trans-iso-generic_1500_05_+0.0_0.56_0001_0.00_model.txt'
	radius_star = radius of the star in solar radii
	radius_planet = radius of the planet in jupiter radii
	gravity_planet = gravity of the planet in cgs units
	temperature_planet = temperature of the planet in K
	mean_molecular_weight = scale factor for the atmospheric mean molecular weight. This will be multiplied by the mass of a proton and used to calculate the planetary scale height. The default is set to 2.3 for H/He-dominated atmosphere.
	
	OUTPUTS:
	model_wavelength = wavlength array for the model (microns)
	transit_depth = (Rp/R*)^2 model
	"""


	# SET UP THE CONSTANTS
	kb = 1.380658E-16 # gm*cm^2/s^2 * Kelvin
	mu_model = 1.6726E-24 * 2.3 #g  cgs  Hydrogen + Helium Atmosphere
	mu_planet = 1.6726E-24 * mean_molecular_weight #g  cgs rescaled atmospheric mean molecular weight
	tau = 0.56 # optical depth
	rsun = 69580000000 # cm
	rjup = 6991100000 # cm

	# Read in the file from the directory specified by the user
	# define the columns of data from the file read in 
	model_filename = os.path.join(model_dir,model_file)
	print(model_filename)

	model_wav, model_rp = np.loadtxt(model_filename, unpack=True)

	model_rprs = np.sqrt(model_rp) * (radius_planet/radius_star)

	rstar = radius_star * rsun
	rpl = radius_planet * rjup
	rp2 = rpl

	# Extract the starting temperature and gravity from the input filename
	temperature_model = float(model_file.split('_')[1]) # Pull the temperature of the model from the filename
	gravity_model = float(model_file.split('_')[2]) # Pull the gravity of the model from the filename
	gravity_model = gravity_model*1e2 # gravity of the selected model rescaled to cm

	# Calculate the transmission spectrum based on the file parameters
	#	 to get the baseline radius
	h1  = (kb * temperature_model) / (mu_model * gravity_model) # Calculate the scale height of the model atmosphere
	r1  = np.sqrt(model_rp) * rsun #cm - observed radius
	rp1 = rjup #cm - models use Rp = Rjup for bulk radius
	z1  = r1 - rp1 #cm
	epsig1 = tau * np.sqrt((kb * temperature_model * mu_model * gravity_model) / (2. * np.pi * rp1)) * np.exp(z1 / h1)

	# Rescale the model atmosphere based on the input parameters
	h2 = (kb * temperature_planet) / (mu_planet * gravity_planet)
	z2 = h2 * np.log(epsig1 / tau * np.sqrt((2. * np.pi * rp2)/(kb * temperature_planet * mu_planet * gravity_planet)))
	r2 = z2 + rp2

	# Resort the data so it goes from lowest wavelength to highest
	srt = np.argsort(model_wav)
	model_wavelength = model_wav[srt]
	radius_sort = r2[srt]

	transit_depth = (radius_sort / rstar)**2 # Output is in (Rp/Rs)^2

	# These are the OUTPUTS
	return model_wavelength, transit_depth














if __name__ == '__main__':
	model_dir = './Data'
	model_file = 'trans-iso-generic_0400_05_+0.0_0.35_0001_0.06_model.txt'
	# model_file = 'trans-iso-generic_1500_05_+0.0_0.56_0001_0.00_model.txt'
	
	# Set up your input values for the planet system
	radius_star = 1.100 # Rstar
	radius_planet = 1.200 # RJ
	gravity_planet = 610 # cgs
	temperature_planet = 450 # K
	# mean_molecular_weight = 2.3 
	# If mean_molecular_weight is not set it will use 2.3 as a defult (H+He atmosphere)
	
	# result = gg_rescale(model_dir, model_file, rs, rp, gp, tp)
	# This will give you a tuple with wavelength and transit depth.
	# Or do the following to put them into your own arrays
	wav_default_mmw, rprs_default_mmw = gg_rescale(model_dir, model_file, radius_star, radius_planet, gravity_planet, temperature_planet)



	# Run 2 with your own input mean molecular weight
	radius_star = 1.100 # Rstar
	radius_planet = 1.200 # RJ
	gravity_planet = 610 # cgs
	temperature_planet = 450 # K
	mean_molecular_weight = 5 # Rescale the atmosphere to your own mean molecular weight.

	wav_set_mmw, rprs_set_mmw = gg_rescale(model_dir, model_file, radius_star, radius_planet, gravity_planet, temperature_planet, mean_molecular_weight)

	
	# Plot the default mean molecular weight atmosphere against the scaled one
	plt.plot(wav_default_mmw, rprs_default_mmw)
	plt.plot(wav_set_mmw, rprs_set_mmw, color='red')
	plt.xlim(0.3,15)
	plt.xscale('log')
	plt.tight_layout()
	plt.show()



