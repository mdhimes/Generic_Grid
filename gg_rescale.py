# ggg_rescale.py
import sys, os
import itertools
import numpy as np 
import scipy.interpolate as si
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
    # from the documentation in the full grid:
    #    (use mean_mol_weight = 2.3, 2.5, 3.2, 4.0 and 5.5 for models closest to 1, 10, 50, 100 and 200 times solar metallicity models, respectively.) 
	mus = {'-1.0' : 2.3, '+0.0' : 2.3, '+1.0' : 2.5, '+1.7' : 3.2, '+2.0' : 4.0, '+2.3' : 5.5}

	# SET UP THE CONSTANTS
	kb = 1.380658E-16 # gm*cm^2/s^2 * Kelvin
	mu_model = 1.6726E-24 * mus[model_file.split('_')[3]] #g  cgs  Hydrogen + Helium Atmosphere
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


def gg_rescale_plus(model_dir, radius_star, radius_planet, gravity_planet, temperature_planet, mean_molecular_weight, C_to_O_planet, haze_planet, cloud_planet, powabunga=6):
	"""
	gg_rescale_plus
	This function takes a directory with the ATMO grid files, and some desired planetary parameters.
	The grid models are used to compute an inverse-distance-weighted spectrum for the given parameters.

	INPUTS:
	model_dir = string for the directory containing the file you want to scale
		e.g. '/Users/name/Documents/Generic_Goyal_Grid/Rainout_condensation/'
	radius_star = radius of the star in solar radii
	radius_planet = radius of the planet in jupiter radii
	gravity_planet = gravity of the planet in units of m/s
	temperature_planet = temperature of the planet in K
	mean_molecular_weight = scale factor for the atmospheric mean molecular weight. This will be multiplied by the mass of a proton and used to calculate the planetary scale height. The default is set to 2.3 for H/He-dominated atmosphere.
	C_to_O_planet = C/O for the planet
	haze_planet = planetary haze parameter
	cloud_planet = planetary cloud parameter
	powabunga = inverse-distance weighting exponent. Weight is 1/dist**powabunga.
				the default, 6, sets this exponent to the number of grid parameters.

	OUTPUTS:
	model_wavelength = wavlength array for the model (microns)
	transit_depth = (Rp/R*)^2 model
	"""
	# Grid parameters
	temps_str = np.array(['0300', '0400', '0500', '0600', '0700', '0800', 
	                      '0900', '1000', '1100', '1200', '1300', '1400', 
	                      '1500', '1600', '1700', '1800', '1900', '2000', 
	                      '2100', '2200', '2300', '2400', '2500', '2600'])
	gravs_str = np.array(['05', '10', '20', '50'])
	logZs_str = np.array(['-1.0', '+0.0', '+1.0', '+1.7', '+2.0', '+2.3'])
	CtoOs_str = np.array(['0.35', '0.56', '0.70', '1.00'])
	hazes_str = np.array(['0001', '0010', '0150', '1100'])
	cloud_str = np.array(['0.00', '0.06', '0.20', '1.00'])
	temps = temps_str.astype(int)
	gravs = gravs_str.astype(int)
	logZs = logZs_str.astype(float)
	CtoOs = CtoOs_str.astype(float)
	hazes = hazes_str.astype(int)
	cloud = cloud_str.astype(float)
    # from the documentation in the full grid:
    #    (use mean_mol_weight = 2.3, 2.5, 3.2, 4.0 and 5.5 for models closest to 1, 10, 50, 100 and 200 times solar metallicity models, respectively.) 
	mus_dict = {'-1.0' : 2.3, '+0.0' : 2.3, '+1.0' : 2.5, '+1.7' : 3.2, '+2.0' : 4.0, '+2.3' : 5.5}
	mus = np.array([2.3, 2.3, 2.5, 3.2, 4.0, 5.5]) # for each of the logZ vals
	mu_interp = si.interp1d(mus, logZs) # to find logZ value from user input mean mol weight
	logZ_planet = mu_interp(mean_molecular_weight)
	# For distance weighting later
	temps_unit = scale(temps, temps.min(), temps.max(), [0, len(temps)-1])
	gravs_unit = scale(gravs, gravs.min(), gravs.max(), [0, len(gravs)-1])
	logZs_unit = scale(logZs, logZs.min(), logZs.max(), [0, len(logZs)-1])
	CtoOs_unit = scale(CtoOs, CtoOs.min(), CtoOs.max(), [0, len(CtoOs)-1])
	hazes_unit = scale(hazes, hazes.min(), hazes.max(), [0, len(hazes)-1])
	cloud_unit = scale(cloud, cloud.min(), cloud.max(), [0, len(cloud)-1])
	temperature_planet_unit = scale(temperature_planet, temps.min(), temps.max(), [0, len(temps)-1])
	gravity_planet_unit     = scale(gravity_planet,     gravs.min(), gravs.max(), [0, len(gravs)-1])
	logZ_planet_unit        = scale(logZ_planet,        logZs.min(), logZs.max(), [0, len(logZs)-1])
	C_to_O_planet_unit      = scale(C_to_O_planet,      CtoOs.min(), CtoOs.max(), [0, len(CtoOs)-1])
	haze_planet_unit        = scale(haze_planet,        hazes.min(), hazes.max(), [0, len(hazes)-1])
	cloud_planet_unit       = scale(cloud_planet,       cloud.min(), cloud.max(), [0, len(cloud)-1])

	# SET UP THE CONSTANTS
	kb = 1.380658E-16 # gm*cm^2/s^2 * Kelvin
	mu_planet = 1.6726E-24 * mean_molecular_weight #g  cgs rescaled atmospheric mean molecular weight
	tau  = 0.56 # optical depth
	rsun = 69580000000 # cm
	rjup = 6991100000 # cm

	# Determine which model files are needed
	fmodels = []
	# Closest indices
	it = np.argmin(np.abs(temperature_planet - temps))
	ig = np.argmin(np.abs(gravity_planet - gravs))
	iZ = np.argmin(np.abs(logZ_planet - logZs))
	ir = np.argmin(np.abs(C_to_O_planet - CtoOs))
	ih = np.argmin(np.abs(haze_planet - hazes))
	ic = np.argmin(np.abs(cloud_planet - cloud))

	# Next closest indices, if applicable
	if temps[it] - temperature_planet != 0:
		it2 = np.argmin(np.abs(np.delete(temps, it) - temperature_planet))
		if it <= it2:
			it2 += 1
		it2 = [it2]
	else:
		it2 = []
	it = [it]
	if gravs[ig] - gravity_planet != 0:
		ig2 = np.argmin(np.abs(np.delete(gravs, ig) - gravity_planet))
		if ig <= ig2:
			ig2 += 1
		ig2 = [ig2]
	else:
		ig2 = []
	ig = [ig]
	if logZs[iZ] - logZ_planet != 0:
		iZ2 = np.argmin(np.abs(np.delete(logZs, iZ) - logZ_planet))
		if iZ <= iZ2:
			iZ2 += 1
		iZ2 = [iZ2]
	else:
		iZ2 = []
	iZ = [iZ]
	if CtoOs[ir] - C_to_O_planet != 0:
		ir2 = np.argmin(np.abs(np.delete(CtoOs, ir) - C_to_O_planet))
		if ir <= ir2:
			ir2 += 1
		ir2 = [ir2]
	else:
		ir2 = []
	ir = [ir]
	if hazes[ih] - haze_planet != 0:
		ih2 = np.argmin(np.abs(np.delete(hazes, ih) - haze_planet))
		if ih <= ih2:
			ih2 += 1
	else:
		ih2 = []
	ih = [ih]
	if cloud[ic] - cloud_planet != 0:
		ic2 = np.argmin(np.abs(np.delete(cloud, ic) - cloud_planet))
		if ic <= ic2:
			ic2 += 1
		ic2 = [ic2]
	else:
		ic2 = []
	ic = [ic]

	gravity_planet *= 100 # m/s --> cm/s

	specs = []
	dists = []
	for T, g, Z, r, h, c in itertools.product(it+it2, 
                                              ig+ig2, 
                                              iZ+iZ2, 
                                              ir+ir2, 
                                              ih+ih2, 
                                              ic+ic2):
		fmodel = 'trans-iso-generic_'+temps_str[T]+'_'+gravs_str[g]+'_'+logZs_str[Z]+'_'+CtoOs_str[r]+\
                 '_'+hazes_str[h]+'_'+cloud_str[c]+'_model.txt'
		model_file = fmodel
		model_filename = os.path.join(model_dir, fmodel)
		fmodels.append(model_filename)

		# Read in the file from the directory specified by the user
		# define the columns of data from the file read in 
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
		mu_model      = 1.6726E-24 * mus_dict[model_file.split('_')[3]] #g  cgs
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
		specs.append(transit_depth)
		dist = ((temperature_planet_unit - temps_unit[T])**2 + 
                (gravity_planet_unit     - gravs_unit[g])**2 + 
                (logZ_planet_unit        - logZs_unit[Z])**2 + 
                (C_to_O_planet_unit      - CtoOs_unit[r])**2 + 
                (haze_planet_unit        - hazes_unit[h])**2 + 
                (cloud_planet_unit       - cloud_unit[c])**2)**0.5
		dists.append(dist)

	if len(specs) != len(dists):
		raise Exception("This should never happen.")

	# Now, weight by inverse distance
	if len(specs) == 1:
		transit_depth = specs[0]
	else:
		specs = np.asarray(specs)
		dists = np.asarray(dists)
		transit_depth = 0
		# normalize the weights
		weights = (1./dists**powabunga) / np.sum(1./dists**powabunga)
		if np.sum(weights) - 1 >= 1e-7:
			raise ValueError("Weights do not sum to 1 within the allowed tolerance.")
		for i in range(len(specs)):
			transit_depth += weights[i] * specs[i]

	# These are the OUTPUTS
	return model_wavelength, transit_depth


def scale(val, vmin, vmax, scalelims):
    """
    Scales a value according to min/max values and scaling limits.

    Inputs
    ------
    val      : array. Values to be scaled.
    vmin     : array. Minima of `val`.
    vmax     : array. Maxima of `val`.
    scalelims: list, floats. [min, max] of range of scaled data.

    Outputs
    -------
    Array of scaled data.
    """
    return (scalelims[1] - scalelims[0]) * (val - vmin) / \
           (vmax - vmin) + scalelims[0]












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



