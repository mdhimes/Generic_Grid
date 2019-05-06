# Generic_Grid
Scaling function for the generic ATMO grid presented in Goyal, J., et al., 2018b, MNRAS, accepted for publication October 30th 2018.

The full grid of models can be found at https://drive.google.com/drive/folders/1ZFbkPdqg37_Om7ECSspSpEp5QrUMfA9J?usp=sharing

The data folder contains example models for this function. The name structure is kept consistent between models for better reading.

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
	
