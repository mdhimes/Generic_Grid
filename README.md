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
		where the file naming convention is trans-iso-generic_temperature_gravity_metallicity_CtoO_haze_cloud_model.txt

	rs = radius of the star in solar radii
	rp = radius of the planet in jupiter radii
	gp = gravity of the planet in cgs units
	tp = temperature of the planet in K

	OUTPUTS:
	wav = wavlength array for the model (microns)
	rprs = (Rp/R*)^2 model
	
