# EE367_Final_Project
Author: Nicole Bienert

Purpose: Use ADMM to invert for conductivity (and reflection coeffecient) 
	in a 2D cross section of a glacier

Instructions: 
Inverting for only conductivity
	Run the code titled "ADMMTV_ver1_6_5_run_me"
		You can change which temperature distribution you want to 
			simulate by commenting or uncommenting the corresponding 
			distribution (lines 38-58)
		You can change the number of transects to be used by commenting 
			or uncommenting the corresponding lines (lines 59-99)

Applying the Winebrenner approach to invert for average conductivity
	Run the code titled WINEBRENNER2004

Inverting for reflection coeffecient and conductivity: 
	I made several attempts at this and could not get it working, but the 
		codes run and give simulated outputs.
	ADMMTV_ver2_1: I would suggest running this one. This one implements my 
			modified version of the ADMM framework. x includes 
			both the conductivities and reflection coeffecients. There 
			are additional terms in the minimization and constraint, 
			promoting sparseness in the reflection coeffecient and applying
			the TV prior to the conductivity
	ADMMTV_ver2_3: Sparseness didn't seem to be working, so this one uses
			Dx as the prior on reflection coeffecient
	ADMMTV_ver2_4: This one attempts to include internal layers for additional 
			information
	
