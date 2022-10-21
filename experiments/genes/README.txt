Preprocessing scripts for data from Kemmeren et al. 2014 - "Large-Scale Genetic Perturbations Reveal Regulatory Networks and an Abundance of Gene-Specific Repressors". 

Scripts:
	makeData.R
	query_sgd.py (script to create or update sgd_naming.json)

Output:
	Kemmeren.HDF5
		obs 		observational data
		int		interventional data
		plateobs 	plate extraction for observational data
		plateint 	plate extraction for interventional data
		intpos		intervention position in genenames
		p 		number of total variables
		nObs		number of observations
		nInt		number of interventions
		genenames	naming for variables
		mutantnames	naming for interventions (subset of genenames)
	Kemmeren.RData
		/obser/data	observational data
		/obser/extplate	 plate extraction for observational data 
		/inter/data	interventional data
		/inter/extplate	 plate extraction for interventional data
		/ids/genes	naming for variables
		/ids/mutants	naming for interventions (subset of genenames)
	sgd_naming.json (contains unique gene naming queried from SGD)