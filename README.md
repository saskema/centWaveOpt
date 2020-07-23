## centWaveOpt
This script is part of the supplementary material for the paper "Automated Optimization of XCMS Parameters for Improved Peak Picking of LC-MS Data using the Coefficient of Variation and Parameter Sweeping for Untargeted Metabolomic" by Sascha K. Manier et al., Drug Testing and Analysis, 2018.
Use it to initially optimise for your peak XCMS parameters. In a first step the parameters will be used in a certain range and replace the XCMS-Online Q-Exactive preset parameters (Tautenhahn et al., Anal Chem., 2012) one by one. In a second step the top three parameters (except prefilter 1 and 2) will be combined in every possible way and the above described process will be repeated to exclude cross effects.

## How to use it
Place your XCMS readable files in a folder called "QC" and place it in your working directory. Specify your file ending in the beginning of the script and run it. The script will create a folder in your working directory called "centWaveOptResults" where it will save all the parameter sweeping results in CSV files.
You will also find a file called "parameter.csv" that contains the proposed parameter values sorted by the order XCMS lists the parameters. If a value is "Inf" the corresponding median cv did not change during parameter sweeping.

## Dependencies
This script requires tidyverse and XCMS. If not installed you can use the following lines to install them:

	install.packages('tidyverse')
	source('https://bioconductor.org/biocLite.R')
	biocLite('xcms')


