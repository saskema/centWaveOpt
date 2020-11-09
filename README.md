# Supplementary Material
The code provided by this repository is part of the supplementary material for the paper "Automated Optimization of XCMS Parameters for Improved Peak Picking of LC-MS Data using the Coefficient of Variation and Parameter Sweeping for Untargeted Metabolomic" by Sascha K. Manier et al., Drug Testing and Analysis, 2018. See the following sections for detailed information about it.

__There is now a "devel" branch that contains a modified version of ___centWaveOpt.R___ in order to evaluate each results table automatically. However, this will not be merged into the master branch to preserve the script as it was used in the corresponding publication.__

# Dependencies
The code deployed by this repository require the following packages:

- XCMS (Bioconductor)
- CAMERA (Bioconductor)
- tidyverse (CRAN)
- stringr (CRAN)
- ggsignif (CRAN)

If not installed you can use the following lines to install them:

	install.packages(c("tidyverse", "stringr", "ggsignif"))
	source('https://bioconductor.org/biocLite.R')
	biocLite(c("xcms", "CAMERA"))

# centWaveOpt.R
This code was used to optimise the parameters of XCMS. In a first step the parameters were applied in a certain range and replaced the XCMS-Online Q-Exactive preset parameters (Tautenhahn et al., Anal Chem., 2012) one by one. In a second step the top three parameters (except prefilter 1 and 2) were combined in every possible way and the above described process will was repeated to exclude cross effects.

## How it works
XCMS readable files were placed in a folder called "QC" in the working directory. File endings were specified in the beginning of the script. The scrip creates a folder in the working directory called "centWaveOptResults" where it will save all the parameter sweeping results in CSV files.

## Recommendation
Opening each and every CSV file is somewhat cumbersome and time consuming. However during the course of the research this made sense for me since I needed to see the results in detail. If someone uses this for his own research he is suggested to condense the results into one text file containing only the necessary parameters. Since this is the supplementary data to the above mentioned publication I will keep this file at the state when the paper was created.

# XCMS.R
XCMS.R was used for preprocessing mzXML files. The Results are returned in a folder called "XCMSResults" as "features.csv".

# replicatesEval.R
replicatesEval.R was used to genereate boxplots (Figure 2 and Figure S1-5). Results of XCMS.R were collected and renamed according to the corresponding experiment (e.g. "Mix\_High\_Result\_A.csv").

# pseudoMetabEval.R
This script was used to evaluate pseudometabolomic experiments. After running "XCMS.R" two columns were inserted in the result file. The column "category" was inserted as the second column and contained information of the feature if it was unknown, an adduct, an isotope, or a parent compound. The column "identity" was inserted as the third column of the table and contained information about the the drug from which the feature derives according to the procedure described in the publication. After this preparation the script was run. The results were used to generate Figures 3 and 4, as well as the Figures S6 and S7.
