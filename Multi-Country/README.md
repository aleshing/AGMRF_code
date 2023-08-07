# Multi-Country Application
This directory contains code used to reproduce the application estimating U5MR
at the subnational level across multiple countries.

## Downloading Data
Reproducing the application requires access to data from the 2016, 2016, 2014, 
2015, 2015, and 2016 Standard DHS surveys from Burundi, Ethiopia, Kenya, Rwanda,
Tanzania, and Uganda respectively. Instructions for gaining 
access to these surveys are available at https://dhsprogram.com/.

Once one has access to the data from these surveys, they should download the 
Births Recode and the Geographic data sets. 
These should be the following files:
 - Burundi 2016: BUBR70DT.zip, BUGE71FL.zip
 - Ethiopia 2016: ETBR71DT.zip, ETGE71FL.zip 
 - Kenya 2014: KEBR72DT.zip, KEGE71FL.zip
 - Rwanda 2015: RWBR70DT.zip, RWGE72FL.zip
 - Tanzania 2015: TZBR7ADT.zip, TZGE7AFL.zip
 - Uganda 2016: UGBR7BDT.zip, UGGE7AFL.zip

Unzip these files and place the folders BUBR70DT, BUGE71FL, ETBR71DT, ETGE71FL, 
KEBR72DT, KEGE71FL, RWBR70DT, RWGE72FL, TZBR7ADT, TZGE7AFL, UGBR7BDT, UGGE7AFL 
into the Data directory.

## Running the code
To reproduce the results from the application, set the working directory to the
Code directory, and run the following scripts in the following order:
 - direct_estimates.R
 - smoothed_direct_separate.R
 - smoothed_direct_estimates.R
 - process_results.R
 - plot_results.R

To reproduce results for the more general adaptive bym2, run the following 
scripts in the following order:
 - direct_estimates.R
 - smoothed_direct_estimates_general.R
 - process_results_general.R
 
If direct_estimates.R has already been run, it does not need to be rerun.