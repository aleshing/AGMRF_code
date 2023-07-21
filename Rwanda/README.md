# Rwanda Application
This directory contains code used to reproduce the application estimating U5MR
at the national level in Rwanda. 

## Downloading Data
Reproducing the application requires access to data from the 1992, 2000, 2005, 
2008, 2010, and 2015 Standard DHS surveys from Rwanda. Instructions for gaining 
access to these surveys are available at https://dhsprogram.com/.

Once one has access to the data from these surveys, they should download the 
Births Recode data sets. These should be files RWBR21DT.zip, RWBR41DT.zip, 
RWBR53DT.zip, RWBR5ADT.zip, RWBR61DT.zip, RWBR70DT.zip for the 1992, 2000, 2005, 
2008, 2010, 2015 surveys respectively. Unzip these files and place the 
folders RWBR21DT, RWBR41DT, RWBR53DT, RWBR5ADT, RWBR61DT, RWBR70DT into the 
Data directory.

## Running the code
To reproduce the results from the application, set the working directory to the
Code directory, and run the following scripts in the following order:
 - direct_estimates.R
 - smoothed_direct_estimates.R
 - plot_results.R

To reproduce results when including a separate intercept for conflict years,
run the following scripts in the following order:
 - direct_estimates.R
 - smoothed_direct_estimates_conflict_int.R
 - plot_results_conflict_int.R
 
If direct_estimates.R has already been run, it does not need to be rerun.


