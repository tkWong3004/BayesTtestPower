# On a Generalizable Approach for Sample Size Determination in Bayesian t Tests 
This repository contains codes for running the Shiny app:

  `app.r` The R script of the app.

  `onesample.r` , `twosample.r`, `boost_noncentralt.cpp`  and `pt.cpp` contains the needed functions for conducting BFDA behind the app.

The codes for generating the tables and figures:

  `probabilities.RData` The dataset contain the calculated power for generating the plots. 

  `Figures_and_Tables.R` The scripts for producing the tables and figures. 

### Dependencies

This repository relies on the following R packages:

`rootSolve` version ‘1.8.2.4’: Used for solving nonlinear equations and other root-finding tasks.

`shiny` version ‘1.8.1.1’: For creating interactive web applications in R.

`Rcpp` version ‘1.0.12’: To interface with C++ code for improved performance.

`BH` version ‘1.84.0.0’: A header-only package that provides the Boost C++ libraries.

# How to run the shiny app locally?
 
1.  Download and put all the files `app.r` `onesample.r` , `twosample.r`, `boost_noncentralt.cpp`  and `pt.cpp` in the same working directory.
   
2a.  Open `app.r` and run the app.
  
2b.  Alternatively, `ScriptForBFDA.r` is availiable for conducting BFDA. 
  

  

  
