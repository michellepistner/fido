# fido 1.0.0

More changes to make the final push for CRAN! Changes include:

* Updating the description file to match CRAN's standards
* Removing the dependency to driver
* Linking the one remote package (MicrobeDS) to a separate CRAN-like Github repo using the "drat" package. This is to match CRAN's policy against no use of "Remotes".
* Precomputing the longer running vignettes, cutting R CMD check times in half. The original vignettes can be found in the vignettes folder with a ".orig" extension for reproducibility.

# fido 0.1.14

* More changes to prepare for version 0.2 (and CRAN soon), including:
* Lots of cleaning of small warnings and notes
* Fixing erorrs with the base plot function
* Cleaning up some manual entries
* Fixing the errors posted on Github (to_ilr among others)

# fido 0.1.13

* tons of tiny changes to prepare for version 0.2 (and ultimately CRAN) featured changes include:
* plot and coef methods for *fit objects now abide by CRAN generic naming scheme - may cause
  some problems if prior code had positional arguments to these functions
* now being tested using travis integration (vignettes not tested there due to time constraints)
* merged Kim's fixes to maltipoo code (will almost certainly change maltipoo results; there
  had been a bug in prior versions)
* Internally stored mallard and mallard_family data now import as `mallard` and `mallard_family` 
  rather than both importing as `ps`. Solved a R CMD CHECK warning. 
* Lots of tiny updates to documentation 

# fido 0.1.12

* orthus added for jointly fitting multinomial and gaussian data (e.g., 16S microbiome and metabolomics data)

# fido 0.1.11

* Kim fixed windows installation problems (Thanks Kim!)
* Fixed bug in predict that was effecting ppc
* Fixed error when pars!=NULL in summary
* lambda_to_iqlr now had default behavior when focus.cov=NULL
* Fixed issue with Xi=NULL in pibblefit when trying to transform
* small changes to make pibble wrapper around c++ functions faster

# fido 0.1.10

* basset added for fitting non-linear regression using fido
* numerous changes to make installation easier
* Added a `NEWS.md` file to track changes to the package.
