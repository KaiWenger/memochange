## Test environments
* local OS X install, R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Kai Wenger <kai.wenger@gmx.de>'

  The note occurs since this is my first submission.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Resubmission
This is a resubmission. In this version we have:

* Added more details to the description text 

* Added references describing the theoretical background of methods in the package in the Decription field of the DESCRIPTION file

* Added \value tags to the Rd files in which we explain the returned objects

* Deleted par() in all of our examples