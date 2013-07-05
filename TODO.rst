* Check the age of the Universe at observation redshift.
* Deal with observation without values in some filters.
* Code the theoretical SED serialisation and the database to store them.
* When a default parameter is an array, use the "eval ..." the the generated
  ini file.
* Change metallicity numbering (ex. 0.02 for solar).
* Code the creation of mock catalogues.
* Create some plot to check the fit quality:
  - plot giving for each band the ratio between the observed flux and the best
    sed flux.
  - colour plots (colours chosen by the user) comparing the colour ranges
    covered by the simulation and the colours of the observed galaxies.)
* Write documentation.
* Code IGM effect module.
* Make the IGM effect a module separated from the main SED creation modules
  and applied according to the redshift of the observation.
* Make the IGM module work with a list (or range) of redshift to allow the
  use of pcigale to compute photometric redshifts.
