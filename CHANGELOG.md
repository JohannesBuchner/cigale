# Change Log

## 0.9.0 (2016-04-04)
### Added
- When using the `savefluxes` module, all the output parameters were saved. This is not efficient when the user is only interested in some of the output parameters but not all. We introduce the "variables" configuration parameter for `savefluxes` to list the output parameters the user wants to save. If the list is left empty, all parameters are saved, preserving the current behaviour. This should increase the speed substantially when saving memory. (Médéric Boquien)
- Similarly to the `savefluxes` module, in the `pdf_analysis` module if the list of physical properties is left empty, all physical parameters are now analysed. (Médéric Boquien)
- It is now possible to pass the parameters of the models to be computed from a file rather than having to indicate them in pcigale.ini. This means that the models do not necessarily need to be computed on a systematic grid of parameters. The name of the file is passed as an argument to the parameters\_file keyword in pcigale.ini. If this is done, the creation\_modules argument is ignored. Finally, the file must be formatted as following: each row is a different model and each column a different parameter. They must follow the naming scheme: module\_name.parameter\_name, that is "bc03.imf" for instance. (Médéric Boquien)
- Addition of the `schreiber2016` SED creation module implementing the Schreiber et al. (2016) dust models. (Laure Ciesla)
- The physical parameters provided in pcigale.ini were not checked at startup against what the modules could accept. This could lead to a runtime crash if an unexpected value was passed to the module. Now the parameters are checked at startup. If an issue is found, it is indicated and the user is asked to fix it before launching cigale again. The validation file is build at the same time as pcigale.ini. (Médéric Boquien)
- There is a new attenuation module (`dustatt_2powerlaws`) based on a double power law: a birth cloud power law applied only to the young star population and an ISM power law applied to both the young and the old star population. (Yannick Roehlly & Véronique Buat)

### Changed
- The estimates of the physical parameters from the analysis of the PDF and from the best fit were recorded in separate files. This can be bothersome when trying to compare quantities from different files. Rather, we generate a single file containing all quantities. The ones estimated from the analysis of the PDF are prefixed with "bayes" and the ones from the best fit with "best". (Médéric Boquien)
- To homogenize input and output files, the "observation\_id" has been changed to "id" in the output files. (Médéric Boquien)
- The output files providing estimates of the physical properties are now generated both in the form of text and FITS files. (Médéric Boquien)
- When using the `dustatt_calzleit module`, choosing ẟ≠0 leads to an effective E(B-V) different from the one set by the user. Now the E(B-V) will always correspond to the one specified by the user. This means that at fixed E(B-V), A(V) depends on ẟ. (Médéric Boquien)
- The pcigale-mock tool has been merged into pcigale-plots; the mock plots can be obtained with the "mock" command.
- The `sfhdelayed` module is now initialised with _init_code() to be consistent with the way things are done in other modules. This should give a slight speedup under some sircumstances too. (Médéric Boquien)
- In `sfhfromfile`, the specification of the time grid was vague and therefore could lead to incorrect results if it was not properly formatted by the end user. The description has been clarified and we now check that the time starts from 0 and that the time step is always 1 Myr. If it is not the case we raise an exception. (Médéric Boquien)
- When the redshift is not indicated in pcigale.ini, the analysis module fills the list of redshifts from the redshifts indicated in the input flux file. This is inefficient as analysis modules should not have to modify the configuration. Now this is done when interpreting pcigale.ini before calling the relevant analysis module. As a side effect, "pigale check" now returns the total number of models that cigale will compute rather than the number of models per redshift bin. (Médéric Boquien)
- The optionally saved spectra in the `pdf_analysis` and `savefluxes` modules were saved in the VO-table format. The most important downside is that it is very slow to write to, which proves to be a major bottleneck in the computing speed. To overcome this issue, we rather save the spectra using the FITS formation. Instead of having one file containing the spectra (including the various components) and the SFH in a single file, now we have one file for the spectra and one file for the SFH.

### Fixed
- To estimate parameters in log, pcigale determines which variables end with the "\_log" string and removed it to compute the models. However in some circumstances, it was overzealous. This has been fixed. (Médéric Boquien)
- When estimating a parameter in log, these were not scaled appropriately and taken in log when saving the related χ² and PDF. (Médéric Boquien)
- In the presence of upper limits, correct the scaling factor of the models to the observations before computing the χ², not after. (Médéric Boquien)
- When called without arguments, pcigale-plots would crash and display the backtrace. Now it displays the a short help showing how to use it. (Médéric Boquien)
- For `sfh2exp`, when setting the scale of the SFH with sfr0, the normalisation was incorrect by a factor exp(-1/tau\_main). (Médéric Boquien)
- The mass-dependent physical properties are computed assuming the redshift of the model. However because we round the observed redshifts to two decimals, there can be a difference of 0.005 in redshift between the models and the actual observation if CIGALE computes the list of redshifts itself. At low redshift, this can cause a discrepancy in the mass-dependent physical properties: ~0.35 dex at z=0.010 vs 0.015 for instance. Therefore we now evaluate these physical quantities at the observed redshift at full precision. (Médéric Boquien, issue reported by Samir Salim)
- In the `sfhfromfile` module, an extraneous offset in the column index made that it took the previous column as the SFR rather than the selected column. (Médéric Boquien)
- In `sfhfromfile`, if the SFR is made of integers cigale crashed. Now we systematically convert it to float. (Médéric Boquien)
- The order of the parameters for the analysis modules would change each time a new pcigale.ini was generated. Now the order is fixed. (Médéric Boquien)
- In the output the sfh.age parameter would correspond to the input value minus 1. Now both values are consistent with one another. (Laure Ciesla & Médéric Boquien)
- In rare circumstances requiring a specific distribution of redshifts the integration of the spectrum in some filters was not done correctly, inducing relative errors of ~10¯⁵-10¯⁶. (Médéric Boquien)
- The absorption of the Lyman continuum from old stars tended to be overestimated leading to some “negative fluxes” for the Lyman continuum. (Médéric Boquien)

### Optimised
- Prior to version 0.7.0, we needed to maintain the list of redshifts for all the computed models. Past 0.7.0 we just infer the redshift from a list unique redshifts. This means that we can now discard the list of redshifts for all the models and only keep the list of unique redshifts. This saves ~8 MB of memory for every 10⁶ models. the models should be computed slightly faster but it is in the measurement noise. (Médéric Boquien)
- The `sfhfromfile` module is now fully initialised when it is instantiated rather than doing so when processing the SED. This should be especially sensitive when processing different SED. (Médéric Boquien)
- We do not store the time grid in the SED anymore given that we assume it starts at 0 Myr with steps of 1 Myr, we can easily reconstruct to save it if needed. It should save a little bit of memory and it should go a little bit faster. (Médéric Boquien)
- To compute the stellar spectrum of the young component, do not pass the full SFH with the old part set to 0. Rather, only pass the corresponding part of the SFH. This nearly doubles the computing speed of the stellar spectrum (Médéric Boquien)
- Computers are much better at multiplying than at dividing. Therefore to correct the spectral emission when redshifting we multiply by 1/(1+z) rather than dividing by 1+z. (Médéric Boquien)

---

## 0.8.1 (2015-12-07)
### Fixed
-  To estimate parameters in log, pcigale determines which variables end with the "\_log" string and removed it to compute the models. However in some circumstances, it was overzealous. This has been fixed. (Médéric Boquien)

## 0.8.0 (2015-12-01)
### Added
- The evaluation of the parameters is always done linearly. This can be a problem when estimating the SFR or the stellar mass for instance as it is usual to estimate their log rather. Because the log is non-linear, the likelihood-weighted mean of the log is not the log of the likelihood-weighted mean. Therefore the estimation of the log of these parameters has to be done during the analysis step. This is now possible. The variables to be analysed in log just need to be indicated with the suffix "_log", for instance "stellar.m\_star\_log". (Médéric Boquien, idea suggested by Samir Salim)

### Fixed
- Running the scripts in parallel trigger a deadlock on OS X with python 3.5. A workaround has been implemented. (Médéric Boquien)
- When no dust emission module is used, pcigale genconf complains that no dust attenuation module is used. Correctly specify dust emission and not attenuation. (Médéric Boquien and Laure Ciesla)
- Allowing more flexibility to read ASCII files broke the handling of FITS files. It is now fixed. (Yannick Roehlly)

### Changed
- The attenuation.ebvs\_main and attenuation.ebvs\_old parameters are no longer present as they were duplicates of attenuation.E\_BVs.stellar.old and attenuation.E\_BVs.stellar.young (that are still available).

---

## 0.7.0 (2015-11-19)
### Added
- The pcigale-mock utility has been added to generate plots comparing the exact and pcigale-estimated parameters. This requires pcigale to be run beforehand with the `pdf_analysis` module and the mock\_flag option set to True. (Denis Burgarella and Médéric Boquien)
- The pcigale-filter utility has been added to easily list, plot, add, and remove filters without having the rebuild the database entirely. (Médéric Boquien)
- It is now possible to analyse the flux in a band as a regular parameter. It can be useful for flux predictions. (Yannick Roehlly)
- The redshift can be a now used as a free parameter, enabling pcigale to estimate the photometric redshift. (Médéric Boquien)
- When running "pcigale genconf", the list of modules is automatically checked against the list of official modules. If modules are missing, information is printed on the screen indicating the level of severity (information, warning, or error) and the list of modules that can be used. (Médéric Boquien)

### Changed
- The galaxy\_mass parameter was very ambiguous. In reality it corresponds to the integral of the SFH. Consequently it has been renamed sfh.integrated. (Médéric Boquien)
- In the Calzetti attenuation module, add a warning saying that is the power law slope is different than 0, E(B-V) will no longer be the real one. (Yannick Roehlly)
- Add "B\_B90" to the list of computed attenuation so that users can calculate the effective E(B-V). (Yannick Roehlly)
- Computing the parameters and their uncertainties through the histogram of the PDF is slow and can introduce biases in some cases. Rather, now the estimated values of the parameters and the corresponding uncertainties are simply computed from the weighted mean and standard deviation of the models that are at least 0.1% as likely as the best model to reproduce the observations. The differences in the estimates are very small except when very few models are used. (Médéric Boquien)
- Magic values to indicate invalid values (e.g. values lower than -99) are difficult to handle safely. They have been replaced with NaN wherever appropriate. The logic of the input flux file stays the same for the time being but the magic values are converted internally after reading it. Users are advised to replace magic values with NaN. The output files now use NaN instead of magic number to indicate invalid values. (Médéric Boquien)
- Rename the the AGN faction added by dale2014 module from agn.fracAGN to agn.fracAGN\_dale2014 to avoid conflict with fritz2006 module. (Yannick Roehlly)
- Remove the radio component from the dale2014 model so that it can be used with the more flexible radio module, courtesy Daniel Dale. (Laure Ciesla and Médéric Boquien)

### Fixed
- The SFH is modelled using discrete star formation episodes every 1 Myr. This means that as the SFH is not really continuous (the input single stellar population do not allow us to compute that properly), we should not integrate SFH(t), but simply sum SFH(t) as t is discrete. In most cases the difference is very small. The only case where it makes a difference is for a very rapidly varying SFH, for instance taking τ=1 Myr. (Médéric Boquien)
- Ensure that the flux can be computed even if the redshifting module has not been applied. By default in that case we assume a distance of 10 parsecs. While in practice it should never happen as the redshifting module is mandatory, this can be more important when using pcigale as a library. (Médéric Boquien and Yannick Roehlly)
- When called without arguments, pcigale would crash. Now it displays a brief message to remind how it should be invoked. (Médéric Boquien)
- Raise an exception instead of crash when an unknown IMF is requested. (Médéric Boquien)
- When the error column for a flux was present in the observation table but not in the used column list (when the user prefers to use a default error instead of the real one), the error on the flux was set to 0. (Yannick Roehlly)
- For some reason a point in the GALEX FUV filter has a negative transmission. That should not happen. After comparison with the filter on the GALEX website it has been set to 0. (Médéric Boquien)
- Shorten the left and right 0 parts of the pseudo D4000 filter so that it can be applied on smaller spectra. (Yannick Roehlly)
- To compute the reduced χ², we need to divide by the number of bands-1 (and not the number of bands). We do that because we consider that the models depend on one meta-parameter. (Médéric Boquien)
- The test to determine when to take into account upper limits did not work according the specifications. Now upper limits are always taken into account when they should. (Médéric Boquien)
- The nebular emission could be counted in excess in the dust luminosity as both lines and the Lyman continuum could be attenuated. Now we do not extend the attenuation under 91 nm. Also, a new component as been added, taking specifically the Lyman continuum absorption by gas, allowing to conserve the information about the intrinsic stellar Lyman continuum if need be. (Yannick Roehlly and Médéric Boquien)
- The Flambda table in the VO-table export does not reflect the fact that it stores luminosity densities. Accordingly, it has been renamed Llambda. (Yannick Roehlly and Médéric Boquien)
- When the flux file contains a mix of spaces and tabulations as column separators, pcigale discards the header and takes the first data line as the header. Now pcigale properly handles such a combination. Bug reported by Paola Santini. (Médéric Boquien)

### Optimised
- Major speedup to build the database by inserting multiple models in the database at once rather one model at a time. On an SSD, the total run time of "python setup.py build" goes from 5m20s to 2m42s. The speedup should be even more spectacular on a rotating hard drive. (Médéric Boquien)
- Memory usage reduction using in-place operations (e.g. a/=2 rather than a=a/2, saving the creation of a temporary array the size of a) where possible. (Médéric Boquien)
- The creation and handling of mock catalogues has been streamlined. (Médéric Boquien)
- Slight speedup using np.full() where possible to create an array with all elements set to the same value. (Médéric Boquien)
- Computing the scaling factors and the χ² in one step over the entire grid of models is very memory-intensive, leading to out-of-memory issues when using multiple cores. Rather, let's compute them band by band, as this avoids the creation of large temporary arrays, while keeping the computation fast. (Médéric Boquien).
- Each core copied the subset of models corresponding to the redshift of the object to be analysed. This is a problem as it can strongly increase memory usage with the number of cores, especially when there are many models and just one redshift. Rather than making a copy, we use a view, which not only saves a considerable amount of memory but is also faster as there is no need to allocate new, large arrays. This is made possible as models are regularly ordered with redshift. (Médéric Boquien)
- Various minor optimisations. (Médéric Boquien)

---

## 0.6.0 (2015-09-07)
### Added
- New module to compute a star formation history as described in Buat et al. 2008. (Yannick Roehlly)
- New module to compute a periodic SFH. Each star formation episode can be exponential, "delayed", or rectangular. (Médéric Boquien and Denis Burgarella)
- New module performing a quench on the star formation history. (Yannick Roehlly)
- New module to compute the emission of a modified black body. (Denis Burgarella)
- New module to compute the physical properties measured on the emission spectrum (e.g. spectral indices, ultraviolet slope, etc.). (Denis Burgarella)
- New pseudo filters to compute line fluxes and spectral indices. (Yannick Roehlly)
- The dust masses are now computed for the draine2007 and draine2014 modules. (Médéric Boquien)

### Changed
- The nebular\_lines\_width parameter of the nebular module is now called lines\_width as the nebular prefix was redundant. (Médéric Boquien)
- Prefix the ouput variables of SFH-related modules with "sfh" to facilitate their identification in the output files. (Médéric Boquien)
- Prefix the output variables of the fritz2006 AGN module with "agn" to facilitate their identification in the output files. (Médéric Boquien)
- Prefix the redshift with "universe". (Médéric Boquien)
- With pcigale-plot, draw the spectra only to λ=50 cm as the models do not extend much further and there is very rarely any observation beyond λ=21 cm. (Médéric Boquien)
- As pcigale is getting much faster, display the number of computed models every 250 models rather than every 100 models. (Médéric Boquien)
- Give default values for the dl2014, sfh\_buat, and sfhdelayed modules to allow for quick test runs. (Médéric Boquien)
- Now pcigale-plots plots errors on upper limits. (Denis Burgarella)

### Fixed
- When plotting, round the redshift to two decimals to match the redshift of the model. (Médéric Boquien)
- Ensure that all the input parameters of the nebular module are also output so it is possible to analyse them. (Médéric Boquien)
- Properly take the Lyman continuum photons absorbed by dust into account to compute the dust emission. (Médéric Boquien)
- Improve the readability of the pcigale-plots generated spectra by drawing the observed fluxes on top of other lines. (Médéric Boquien)
- The nebular components are now plotted with pcigale-plots. (Médéric Boquien)
- When a filter that was not in the database was called, pcigale would crash ungracefully as the exception invoked does not exist in Python 3.x. Now use a Python 3.x exception. (Médéric Boquien)
- The dustatt\_powerlaw module could not identify which physical components to attenuate when the nebular module was called. (Médéric Boquien)
- The displayed counter for the number of objects already analysed could be slightly offset from the number of models actually computed. (Médéric Boquien)
- Change the method to compute the χ² in the presence of upper limits as the previous method did not always converge. (Denis Burgarella and Médéric Boquien)

### Optimised
- When plotting, do not recompute the luminosity distance, which is very slow, but rather get the one computed during the analysis and that is given in the output files. (Médéric Boquien)
- Adding new physical components with a wavelength sampling different than that of the pre-existing grid is slow as a common grid has to be computed and all components interpolated over it. The nebular lines are especially bad for that, owing to the fact that each line is sampled with 19 points. This is excessive as sampling over 9 points barely changes the fluxes while speeding up the computation of the models by ~20%. (Médéric Boquien)
- Rather than resampling the filter transmission on the wavelength grid every time the flux is computed in a given filter, put the resampled filter into a cache. (Médéric Boquien)
- Rather than recomputing every time the merged wavelength grid from two different wavelength grids (for instance when adding a new physical component or when integrating the spectrum in a filter), put the results in a cache. (Médéric Boquien)
- Before adding a new component to a SED, we first copy the original SED without that component from the cache. This copy can be very slow when done automatically by python. We rather do this copy manually, which is much faster. (Médéric Boquien)
- When adding a new physical component with a different wavelength sampling, rather than reinterpolating all the components over the new grid, compute the interpolation only for new wavelengths. (Médéric Boquien)
- Various minor optimisations. (Médéric Boquien)
- When computing the flux in filters, np.trapz() becomes a bottleneck of the code. A large part of the time is actually spent on safeguards and on operations for nD arrays. However here we only have 1D arrays and some variables can be cached, which allows some optimisations to compute fluxes faster. (Médéric Boquien)
- The output parameters of a model were stored in an ordered dictionary. While convenient to keep the order of insertion it is very slow as it is implemented in pure Python for versions up to 3.4. Rather we use a regular dictionary and we reorder the parameters alphabetically. (Médéric Boquien)
- To store the SED in memory and retrieve them later, we index them with the list of parameters used to compute them. We serialise those using JSON. However JSON is slow. As these data are purely internal, rather use marshal, which is much faster than JSON. (Médéric Boquien)

---

## 0.5.1 (2015-04-28)
### Changed
- Set the default dale2014 AGN fraction to 0 to avoid the accidental inclusion of AGN. (Denis Burgarella)
- Modify the name of the averaged SFR: two averaged SFRs over 10 (sfh.sfr10Myrs) and 100Myrs (sfh.sfr100Myrs). (Denis Burgarella)
- Improve the documentation of the `savefluxes` module. (Denis Burgarella)

### Fixed
- Correction of the x-axis limits. (Denis Burgarella)
- Fix the detection of the presence of the agn.fritz2006\_therm in pcigale-plots. (Denis Burgarella)
- Correct the wavelength in the SCUBA 450 μm filter. (Denis Burgarella)
- Install the ancillary data required to make plots. (Yannick Roehlly)

## 0.5.0 (2015-04-02)

---

## 0.4.0 (2014-10-09)

---

## 0.3.0 (2014-07-06)

---

## 0.2.0 (2014-06-10)

---

## 0.1.0 (2014-05-26)
