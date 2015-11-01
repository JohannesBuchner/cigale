# Change Log

## [Unreleased]
### Added
- The pcigale-mock utility has been added to generate plots comparing the exact and pcigale-estimated parameters. This requires pcigale to be run beforehand with the pdf_analysis module and the mock_flag option set to True. (Denis Burgarella and Médéric Boquien)
- The pcigale-filter utility has been added to easily list, plot, add, and remove filters without having the rebuild the database entirely. (Médéric Boquien)
- It is now possible to analyse the flux in a band as a regular parameter. It can be useful for flux predictions. (Yannick Roehlly)
- The redshift can be a now used as a free parameter, enabling pcigale to estimate the photometric redshift. (Médéric Boquien)

### Changed
- The galaxy_mass parameter was very ambiguous. In reality it corresponds to the integral of the SFH. Consequently it has been renamed sfh.integrated. (Médéric Boquien)
- In the Calzetti attenuation module, add a warning saying that is the power law slope is different than 0, E(B-V) will no longer be the real one. (Yannick Roehlly)
- Add "B_B90" to the list of computed attenuation so that users can calculate the effective E(B-V). (Yannick Roehlly)
- Computing the parameters and their uncertainties through the histogram of the PDF is slow and can introduce biases in some cases. Rather, now the estimated values of the parameters and the corresponding uncertainties are simply computed from the weighted mean and standard deviation of the models that are at least 0.1% as likely as the best model to reproduce the observations. The differences in the estimates are very small except when very few models are used. (Médéric Boquien)
- Magic values to indicate invalid values (e.g. values lower than -99) are difficult to handle safely. They have been replaced with NaN wherever appropriate. The logic of the input flux file stays the same for the time being but the magic values are converted internally after reading it. Users are advised to replace magic values with NaN. The output files now use NaN instead of magic number to indicate invalid values. (Médéric Boquien)
- Rename the the AGN faction added by dale2014 module from agn.fracAGN to agn.fracAGN_dale2014 to avoid conflict with fritz2006 module. (Yannick Roehlly)

### Fixed
- The SFH is modelled using discrete star formation episodes every 1 Myr. This means that as the SFH is not really continuous (the input single stellar population do not allow us to compute that properly), we should not integrate SFH(t), but simply sum SFH(t) as t is discrete. In most cases the difference is very small. The only case where it makes a difference is for a very rapidly varying SFH, for instance taking τ=1 Myr. (Médéric Boquien)
- Ensure that the flux can be computed even if the redshifting module has not been applied. By default in that case we assume a distance of 10 parsecs. While in practice it should never happen as the redshifting module is mandatory, this can be more important when using pcigale as a library. (Médéric Boquien and  Yannick Roehlly)
- When called without arguments, pcigale would crash. Now it displays a brief message to remind how it should be invoked. (Médéric Boquien)
- Raise an exception instead of crash when an unknown IMF is requested. (Médéric Boquien)
- When the error column for a flux was present in the observation table but not in the used column list (when the user prefers to use a default error instead of the real one), the error on the flux was set to 0. (Yannick Roehlly)
- For some reason a point in the GALEX FUV filter has a negative transmission. That should not happen. After comparison with the filter on the GALEX website it has been set to 0. (Médéric Boquien)
- Shorten the left and right 0 parts of the pseudo D4000 filter so that it can be applied on smaller spectra. (Yannick Roehlly)
- To compute the reduced χ², we need to divide by the number of bands-1 (and not the number of bands). We do that because we consider that the models depend on one meta-parameter. (Médéric Boquien)
- The test to determine when to take into account upper limits did not work according the specifications. Now upper limits are always taken into account when they should. (Médéric Boquien)

### Optimised
- Major speedup to build the database by inserting multiple models in the database at once rather one model at a time. On an SSD, the total run time of "python setup.py build" goes from 5m20s to 2m42s. The speedup should be even more spectacular on a rotating hard drive. (Médéric Boquien)
- Memory usage reduction using in-place operations (e.g. a/=2 rather than a=a/2, saving the creation of a temporary array the size of a) where possible. (Médéric Boquien)
- The creation and handling of mock catalogues has been streamlined. (Médéric Boquien)
- Slight speedup using np.full() where possible to create an array with all elements set to the same value. (Médéric Boquien)
- Computing the scaling factors and the χ² in one step over the entire grid of models is very memory-intensive, leading to out-of-memory issues when using multiple cores. Rather, let's compute them band by band, as this avoids the creation of large temporary arrays, while keeping the computation fast. (Médéric Boquien).
- Each core copied the subset of models corresponding to the redshift of the object to be analysed. This is a problem as it can strongly increase memory usage with the number of cores, especially when there are many models and just one redshift. Rather than making a copy, we use a view, which not only saves a considerable amount of memory but is also faster as there is no need to allocate new, large arrays. This is made possible as models are regularly ordered with redshift. (Médéric Boquien)
- Various minor optimisations. (Médéric Boquien)
