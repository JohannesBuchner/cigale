#########################
How to use pcigale script
#########################

pcigale comes with a script (a programme callable from the command line) that
will help you to run the complete analysis of a set of galaxy with known
redshift and fluxes in various filters.

Generation of a minimal configuration files
===========================================

First, you need a minimal configuration file. Create a directory where you
will proceed with your analysis. In a shell, go inside this directory and
launch ``pcigale init``. This will create a *pcigale.ini* file in the
directory with that content:

.. code-block:: cfg
  :emphasize-lines: 4,9,14,17

  # File containing the observation data to be fitted. Each flux column
  # must have the name of the corresponding filter, the error columns are
  # suffixed with '_err'. The values must be in mJy.
  data_file =

  # Order of the modules use for SED creation. Available modules : bc03,
  # dh2002, dl2007, dustatt_calzleit, dustatt_powerlaw, igmattenuation,
  # lines, loadfile, m2005, sfh2exp, sfhfromfile.
  sed_modules = ,

  # Module used to redshift the SED before the integration in the filters.
  # This is a SED creation module that accepts a 'redshift' parameter (see
  # the documentation).
  redshift_module =

  # Method used for statistical analysis. Available methods: psum.
  analysis_method =

You need to edit this file to complete:

* *data_file* is the name of the file containing you observations. This is a
  FITS or VO-Table file containing the observations to be analysed. It must
  contain a column named *id* with an identifier for each source, a column
  named *redshift* and columns containing the fluxes in *mJy* in various
  filters. The flux columns must be named exactly as the corresponding filters
  in the pcigale database. When a column is the error associated with a
  filter, its name must prefix the filter name with *_err*.

* *sed_modules* is the list of pcigale SED creation modules, **in the right
  order**, that will be used to create the theoretical SED (see the various
  module documentations). The comments above the entry list the available
  modules.

* *redshift_module* is the name of the module that will be used to redshift
  and add IGM attenuation to the SED before computing the fluxes that will be
  compared to observed ones.

* *analysis_method* is the name of the statistical analysis module that will
  be used. Here again, the comments list the available modules.

Here is a sample completed minimal configuration file:

.. code-block:: cfg

  # File containing the observation data to be fitted. Each flux column
  # must have the name of the corresponding filter, the error columns are
  # suffixed with '_err'. The values must be in mJy.
  data_file = data.fits

  # Order of the modules use for SED creation. Available modules : bc03,
  # dh2002, dl2007, dustatt_calzleit, dustatt_powerlaw, igmattenuation,
  # lines, loadfile, m2005, sfh2exp, sfhfromfile.
  sed_modules = sfh2exp, m2005, dustatt_calzleit, dh2002

  # Module used to redshift the SED before the integration in the filters.
  # This is a SED creation module that accepts a 'redshift' parameter (see
  # the documentation).
  redshift_module = igmattenuation

  # Method used for statistical analysis. Available methods: psum.
  analysis_method = psum

Generation of the main configuration file
=========================================

You can now run the ``pcigale genconf`` command. pcigale will read the data
file to identify which columns contain filter fluxes and errors. It will then
pick each module and add their configuration section to the ini file, that
will become something like this:

.. code-block:: cfg

   # File containing the observation data to be fitted. Each flux column
   # must have the name of the corresponding filter, the error columns are
   # suffixed with '_err'. The values must be in mJy.
   data_file = data.fits

   # Order of the modules use for SED creation. Available modules : bc03,
   # dh2002, dl2007, dustatt_calzleit, dustatt_powerlaw, igmattenuation,
   # lines, loadfile, m2005, sfh2exp, sfhfromfile.
   sed_modules = sfh2exp, m2005, dustatt_calzleit, dh2002

   # Module used to redshift the SED before the integration in the filters.
   # This is a SED creation module that accepts a 'redshift' parameter (see
   # the documentation).
   redshift_module = igmattenuation

   # Method used for statistical analysis. Available methods: psum.
   analysis_method = psum

   # List of the columns in the observation data file to use for the
   # fitting.
   column_list = WFI_U38, WFI_U38_err, WFI_U, WFI_U_err, WFI_B, WFI_B_err,
                 WFI_V, WFI_V_err, WFI_R, WFI_R_err, WFI_I, WFI_I_err, WFI_z,
                 WFI_z_err, IRAC1, IRAC1_err, IRAC2, IRAC2_err, IRAC3,
                 IRAC3_err, IRAC4, IRAC4_err, MIPS1, MIPS1_err, PACS_green,
                 PACS_green_err


   # Configuration of the SED creation modules.
   [sed_creation_modules]

     [[sfh2exp]]
       # Age of the oldest stars in the galaxy in Myr. The precision is 1 Myr.
       age =
       # e-folding time of the late starburst population model in Myr.
       tau_burst =
       # Mass fraction of the late burst population.
       f_burst =
       # e-folding time of the main stellar population model in Myr.
       tau_main =
       # Age of the late burst in Myr. Precision is 1 Myr.
       burst_age =

     [[m2005]]
       # Metallicity Z.
       metallicity =
       # Age [Myr] of the separation between the young and the old star
       # populations. The default value in 10^7 years (10 Myr). Set to 0 not to
       # differentiate ages (only an old population).
       separation_age = 10
       # Initial mass function: salp (Salpeter) or krou (Krupa)
       imf =

     [[dustatt_calzleit]]
       # Width (FWHM) of the UV bump in nm.
       uv_bump_width =
       # Reduction factor for the E(B-V)* of the old population compared to the
       # young one (<1).
       E_BVs_old_factor =
       # Central wavelength of the UV bump in nm.
       uv_bump_wavelength = 217.5
       # Slope delta of the power law modifying the attenuation curve.
       powerlaw_slope =
       # E(B-V)*, the colour excess of the stellar continuum light for the
       # young population.
       E_BVs_young =
       # Name of the contribution containing the spectrum of the young
       # population.
       young_contribution_name = m2005_young
       # List of the filters for which the attenuation will be computed.
       filters = V_B90, FUV
       # Name of the contribution containing the spectrum of the old
       # population. If it is set to 'None', only one population is considered.
       old_contribution_name = m2005_old
       # Amplitude of the UV bump in nm.
       uv_bump_amplitude =

     [[dh2002]]
       # List of attenuation value names (in the SED's info dictionary). A new
       # re-emission contribution will be added for each one.
       attenuation_value_names =
       # Alpha slope.
       alpha =


   # Set the 'redshift' parameter to None (or delete the line). If there
   # are other parameters, you must give only one value for each.
   [redshift_configuration]
     # If set to true, the cosmological dimming is applied to the fluxes.
     dimming = True
     # Parameter which scales the tau value at each wavelength.
     rtau = 1.0
     # Redshift to apply to the galaxy.
     redshift = 0.0


   # Configuration of the statistical analysis method.
   [analysis_configuration]
     # If true, save the best SED for each observation to a file.
     save_best_sed = False
     # List of the variables (in the SEDs info dictionaries) for which the
     # statistical analysis will be done.
     analysed_variables = (sfr, average_sfr)
     # If true, for each observation and each analysed variable plot the
     # probability density function.
     plot_pdf = False
     # Maximum number of bins used to compute the probability density
     # function. This is only used when saving or printing the PDF. If there
     # are less values, the probability is given for each one.
     pdf_max_bin_number = 50
     # If true, for each observation save a plot of the best SED and the
     # observed fluxes.
     plot_best_sed = False
     # If true, for each observation and each analysed variable plot the
     # value vs reduced chi-square distribution.
     plot_chi2_distribution = False
     # If true, for each observation and each analysed variable save the
     # probability density function.
     save_pdf = False
     # Type of storage used to cache the generate SED.
     storage_type = memory

If you don't want to use all the filters present in you data file, delete the
name of their columns (and error columns) from the *column_list* entry. Also,
if you want to use some filters but not the error, delete the error column
names.

To each module is associated a set of parameters (see
:doc:`/sed_modules/index`). You may give a list of possible values for each
one, separated by commas. pcigale will compute all the possible combination of
parameter values and create the corresponding SED that will be compared to
your observations.

.. warning::

 For some parameters, a *single value* is in fact a list, for instance the
 *analysed_variables*; in that case, use the *tuple* notation, between
 parenthesis.

.. todo::

 In fact the parenthesis notation needs to be implemented.

.. note::

 You may also want to evaluate Python code to generate the list of values. You
 can use a string (between quotation marks) beginning with *eval*. In this
 string, you can use *numpy* with *np*. For instance, if a parameter is set to
 "eval np.arange(100)", its possible will be from 0 to 99 every 1.

Once the configuration is complete, `pcigale check` will give you the number
of SED that will be computed. You can run the analysis with `pcigale run`.
