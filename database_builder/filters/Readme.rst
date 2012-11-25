The .dat files in this directory are the filter transmissions that are added
to the pcigale database. Each filter file has this structure:
- The first line is the filter name beginning with a #;
- The second line is the filter type (energy or photon) beginning with a #;
- The third line is the filter description;
- Everything else is the transmission table with the wavelength in Å and the
  transmission value.

To add a new filter to the database, add a new file with this structure,
taking care to use a name that is not yet used, before launching the build script.
