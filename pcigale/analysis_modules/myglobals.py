# This is a dummy module to store the data that need to be shared between all
# processes. In pratice, once set these data should never change. Note that we
# cannot initialise thoses variables in subprocesses.

# Variables known at init
warehouse = []
creation_modules = []
creation_modules_params = []
analysed_variables = []
save_best_sed = []
save_chi2 = []
save_pdf = []
filters = []
n_models = 0
n_obs = 0

# Variables that need to be initilised after having computed the models
model_fluxes = []
model_variables = []
model_redshifts = []
model_info = []

mass_proportional_info = []
has_sfh = []
info_keys = []

redshifts = []
w_redshifts = []

# Misc
t_begin = 0.

# Variables modified by forked processes
n_computed = []