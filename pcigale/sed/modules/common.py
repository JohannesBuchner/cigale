# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


class SEDCreationModule(object):
    """Abstract class, the pCigale SED creation modules are based on.
    """

    # parameter_list is a dictionary containing all the parameters used by
    # the module. Each parameter name is associate to a tuple (variable type,
    # description [string], default value). Each module must define its
    # parameter list, unless it does not use any parameter. Using None means
    # that there is no description or default value. If None should be the
    # default value, use the 'None' string instead.
    parameter_list = {}

    # out_parameter_list is a dictionary containing all the SED parameters
    # that are added to the SED info dictionary and for which a statistical
    # analysis may be done. Each parameter name is associated with its
    # description. In the SED info dictionary, the parameter name in prefixed
    # with the name of the module plus an underscore (to allow several
    # modules to add a parameter with the same name, for instance a repeated
    # module.)
    out_parameter_list = {}

    # comments is the text that is used to comment the module section in
    # the configuration file. For instance, it can be used to give special
    # instructions for the configuration.
    comments = ""

    def __init__(self, name=None, **kwargs):
        """Instantiate a SED creation module

        A name can be given to the module. This can be useful when a same
        module is used several times with different parameters in the SED
        creation process.

        The module parameters values can be passed as keyworded paramatres.
        """
        self.name = name

        # parameters is a dictionary containing the actual values for each
        # module parameter.
        self.parameters = kwargs

    def _process(self, sed, parameters):
        """Do the actual processing of the module on a SED object

        This method is called with an object and a complete parameters
        dictionary. It is not meant to be called directly but through the
        process method.

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary of the module parameters

        """
        raise NotImplementedError()

    def process(self, sed, parameters=None):
        """Process a SED object with the module

        This method is responsible for checking the module parameters (whether
        they are given in the method call or are taken from parameters class
        attribute) before doing the actual processing (_process method). If a
        parameter is not given but exists in the parameter_list with a default
        value, this value is used.

        The SED object is updated during the process, one must take care of
        copying it before, if needed.

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary
            Dictionary containing the module parameter values, if it is not
            given, the module parameter values are used

        Raises
        ------
        KeyError : when not all the needed parameters are given.

        """

        # If the parameter dictionary is not passed, use the module one
        if not parameters:
            parameters = self.parameters

        # For parameters that are present on the parameter_list with a default
        # value and that are not in the parameters dictionary, we add them
        # with their default value.
        for key in self.parameter_list:
            if (not key in parameters) and (
                    self.parameter_list[key][2] is not None):
                parameters[key] = self.parameter_list[key][2]

        # If the keys of the parameters dictionary are different from the one
        # of the parameter_list dictionary, we raises a KeyError. That means
        # that a parameter is missing (and has no default value) or that an
        # unexpected one was given.
        if not set(parameters.keys()) == set(self.parameter_list.keys()):
            missing_parameters = (set(self.parameter_list.keys())
                                  - set(parameters.keys()))
            unexpected_parameters = (set(parameters.keys())
                                     - set(self.parameter_list.keys()))
            message = ""
            if missing_parameters:
                message += ("Missing parameters: " +
                            ", ".join(missing_parameters) +
                            ".")
            if unexpected_parameters:
                message += ("Unexpected parameters: " +
                            ", ".join(unexpected_parameters) +
                            ".")
            raise KeyError("The parameters passed are different from the "
                           "expected one." + message)

        # TODO: We should also check that all parameters is from the right
        # type.

        # We do the actual processing.
        self._process(sed, parameters)


def get_module(name):
    """Get a SED creation module from its name

    Parameters
    ----------
    module_name : string
        The name of the module we want to get the class. This name can be
        prefixed by anything using a dot, then the part before the dot is
        used to determine the module to load (e.g. 'dh2002_ir.1' will return
        the 'dh2002_ir' module).

    Returns
    -------
    a pcigale.sed.modules.Module instance
    """

    # Determine the real module name by removing the dotted prefix.
    module_name = name.split('.')[0]

    try:
        # TODO Find a better way to do dynamic import
        import_string = 'from . import ' + module_name + ' as module'
        exec import_string
        return module.Module(name=name)
    except ImportError:
        print('Module ' + module_name + ' does not exists!')
        raise
