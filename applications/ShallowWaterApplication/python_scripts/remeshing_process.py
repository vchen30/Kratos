from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMeshing

from json_utilities import *
import json
import os

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiGridProcess(Model, settings["Parameters"])

class MultiGridProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"          : "MainModelPart",
            "maximum_sub_grids"        : 1,
            "level_set_parameters"     : {
                "scalar_variable"          : "FREE_SURFACE_ELEVATION",
                "gradient_variable"        : "FREE_SURFACE_GRADIENT"
            },
            "fix_contour_model_parts"          : []
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.dim = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.maximum_sub_grids = self.settings["maximum_sub_grids"].GetInt()
        
        self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_parameters"]["scalar_variable"].GetString() )
        self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_parameters"]["gradient_variable"].GetString() )

    def ExecuteInitialize(self):
        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.model_part)
        self.find_nodal_h.Execute()

        # Calculate the parameters of automatic remeshing

        # NOTE: Add more model part if interested
        submodelpartslist = self.__generate_submodelparts_list_from_input(self.settings["fix_contour_model_parts"])

        for submodelpart in submodelpartslist:
            for node in submodelpart.Nodes:
                node.Set(KratosMultiphysics.BLOCKED, True)

        self._CreateGradientProcess()


    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        self._ExecuteRefinement()

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass


    def _CreateGradientProcess(self):
        # We compute the scalar value gradient
        self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.model_part, self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)

    def _ExecuteRefinement(self):
        # Calculate the gradient
        self.local_gradient.Execute()
        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

    def __generate_submodelparts_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]
