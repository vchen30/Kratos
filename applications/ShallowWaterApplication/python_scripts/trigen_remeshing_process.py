from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMeshing
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

KratosMultiphysics.CheckForPreviousImport()

from json_utilities import *
import numpy as np
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TriGenRemeshingProcess(Model, settings["Parameters"])

class TriGenRemeshingProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"          : "MainModelPart",
            "maximum_sub_grids"        : 4,
            "level_set_parameters"     : {
                "scalar_variable"          : "FREE_SURFACE_ELEVATION",
                "gradient_variable"        : "FREE_SURFACE_GRADIENT"
            },
            "fix_contour_model_parts"  : [],
            "element_name"             : "ElementName",
            "condition_name"           : "ConditionName",
            "maximum_nodal_h"          : 0.1,
            "maximum_sub_grids"        : 4
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.dim = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.maximum_sub_grids = self.settings["maximum_sub_grids"].GetInt()
        
        self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_parameters"]["scalar_variable"].GetString() )
        self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_parameters"]["gradient_variable"].GetString() )

        self.element_name = self.settings["element_name"].GetString()
        self.condition_name = self.settings["condition_name"].GetString()

        self.maximum_nodal_h = self.settings["maximum_nodal_h"].GetDouble()
        self.maximum_sub_grids = self.settings["maximum_sub_grids"].GetInt()

        self.Mesher = KratosMeshing.TriGenPFEMModeler()

    def ExecuteInitialize(self):
        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        self.node_neigh_search = KratosMultiphysics.FindNodalNeighboursProcess(self.model_part,9,18)
        self.elem_neigh_search = KratosMultiphysics.FindElementalNeighboursProcess(self.model_part, 2, 10)
        self.cond_neigh_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.model_part,2, 10)

        (self.node_neigh_search).Execute()

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
        self._CalculateNodalH()

        h_factor = 0.1
        alpha_shape = 1.2
        node_erase_process = KratosMultiphysics.NodeEraseProcess(self.model_part)
        rem_nodes = True
        add_nodes = True
        (self.Mesher).ReGenerateMesh(self.element_name, self.condition_name, self.model_part, node_erase_process, rem_nodes, add_nodes, alpha_shape, h_factor)

        (self.node_neigh_search).Execute()
        (self.elem_neigh_search).Execute()
        (self.cond_neigh_search).Execute()

    def _CalculateNodalH(self):
        import statistics as stat
        gradient_norm = 0.0
        gradient_eta_values = []
        for node in self.model_part.Nodes:
            gradient_norm = np.linalg.norm(node.GetSolutionStepValue(self.gradient_variable))
            gradient_eta_values.append(gradient_norm)

        # Computing the percentage
        mean = stat.mean(gradient_eta_values)
        stdev = stat.stdev(gradient_eta_values)
    
        for node in self.model_part.Nodes:
            gradient_norm = np.linalg.norm(node.GetSolutionStepValue(self.gradient_variable))
            prob = _normpdf(gradient_norm, mean, stdev)
            new_nodal_h = self.maximum_nodal_h / (np.ceil(prob*self.maximum_sub_grids))
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H, new_nodal_h)

    def __generate_submodelparts_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]

def _linear_interpolation(x, x_list, y_list):
    tb = KratosMultiphysics.PiecewiseLinearTable()
    for i in range(len(x_list)):
        tb.AddRow(x_list[i], y_list[i])
        
    return tb.GetNearestValue(x)

def _normpdf(x, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data = read_external_json(dir_path+"/normal_distribution.json")
    z = (x-mean)/sd
    z_list = data["Z"]
    prob_list = data["Prob"]
    if (z > 0):
        prob = _linear_interpolation(z, z_list, prob_list)
    else:
        prob = 1.0 - _linear_interpolation(-z, z_list, prob_list)
    return prob


def _normvalf(prob, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data = read_external_json(dir_path+"/normal_distribution.json")
    z_list = data["Z"]
    prob_list = data["Prob"]
    if (prob >= 0.5):
        z = _linear_interpolation(prob, prob_list, z_list)
    else:
        z = - _linear_interpolation(1.0 - prob, prob_list, z_list)
    x = z * sd + mean
    return x

