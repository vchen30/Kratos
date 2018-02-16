from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMeshing

KratosMultiphysics.CheckForPreviousImport()

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
            "maximum_sub_grids"        : 1,
            "level_set_parameters"     : {
                "scalar_variable"          : "FREE_SURFACE_ELEVATION",
                "gradient_variable"        : "FREE_SURFACE_GRADIENT"
            },
            "fix_contour_model_parts"  : [],
            "element_name"             : "ElementName",
            "condition_name"           : "ConditionName"
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

        self.Mesher = KratosMeshing.TriGenPFEMModeler()

    def ExecuteInitialize(self):
        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        self.node_neigh_search = KratosMultiphysics.FindNodalNeighboursProcess(self.model_part,9,18)
        self.elem_neigh_search  = KratosMultiphysics.FindElementalNeighboursProcess(self.model_part, 2, 10)
        self.cond_neigh_search  = KratosMultiphysics.FindConditionsNeighboursProcess(self.model_part,2, 10)

        (self.node_neigh_search).Execute()

        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.model_part)
        #~ self.find_nodal_h.Execute()

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
        #~ self.local_gradient.Execute()
        # Recalculate NODAL_H
        #~ self.find_nodal_h.Execute()

        for node in self.model_part.Nodes:
            node.SetsolutionStepValue(KratosMultiphysics.NODAL_H, 0, 0.1)

        h_factor = 0.1;
        alpha_shape = 1.2;
        node_erase_process = KratosMultiphysics.NodeEraseProcess(self.model_part);
        rem_nodes = False
        add_nodes = True
        (self.Mesher).ReGenerateMesh(self.element_name, self.condition_name, self.model_part, node_erase_process, rem_nodes, add_nodes, alpha_shape, h_factor)

        (self.fluid_neigh_finder).Execute()
        (self.elem_neigh_finder).Execute()
        (self.cond_neigh_finder).Execute()

    def __generate_submodelparts_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]

