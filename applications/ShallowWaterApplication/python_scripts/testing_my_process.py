from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
#~ import KratosMultiphysics.MeshingApplication as KratosMesshing
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TestingMyProcess(Model, settings["Parameters"])

class TestingMyProcess(KratosMultiphysics.Process):
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

        #~ self.Mesher = KratosMeshing.TriGenPFEMModeler()
