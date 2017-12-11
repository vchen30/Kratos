from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

#Base class to develop other solvers
class ModelManager(object):
    """The base class for solid mechanic model build process.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs.

    """
    def __init__(self, custom_settings):
        
        default_settings = KratosMultiphysics.Parameters("""
        {
           "model_name": "solid_domain",
           "dimension": 3,
           "bodies_list": [],
           "domain_parts_list": [],
           "processes_parts_list": [],
           "output_model_part_name": "output_domain", 
           "computing_model_part_name": "computing_domain",
           "import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
           },
           "variables":[],
           "dofs": []
         }
        """)

        
        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Defining the model_part
        self.main_model_part = KratosMultiphysics.ModelPart(self.settings["model_name"].GetString())   

        #TODO: replace this "model" for real one once available in kratos core
        self.model = {self.settings["model_name"].GetString() : self.main_model_part}
        
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, self.settings["dimension"].GetInt())

        computing_model_part = self.settings["computing_model_part_name"].GetString()
        self.main_model_part.CreateSubModelPart(computing_model_part)
        self.model.update({computing_model_part: self.main_model_part.GetSubModelPart(computing_model_part)})

        #output_model_part = self.settings["output_model_part_name"].GetString()
        #self.main_model_part.CreateSubModelPart(output_model_part)
        #self.model.update({output_model_part: self.main_model_part.GetSubModelPart(output_model_part)})

        
        # Legacy
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.settings["dimension"].GetInt())

        # Variables and Dofs settings
        self.nodal_variables = []
        self.dof_variables   = []
        self.dof_reactions   = [] 
        
        
    def ImportModel(self):

        self._add_variables()
                
        print("::[Model_Manager]:: Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["import_settings"]["input_filename"].GetString()
        
        if(self.settings["import_settings"]["input_type"].GetString() == "mdpa"):            
            # Import model part from mdpa file.
            print("   Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa ")
            sys.stdout.flush()
          
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            # print("   Finished reading model part from mdpa file ")
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            self._add_dofs()
            
        elif(self.settings["import_settings"]["input_type"].GetString() == "rest"):
            # Import model part from restart file.
            restart_path = os.path.join(problem_path, self.settings["import_settings"]["input_filename"].GetString() + "__" + self.settings["import_settings"]["input_file_label"].GetString() )
            if(os.path.exists(restart_path+".rest") == False):
                raise Exception("Restart file not found: " + restart_path + ".rest")
            print("   Loading Restart file: ", restart_path + ".rest ")
            # set serializer flag
            serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = KratosMultiphysics.Serializer(restart_path, serializer_flag)
            serializer.Load(self.main_model_part.Name, self.main_model_part)

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
            #I use it to rebuild the contact conditions.
            load_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] +1;
            self.main_model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step
            # print("   Finished loading model part from restart file ")            

        else:
            raise Exception("Other input options are not yet implemented.")

        
        print ("::[Model_Manager]:: Finished importing model part.")
            
    def ExportModel(self):
        name_out_file = self.settings["import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        # Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)
       

    def GetProcessInfo(self):
        return self.main_model_part.ProcessInfo

    def GetModel(self):
        return self.model
 
    def GetMainModelPart(self):
        return self.main_model_part
                
    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def GetOutputModelPart(self):
        #return self.main_model_part.GetSubModelPart(self.settings["output_model_part_name"].GetString())
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())
        
    def SaveRestart(self):
        pass #one should write the restart file here
    
    def SetVariables(self, variables):
        self.nodal_variables = self.nodal_variables + variables
    
    
    #### Model manager internal methods ####

    def _add_variables(self):

        self._set_variables()
        
        self.nodal_variables = self.nodal_variables + self.dof_variables + self.dof_reactions 

        self.nodal_variables = [self.nodal_variables[i] for i in range(0,len(self.nodal_variables)) if self.nodal_variables[i] != 'NOT_DEFINED']

        for variable in self.nodal_variables:            
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))
            #print(" Added variable ", KratosMultiphysics.KratosGlobals.GetVariable(variable),"(",variable,")")
            
        print("::[Model_Manager]:: General Variables ADDED")
                                                              
        
    def _add_dofs(self):
        AddDofsProcess = KratosSolid.AddDofsProcess(self.main_model_part, self.dof_variables, self.dof_reactions)
        AddDofsProcess.Execute()
                
        print("::[Model_Manager]:: DOF's ADDED")

    def _set_variables(self):
                
        # Add displacements
        self.dof_variables = self.dof_variables + ['DISPLACEMENT']
        self.dof_reactions = self.dof_reactions + ['REACTION'] 
        
        # Add dynamic variables
        self.nodal_variables = self.nodal_variables + ['VELOCITY','ACCELERATION']
        
        # Add specific variables for the problem conditions
        self.nodal_variables = self.nodal_variables + ['VOLUME_ACCELERATION','POSITIVE_FACE_PRESSURE','NEGATIVE_FACE_PRESSURE','POINT_LOAD','LINE_LOAD','SURFACE_LOAD']
        
        # Add nodal force variables for component wise calculation
        if( self.settings.Has("component_wise") ):
            if self.settings["component_wise"].GetBool():
                self.nodal_variables = self.nodal_variables + ['INTERNAL_FORCE','EXTERNAL_FORCE']
 
        # Add rotational variables
        if self._check_input_dof("ROTATION"):                    
            # Add specific variables for the problem (rotation dofs)
            self.dof_variables = self.dof_variables + ['ROTATION']
            self.dof_reactions = self.dof_reactions + ['TORQUE']

            self.nodal_variables = self.nodal_variables + ['ANGULAR_VELOCITY','ANGULAR_ACCELERATION']
            # Add specific variables for the problem conditions
            self.nodal_variables = self.nodal_variables + ['POINT_MOMENT']
            # Add large rotation variables
            self.nodal_variables = self.nodal_variables + ['STEP_DISPLACEMENT','STEP_ROTATION','DELTA_ROTATION']

            
        # Add pressure variables
        if self._check_input_dof("PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['PRESSURE']
            self.dof_reactions = self.dof_reactions + ['PRESSURE_REACTION']
            #if not self.settings["stabilization_factor"].IsNull():
            #    self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Add contat variables
        if self._check_input_dof("LAGRANGE_MULTIPLIER"):
            # Add specific variables for the problem (contact dofs)
            self.dof_variables = self.dof_variables + ['LAGRANGE_MULTIPLIER_NORMAL']
            self.dof_reactions = self.dof_reactions + ['LAGRANGE_MULTIPLIER_NORMAL_REACTION']        
            
    def _check_input_dof(self, variable):
        dofs_list = self.settings["dofs"]
        for i in range(0, dofs_list.size() ):
            if dofs_list[i].GetString() == variable:
                return True
        return False
    
    def _execute_after_reading(self):
        # The computing_model_part is labeled 'KratosMultiphysics.ACTIVE' flag (in order to recover it)
        self.computing_model_part_name = "computing_domain" 

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["domain_parts_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_parts_list"])
       
        if( self.settings.Has("bodies_list") ):
            params.AddValue("bodies_list",self.settings["bodies_list"])

        # CheckAndPrepareModelProcess creates the computating_model_part
        import check_and_prepare_model_process_solid
        check_and_prepare_model_process_solid.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Get the list of the model_part's in the object Model
        for i in range(self.settings["domain_parts_list"].size()):
            part_name = self.settings["domain_parts_list"][i].GetString()
            if( self.main_model_part.HasSubModelPart(part_name) ):
                self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
                
        for i in range(self.settings["processes_parts_list"].size()):
            part_name = self.settings["processes_parts_list"][i].GetString()
            if( self.main_model_part.HasSubModelPart(part_name) ):
                self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
