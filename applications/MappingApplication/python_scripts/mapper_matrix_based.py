from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as MappingApplication
import mapper_base

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMapper(main_model_part_origin, main_model_part_destination, custom_settings): # TODO why is this function needed in the baseclass?
    return MapperMatrixBased(main_model_part_origin, main_model_part_destination, custom_settings)


class MapperMatrixBased(mapper_base.MapperBase):
    """The python wrapper of the mappers.
    It returns a mapper defined in the custom_settings
    """
    def __init__(self, main_model_part_origin, main_model_part_destination, custom_settings): 

        self._preprocess_interface() # TODO I don't think that this works, most likely not the derived classes method is called!

        super(MapperMatrixBased, self).__init__(main_model_part_origin, 
                                                main_model_part_destination, 
                                                custom_settings)
    

    #### Specific internal functions ####

    def get_mapper_scheme(self):
        if not hasattr(self, '_mapper_scheme'):
            self._mapper_scheme = self._create_mapper_scheme()
        return self._mapper_scheme


    def get_mapping_matrix_builder(self):
        if not hasattr(self, '_mapping_matrix_builder'):
            self._mapping_matrix_builder = self._create_mapping_matrix_builder()
        return self._mapping_matrix_builder


    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver


    def get_builder_and_solver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver
    


    #### Private functions ####

    def _read_model_parts(self, main_model_part_origin, main_model_part_destination):
        echo_level = self.settings["echo_level"].GetInt()
        comm_rank_destination = main_model_part_destination.GetCommunicator().MyPID()

        # Get the ModelParts
        name_submodel_part_origin = self.settings["interface_submodel_part_origin"].GetString()
        name_submodel_part_destination = self.settings["interface_submodel_part_destination"].GetString()

        if name_submodel_part_origin == "":
            self.model_part_origin = main_model_part_origin
            if (echo_level >= 3 and comm_rank_destination == 0):
                print("Mapper: Main ModelPart used for Origin-ModelPart")
            
        else:
            self.model_part_origin = main_model_part_origin.GetSubModelPart(name_submodel_part_origin)
            if (echo_level >= 3 and comm_rank_destination == 0):
                print("Mapper: SubModelPart used for Origin-ModelPart")

        if name_submodel_part_destination == "":
            self.model_part_destination = main_model_part_destination
            if (echo_level >= 3 and comm_rank_destination == 0):
                print("Mapper: Main ModelPart used for Destination-ModelPart")
        else:
            self.model_part_destination = main_model_part_destination.GetSubModelPart(name_submodel_part_destination)
            if (echo_level >= 3 and comm_rank_destination == 0):
                print("Mapper: SubModelPart used for Destination-ModelPart")


    def _create_mapper_scheme(self):
        return MappingApplication.MapperScheme()


    def _create_mapping_matrix_builder(self):
        if self.settings["parallel_type"].GetString() = "OpenMP":
            mapping_matrix_builder = MappingApplication.MappingMatrixBuilder()
        else: # MPI
            mapping_matrix_builder = MappingApplication.TrilinosMappingMatrixBuilder()
        return mapping_matrix_builder
    

    def _create_linear_solver(self):
        if self.settings["parallel_type"].GetString() = "OpenMP":
            import linear_solver_factory as linear_solver_factory
        else: # MPI
            import trilinos_linear_solver_factory as linear_solver_factory

        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver


    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["parallel_type"].GetString() = "OpenMP":
            if(self.settings["block_builder"].GetBool() == True):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        else: # MPI            
            trilinos = __import__(KratosMultiphysics.TrilinosApplication) # TODO do this every time?

            if(self.settings["block_builder"].GetBool() == True):
                builder_and_solver = trilinos.TrilinosBlockBuilderAndSolver(linear_solver)
            else:
                builder_and_solver = trilinos.TrilinosResidualBasedEliminationBuilderAndSolver(linear_solver)

        return builder_and_solver


    def _create_mapper_strategy(self):
        mapper_scheme = self.get_mapper_scheme()
        if mortar... :
            builder_and_solver = self.get_builder_and_solver()

        if mortar... :
            mapper_strategy = MappingApplication.MapperStrategy(self.model_part_origin, 
                                                 self.model_part_destination,
                                                 mapper_scheme,
                                                 mapping_matrix_builder,
                                                 builder_and_solver
                                                 )
        else:
            mapper_strategy = MappingApplication.MapperStrategy(self.model_part_origin, 
                                                 self.model_part_destination,
                                                 mapper_scheme,
                                                 mapping_matrix_builder
                                                 )
        return mapper_strategy


    def _preprocess_interface(self):
        """
        Create the Interface ModelPart. This function has to be implemented in the derived classes
        """
        raise Exception("Mapper creation must be implemented in the derived class.")
