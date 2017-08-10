from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as MappingApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMapper(main_model_part_origin, main_model_part_destination, custom_settings): # TODO why is this function needed in the baseclass?
    return MapperBase(main_model_part_origin, main_model_part_destination, custom_settings)


class MapperBase(object):
    """The python wrapper of the mappers.
    It returns a mapper defined in the custom_settings
    """
    def __init__(self, main_model_part_origin, main_model_part_destination, custom_settings): 
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "mapper_type"                           : "",
            "interface_submodel_part_origin"        : "",
            "interface_submodel_part_destination"   : "",
            "parallel_type"                         : "OpenMP"
            "search_radius"                         : -1.0,
            "search_iterations"                     : 3,
            "echo_level"                            : 0
        } )""" )

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self._read_model_parts(model_part_origin, main_model_part_destination)

        self._do_some_more_checks # TODO implement this


        self._create_mapper()



    #### Public functions ####

    def Map(variable_origin, variable_destination):
        # Create the mapper only here?
        self._get_mapper().Map()


    def InverseMap(variable_origin, variable_destination):
        mapper = self._get_mapper()
        if hasattr(mapper, _inverse_mapper):
            mapper._InitializeInverseMapper()
        self._get_inverse_mapper().Map()
    

    def UpdateInterface():
        self._get_mapper().UpdateInterface()
        if hasattr(mapper, _inverse_mapper):
            self._get_inverse_mapper().UpdateInterface()



    #### Specific internal functions ####

    def validate_and_transfer_matching_settings(self, origin_settings, destination_settings):
        """Transfer matching settings from origin to destination.

        !This function was taken from the base mechanical structural solver (StructuralMechanicsApplication)!

        If a name in origin matches a name in destination, then the setting is
        validated against the destination.

        The typical use is for validating and extracting settings in derived classes:

        class A:
            def __init__(self, model_part, a_settings):
                default_a_settings = Parameters('''{
                    ...
                }''')
                a_settings.ValidateAndAssignDefaults(default_a_settings)
        class B(A):
            def __init__(self, model_part, custom_settings):
                b_settings = Parameters('''{
                    ...
                }''') # Here the settings contain default values.
                self.validate_and_transfer_matching_settings(custom_settings, b_settings)
                super().__init__(model_part, custom_settings)
        """
        for name, dest_value in destination_settings.items():
            if origin_settings.Has(name): # Validate and transfer value.
                orig_value = origin_settings[name]
                if dest_value.IsDouble() and orig_value.IsDouble():
                    destination_settings[name].SetDouble(origin_settings[name].GetDouble())
                elif dest_value.IsInt() and orig_value.IsInt():
                    destination_settings[name].SetInt(origin_settings[name].GetInt())
                elif dest_value.IsBool() and orig_value.IsBool():
                    destination_settings[name].SetBool(origin_settings[name].GetBool())
                elif dest_value.IsString() and orig_value.IsString():
                    destination_settings[name].SetString(origin_settings[name].GetString())
                elif dest_value.IsArray() and orig_value.IsArray():
                    if dest_value.size() != orig_value.size():
                        raise Exception('len("' + name + '") != ' + str(dest_value.size()))
                    for i in range(dest_value.size()):
                        if dest_value[i].IsDouble() and orig_value[i].IsDouble():
                            dest_value[i].SetDouble(orig_value[i].GetDouble())
                        elif dest_value[i].IsInt() and orig_value[i].IsInt():
                            dest_value[i].SetInt(orig_value[i].GetInt())
                        elif dest_value[i].IsBool() and orig_value[i].IsBool():
                            dest_value[i].SetBool(orig_value[i].GetBool())
                        elif dest_value[i].IsString() and orig_value[i].IsString():
                            dest_value[i].SetString(orig_value[i].GetString())
                        elif dest_value[i].IsSubParameter() and orig_value[i].IsSubParameter():
                            self.validate_and_transfer_matching_settings(orig_value[i], dest_value[i])
                            if len(orig_value[i].items()) != 0:
                                raise Exception('Json settings not found in default settings: ' + orig_value[i].PrettyPrintJsonString())
                        else:
                            raise Exception('Unsupported parameter type.')
                elif dest_value.IsSubParameter() and orig_value.IsSubParameter():
                    self.validate_and_transfer_matching_settings(orig_value, dest_value)
                    if len(orig_value.items()) != 0:
                        raise Exception('Json settings not found in default settings: ' + orig_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type.')
                origin_settings.RemoveValue(name)
    

    def get_mapper_communicator(self):
        if not hasattr(self, '_mapper_communicator'):
            self._mapper_communicator = self._create_mapper_communicator()
        return self._mapper_communicator

    

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


    def _create_mapper_communicator(self):
        if self.settings["parallel_type"].GetString() = "OpenMP":
            mapper_communicator = MappingApplication.MapperCommunicator()
        else: # MPI
            mapper_communicator = MappingApplication.MapperMPICommunicator()
        return mapper_communicator


    def _create_mapper(self):
        """
        Create the specific mapper. This function has to be implemented in the derived classes
        """
        raise Exception("Mapper creation must be implemented in the derived class.")
