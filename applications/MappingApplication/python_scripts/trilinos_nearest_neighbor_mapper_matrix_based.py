from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as MappingApplication
import mapper_matrix_based

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMapper(main_model_part_origin, main_model_part_destination, custom_settings): 
    return TrilinosNearestNeighborMapperMatrixBased(main_model_part_origin, main_model_part_destination, custom_settings)


class TrilinosNearestNeighborMapperMatrixBased(mapper_matrix_based.MapperMatrixBased):
    """
    This is the ... Mapper
    """
    def __init__(self, main_model_part_origin, main_model_part_destination, custom_settings):
        mapper_settings = KratosMultiphysics.Parameters("""
        {
            "some_settings" : 444
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, mapper_settings)

        # Construct the base mapper.
        super(TrilinosNearestNeighborMapperMatrixBased, self).__init__(main_model_part_origin, 
                                           main_model_part_destination, 
                                           custom_settings)
