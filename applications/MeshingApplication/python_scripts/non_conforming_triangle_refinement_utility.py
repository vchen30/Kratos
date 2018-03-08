from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

class NonConformingTriangleRefinementUtility:
    def __init__(self, origin_model_part, destination_model_part):
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part
    
    