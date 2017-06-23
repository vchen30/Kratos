
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *

def AssignMaterial(Properties):
    prop_id = 1;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStress2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
