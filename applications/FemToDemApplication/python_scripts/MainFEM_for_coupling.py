from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainFemDem

# Python script created to modify the existing one due to the coupling of the DEM app

class FEM_for_coupling_Solution(MainFemDem.Solution):

	def Info(self):
		print("FEM part of the FEMDEM application")