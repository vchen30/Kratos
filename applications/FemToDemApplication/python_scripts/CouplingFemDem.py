from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainDEM_for_coupling as DEM
import MainFEM_for_coupling as FEM


# Main script of the coupled FEM-DEM App
class FEMDEM_Solution:

	def __init__(self):
		self.FEM_Solution = FEM.FEM_Solution()
		self.DEM_Solution = DEM.Solution()


	def Run(self):

		FEM_Solution.Initialize()

		FEM_Solution.RunMainTemporalLoop()

		FEM_Solution.Finalize()
