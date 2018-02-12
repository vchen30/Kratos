from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainDEM_for_coupling as DEM
import MainFEM_for_coupling as FEM
import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.FemToDemApplication   as KratosFemDem
import CouplingFemDem
import math
import os

def Wait():
	input("Press Something")

# Main script of the coupled FEM-DEM Application 3D
class FEMDEM3D_Solution(CouplingFemDem.FEMDEM_Solution):

#============================================================================================================================
	def Info(self):
		print("Coupling of the 3D FEMDEM App")

#============================================================================================================================
	def Initialize(self):
		self.FEM_Solution.Initialize()
		self.DEM_Solution.Initialize()

		self.SpheresModelPart = self.DEM_Solution.spheres_model_part
		self.DEMParameters = self.DEM_Solution.DEM_parameters
		self.DEMProperties = self.SpheresModelPart.GetProperties()[1]

		self.ParticleCreatorDestructor = PCD.FemDemParticleCreatorDestructor(self.SpheresModelPart,
                                                                             self.DEMProperties,
                                                                             self.DEMParameters)

		self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
		self.nodal_neighbour_finder.Execute()

#============================================================================================================================
	def SolveSolutionStep(self):
		
		# Function to perform the coupling FEM <-> DEM
		self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

		#### SOLVE FEM #########################################
		self.FEM_Solution.solver.Solve()
		########################################################

		#self.GenerateDEM()            # we create the new DEM of this time step
		self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
		#self.CheckForPossibleIndentations()
		#self.CheckInactiveNodes()
		#self.UpdateDEMVariables()     # We update coordinates, displ and velocities of the DEM according to FEM


		# TESTING
		#self.FEM_Solution.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.NODAL_STRESS_VECTOR)
		#KratosFemDem.StressToNodesProcess(self.FEM_Solution.main_model_part).Execute()

		#print(self.FEM_Solution.main_model_part.GetNode(25).GetSolutionStepValue(KratosFemDem.NODAL_STRESS_VECTOR))
		#print(self.FEM_Solution.main_model_part.GetNode(25).GetSolutionStepValue(KratosFemDem.EQUIVALENT_NODAL_STRESS))
		# print(self.FEM_Solution.main_model_part.GetNode(25).GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
		#Wait()
		# TESTING

		self.DEM_Solution.InitializeTimeStep()

		self.DEM_Solution.time = self.FEM_Solution.time
		self.DEM_Solution.step = self.FEM_Solution.step

		self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts, self.DEM_Solution.time,self.DEM_Solution.dt,self.DEM_Solution.step)
		self.DEM_Solution.BeforeSolveOperations(self.DEM_Solution.time)

		#### SOLVE DEM #########################################
		self.DEM_Solution.solver.Solve()
		########################################################
		self.DEM_Solution.AfterSolveOperations()

		self.DEM_Solution.DEMFEMProcedures.MoveAllMeshes(self.DEM_Solution.all_model_parts, self.DEM_Solution.time, self.DEM_Solution.dt)

		self.UpdateDEMVariables() # to print DEM with the FEM coordinates

		# DEM GiD print output
		if self.DEM_Solution.step == 1: # always print the 1st step
			self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
			self.DEM_Solution.time_old_print = self.DEM_Solution.time
		else:
			time_to_print = self.DEM_Solution.time - self.DEM_Solution.time_old_print

			if (self.DEM_Solution.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.DEM_Solution.dt):

				self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
				self.DEM_Solution.time_old_print = self.DEM_Solution.time


		self.DEM_Solution.FinalizeTimeStep(self.DEM_Solution.time)

		# Transfer the contact forces of the DEM to the FEM nodes
		self.TransferNodalForcesToFEM()

		self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

		# Update Coupled Postprocess file for Gid (post.lst)
		self.WritePostListFile()

		# Print required info
		self.PrintPlotsFiles()