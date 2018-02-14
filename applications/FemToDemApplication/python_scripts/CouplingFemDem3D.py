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

		self.GenerateDEM()            # we create the new DEM of this time step
		self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
		self.CheckForPossibleIndentations()
		self.CheckInactiveNodes()
		self.UpdateDEMVariables()     # We update coordinates, displ and velocities of the DEM according to FEM


		# TESTING
		#num = 3
		#for i in range(2,num+1):
		#	print(i)
		#print(list(range(2,num))
		#M = [[1,2],[3,4]]
		#print(M[0][0])
		#Wait()
		#M = KratosMultiphysics.Matrix()
		#M.resize(2,2)
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

#============================================================================================================================
	def GenerateDEM(self): # 3D version

		FEM_elements = self.FEM_Solution.main_model_part.Elements

		# Loop Over Elements to find the INACTIVE ones and generate the DEM only once
		for Element in FEM_elements:

			is_active     = True
			DEM_Generated = Element.GetValue(KratosFemDem.DEM_GENERATED)

			if Element.IsDefined(KratosMultiphysics.ACTIVE):
				is_active = Element.Is(KratosMultiphysics.ACTIVE)

			NumberOfDEM = 0         # Number of nodes with DEM Associated
			for node in range(0, 4): # Loop over nodes of the FE
				Node = Element.GetNodes()[node]
				if Node.GetValue(KratosFemDem.IS_DEM) == True:
					NumberOfDEM += 1

			if is_active == False and DEM_Generated == False: # Let's generate the remaining DEM of this FE

				dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
				dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
				dist03  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[3])
				dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])
				dist13  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[3])
				dist23  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[2], Element.GetNodes()[3])

				print("Elemento: ",Element.Id ," tiene ", NumberOfDEM, " DEMs")

				# --------------------- 1ST SCENARIO -----------------------------
				if NumberOfDEM == 0: # we must create 4 DEM

					# Look to the Node 1 ---------------------------------------------
					Radius1 = self.GetMinimumValue(dist01, dist02, dist03) * 0.5
					Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[0])
					Id1 = Element.GetNodes()[0].Id
					self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, Radius1, Id1)
					Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
					Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, Radius1)

					# Look to the Node 2 ---------------------------------------------
					Radius2 = self.GetMinimumValue(dist01-Radius1, dist12*0.5, dist13*0.5) 
					Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[1])
					Id2 = Element.GetNodes()[1].Id
					self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, Radius2, Id2)
					Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
					Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, Radius2)

					# look to the node 3 ---------------------------------------------
					Radius3 = self.GetMinimumValue(dist02-Radius1, dist12-Radius2, dist23*0.5)
					Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[2])
					Id3 = Element.GetNodes()[2].Id
					self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, Radius3, Id3)
					Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
					Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, Radius3)

					# look to the node 4 ---------------------------------------------
					Radius4 = self.GetMinimumValue(dist03-Radius1, dist13-Radius2, dist23-Radius3)
					Coordinates4 = self.GetNodeCoordinates(Element.GetNodes()[3])
					Id4 = Element.GetNodes()[3].Id
					self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates4, Radius4, Id4)
					Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
					Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, Radius4)

					# DEM generated for this Element
					Element.SetValue(KratosFemDem.DEM_GENERATED, True)

				# --------------------- 2ND SCENARIO -----------------------------
				elif NumberOfDEM == 3: # we must create 1 DEM

					localId = 0 # Local Id of the node without DEM
					for index in range(0, 4):
						if Element.GetNodes()[index].GetValue(KratosFemDem.IS_DEM) == False:
							localId = index
							break

					Coordinates = self.GetNodeCoordinates(Element.GetNodes()[localId])
					Id = Element.GetNodes()[localId].Id

					if localId == 0:
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)
						R0 = self.GetMinimumValue(dist01-R1, dist02-R2, dist03-R3)
						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R0, Id)
						Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R0)

					elif localId == 1:
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)
						R1 = self.GetMinimumValue(dist01-R0, dist12-R2, dist13-R3)
						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R1, Id)
						Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R1)

					elif localId == 2:
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)
						R2 = self.GetMinimumValue(dist02-R0, dist12-R1, dist23-R3)
						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R2, Id)
						Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R2)

					elif localId == 3:
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
						R3 = self.GetMinimumValue(dist03-R0, dist13-R1, dist23-R2)
						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R3, Id)
						Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R3)

					# DEM generated for this Element
					Element.SetValue(KratosFemDem.DEM_GENERATED, True)

				# --------------------- 3RD SCENARIO -----------------------------
				elif NumberOfDEM == 1: # we must create 3 DEM

					localId = 0 # Local Id of the node with DEM
					for index in range(0,4):
						if Element.GetNodes()[index].GetValue(KratosFemDem.IS_DEM) == True:
							localId = index
							break

					RadiusOfDem = Element.GetNodes()[localId].GetValue(KratosMultiphysics.RADIUS)

					if localId == 0:
						R1 = self.GetMinimumValue(dist01-RadiusOfDem, dist12*0.5, dist13*0.5)
						R2 = self.GetMinimumValue(dist02-RadiusOfDem, dist12-R1,dist23*0.5)
						R3 = self.GetMinimumValue(dist03-RadiusOfDem, dist13-R1, dist23-R2)

						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id1 = Element.GetNodes()[1].Id
						Id2 = Element.GetNodes()[2].Id
						Id3 = Element.GetNodes()[3].Id

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					elif localId == 1:
						R0 = self.GetMinimumValue(dist01-RadiusOfDem, dist02*0.5, dist03*0.5)
						R2 = self.GetMinimumValue(dist12-RadiusOfDem, dist02-R0, dist23*0.5)
						R3 = self.GetMinimumValue(dist13-RadiusOfDem, dist03-R0, dist23-R2)

						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id0 = Element.GetNodes()[0].Id
						Id2 = Element.GetNodes()[2].Id
						Id3 = Element.GetNodes()[3].Id

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					elif localId == 2:
						R0 = self.GetMinimumValue(dist02-RadiusOfDem, dist01*0.5, dist03*0.5)
						R1 = self.GetMinimumValue(dist12-RadiusOfDem, dist01-R0, dist13*0.5)
						R3 = self.GetMinimumValue(dist23-RadiusOfDem, dist03-R0, dist13-R1)

						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id0 = Element.GetNodes()[0].Id
						Id1 = Element.GetNodes()[2].Id
						Id3 = Element.GetNodes()[3].Id

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					elif localId == 3:
						R0 = self.GetMinimumValue(dist03-RadiusOfDem, dist01*0.5, dist02*0.5)
						R1 = self.GetMinimumValue(dist13-RadiusOfDem, dist01-R0, dist12*0.5)
						R2 = self.GetMinimumValue(dist23-RadiusOfDem, dist02-R0, dist12-R1)

						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Id0 = Element.GetNodes()[0].Id
						Id1 = Element.GetNodes()[2].Id
						Id2 = Element.GetNodes()[2].Id

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

					# DEM generated for this Element
					Element.SetValue(KratosFemDem.DEM_GENERATED, True)	

				# --------------------- 4RD SCENARIO -----------------------------
				elif NumberOfDEM == 2: # we must create 2 DEM

					NodesWithDEMId = []

					for index in range(0,4):
						if Element.GetNodes()[index].GetValue(KratosFemDem.IS_DEM) == True:
							NodesWithDEMId.append(index)

					if NodesWithDEMId[0] == 0 and NodesWithDEMId[1] == 1: # edge 0-1

						# Radius of existing DEM
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)

						# New DEM info
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id2 = Element.GetNodes()[2].Id
						Id3 = Element.GetNodes()[3].Id

						R2 = self.GetMinimumValue(dist02-R2, dist01-R1, dist23*0.5)
						R3 = self.GetMinimumValue(dist03-R0, dist13-R1, dist23-R2)
							
						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					elif NodesWithDEMId[0] == 0 and NodesWithDEMId[1] == 2: # edge 0-2

						# Radius of existing DEM
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)

						# New DEM info
						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id1 = Element.GetNodes()[1].Id
						Id3 = Element.GetNodes()[3].Id

						R1 = self.GetMinimumValue(dist01-R0, dist12-R2, dist13*0.5)
						R3 = self.GetMinimumValue(dist03-R0, dist23-R2, dist13-R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					elif NodesWithDEMId[0] == 0 and NodesWithDEMId[1] == 3: # edge 0-3

						# Radius of existing DEM
						R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)

						# New DEM info
						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Id1 = Element.GetNodes()[1].Id
						Id2 = Element.GetNodes()[2].Id

						R1 = self.GetMinimumValue(dist01-R0, dist13-R3, dist12*0.5)
						R2 = self.GetMinimumValue(dist02-R0, dist23-R3, dist12-R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

					elif NodesWithDEMId[0] == 1 and NodesWithDEMId[1] == 3: # edge 1-3

						# Radius of existing DEM
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)

						# New DEM info
						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
						Id0 = Element.GetNodes()[0].Id
						Id2 = Element.GetNodes()[2].Id

						R0 = self.GetMinimumValue(dist01-R1, dist03-R3, dist02*0.5)
						R2 = self.GetMinimumValue(dist12-R1, dist23-R3, dist02-R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
						Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

					elif NodesWithDEMId[0] == 2 and NodesWithDEMId[1] == 3: # edge 2-3

						# Radius of existing DEM
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
						R3 = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)

						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
						Id0 = Element.GetNodes()[0].Id
						Id1 = Element.GetNodes()[1].Id

						R0 = self.GetMinimumValue(dist02-R2, dist03-R3, dist01*0.5)
						R1 = self.GetMinimumValue(dist12-R2, dist13-R3, dist01-R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
						Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

					elif NodesWithDEMId[0] == 1 and NodesWithDEMId[1] == 2: # edge 1-2

						# Radius of existing DEM
						R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
						R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)

						Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
						Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[3])
						Id0 = Element.GetNodes()[0].Id
						Id3 = Element.GetNodes()[3].Id	

						R0 = self.GetMinimumValue(dist01-R1, dist02-R2, dist03*0.5)
						R3 = self.GetMinimumValue(dist13-R1, dist23-R2, dist03-R0)			

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
						Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

						self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, R3, Id3)
						Element.GetNodes()[3].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					# DEM generated for this Element
					Element.SetValue(KratosFemDem.DEM_GENERATED, True)

				# --------------------- 5RD SCENARIO -----------------------------
				elif NumberOfDEM == 4: # we avoid posible indentations

					R0  = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
					R1  = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
					R2  = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
					R3  = Element.GetNodes()[3].GetValue(KratosMultiphysics.RADIUS)
					Id0 = Element.GetNodes()[0].Id
					Id1 = Element.GetNodes()[1].Id
					Id2 = Element.GetNodes()[2].Id
					Id3 = Element.GetNodes()[3].Id

					# Check the 6 edges of the element
					if R0 + R1 > dist01:
						R0 = self.GetMinimumValue(R0, dist01*0.5, 1000000)
						R1 = dist01 - R0

						# assign the new radius to the DEM nodes
						self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
						self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

					if R0 + R2 > dist02:
						R0 = self.GetMinimumValue(R0, 0.5*dist02, 1000000)
						R2 = dist02 - R0

						# assign the new radius to the DEM nodes 
						self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
						self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

					if R0 + R3 > dist03:
						R0 = self.GetMinimumValue(R0, 0.5*dist03, 1000000)
						R3 = dist03 - R0

						# assign the new radius to the DEM nodes 
						self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
						self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
						Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					if R1 + R2 > dist12:
						R1 = self.GetMinimumValue(R1, 0.5*dist12, 1000000)
						R2 = dist12 - R1

						# assign the new radius to the DEM nodes
						self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
						self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

					if R1 + R3 > dist13:
						R1 = self.GetMinimumValue(R1, 0.5*dist13, 1000000)
						R3 = dist13 - R1

						# assign the new radius to the DEM nodes
						self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
						self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
						Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					if R2 + R3 > dist23:
						R2 = self.GetMinimumValue(R2, 0.5*dist23, 1000000)
						R3 = dist23 - R2

						# assign the new radius to the DEM nodes
						self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
						self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
						Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
						Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

					# DEM generated for this Element
					Element.SetValue(KratosFemDem.DEM_GENERATED, True)



#============================================================================================================================
	def CalculateDistanceBetweenNodes(self, Node1, Node2):

		X1 = Node1.X
		X2 = Node2.X
		Y1 = Node1.Y
		Y2 = Node2.Y
		Z1 = Node1.Z
		Z2 = Node2.Z
		return math.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2) + (Z1 - Z2) * (Z1 - Z2))

#============================================================================================================================
	def GetMinimumValue(self, val1, val2, val3):
		res = val1
		if val2 < val1:
			res = val2

		if val3 < res:
			res = val3

		return res

#============================================================================================================================
	def CheckForPossibleIndentations(self): # Verifies if an element has indentations between its DEM

		FEM_elements = self.FEM_Solution.main_model_part.Elements

		for Element in FEM_elements:

			is_active     = True
			DEM_Generated = Element.GetValue(KratosFemDem.DEM_GENERATED)

			if Element.IsDefined(KratosMultiphysics.ACTIVE):
				is_active = Element.Is(KratosMultiphysics.ACTIVE)

			NumberOfDEM = 0              # Number of nodes with DEM Associated
			NodesWithDEMLocalId = []     # Local Id of the nodes with DEM
			Radius = []
			Ids = []

			for index in range(0, 4): # Loop over nodes of the FE
				Node = Element.GetNodes()[index]
				if Node.GetValue(KratosFemDem.IS_DEM) == True:
					NumberOfDEM += 1
					NodesWithDEMLocalId.append(index)
					Radius.append(Node.GetValue(KratosMultiphysics.RADIUS))
					Ids.append(Node.Id)

			if NumberOfDEM > 1 and is_active == True and DEM_Generated == False:  # Case in which the DEM have been generated by its neighbours

				dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
				dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
				dist03  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[3])
				dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])
				dist13  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[3])
				dist23  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[2], Element.GetNodes()[3])

				Distances = [[0, dist01, dist02, dist03],
							 [dist01, 0, dist12, dist13],
							 [dist02, dist12, 0, dist23],
							 [dist03, dist13, dist23, 0]]

				for index in range(2, NumberOfDEM + 1):

					R1 = Radius[index - 2]
					R2 = Radius[index - 1]

					Id1 = Ids[index - 2]
					Id2 = Ids[index - 1]

					LocalId1 = NodesWithDEMLocalId[index - 2]
					LocalId2 = NodesWithDEMLocalId[index - 1]

					Distance = Distances[LocalId1][LocalId2]

					if R1 + R2 > Distance:

						R1 = self.GetMinimumValue(R1, 0.5*Distance, 1000000)
						R2 = Distance - R1
						self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
						self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
						Element.GetNodes()[LocalId1].SetValue(KratosMultiphysics.RADIUS, R1)
						Element.GetNodes()[LocalId2].SetValue(KratosMultiphysics.RADIUS, R2)


					if NumberOfDEM >= 3 and index == 3: # we add a connection

						R1 = Radius[index - 3]
						R2 = Radius[index - 1]

						Id1 = Ids[index - 3]
						Id2 = Ids[index - 1]

						LocalId1 = NodesWithDEMLocalId[index - 3]
						LocalId2 = NodesWithDEMLocalId[index - 1]

						Distance = Distances[LocalId1][LocalId2]

						if R1 + R2 > Distance:

							R1 = self.GetMinimumValue(R1, 0.5*Distance, 1000000)
							R2 = Distance - R1
							self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
							self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
							Element.GetNodes()[LocalId1].SetValue(KratosMultiphysics.RADIUS, R1)
							Element.GetNodes()[LocalId2].SetValue(KratosMultiphysics.RADIUS, R2)

					if NumberOfDEM == 4 and index == 4: # we add a connection

						R1 = Radius[index - 3]
						R2 = Radius[index - 1]

						Id1 = Ids[index - 3]
						Id2 = Ids[index - 1]

						LocalId1 = NodesWithDEMLocalId[index - 3]
						LocalId2 = NodesWithDEMLocalId[index - 1]

						Distance = Distances[LocalId1][LocalId2]

						if R1 + R2 > Distance:

							R1 = self.GetMinimumValue(R1, 0.5*Distance, 1000000)
							R2 = Distance - R1
							self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
							self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
							Element.GetNodes()[LocalId1].SetValue(KratosMultiphysics.RADIUS, R1)
							Element.GetNodes()[LocalId2].SetValue(KratosMultiphysics.RADIUS, R2)

						# DEM generated for this Element -> 4 DEM
						Element.SetValue(KratosFemDem.DEM_GENERATED, True)








				'''
				# Check the 6 edges of the element
				if R0 + R1 > dist01:
					R0 = self.GetMinimumValue(R0, dist01*0.5)
					R1 = dist01 - R0

					# assign the new radius to the DEM nodes
					self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
					self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
					Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
					Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

				if R0 + R2 > dist02:
					R0 = self.GetMinimumValue(R0, 0.5*dist02)
					R2 = dist02 - R0

					# assign the new radius to the DEM nodes 
					self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
					self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
					Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
					Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

				if R0 + R3 > dist03:
					R0 = self.GetMinimumValue(R0, 0.5*dist03)
					R3 = dist03 - R0

					# assign the new radius to the DEM nodes 
					self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
					self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
					Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
					Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

				if R1 + R2 > dist12:
					R1 = self.GetMinimumValue(R1, 0.5*dist12)
					R2 = dist12 - R1

					# assign the new radius to the DEM nodes
					self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
					self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
					Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
					Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

				if R1 + R3 > dist13:
					R1 = self.GetMinimumValue(R1, 0.5*dist13)
					R3 = dist13 - R1

					# assign the new radius to the DEM nodes
					self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
					self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
					Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)
					Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

				if R2 + R3 > dist23:
					R2 = self.GetMinimumValue(R2, 0.5*dist23)
					R3 = dist23 - R2

					# assign the new radius to the DEM nodes
					self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
					self.DEM_Solution.spheres_model_part.GetNode(Id3).SetSolutionStepValue(KratosMultiphysics.RADIUS, R3)
					Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
					Element.GetNodes()[3].SetValue(KratosMultiphysics.RADIUS, R3)

				# DEM generated for this Element
				Element.SetValue(KratosFemDem.DEM_GENERATED, True)'''

#============================================================================================================================
	def CheckInactiveNodes(self):

		FEM_Elements = self.FEM_Solution.main_model_part.Elements
		FEM_Nodes    = self.FEM_Solution.main_model_part.Nodes
		DEM_Nodes    = self.SpheresModelPart

		for Element in FEM_Elements:

			is_active = True

			if Element.IsDefined(KratosMultiphysics.ACTIVE):
				is_active = Element.Is(KratosMultiphysics.ACTIVE)

			if is_active == True:

				for i in range(0, 4): # Loop over nodes of the element
					node = Element.GetNodes()[i]
					NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
					NumberOfActiveElements += 1
					node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, NumberOfActiveElements)


		NumberOfActiveElements = 0	
		for node in FEM_Nodes:
			NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)

			if NumberOfActiveElements == 0 and node.GetValue(KratosFemDem.INACTIVE_NODE) == False:

				Id = node.Id
				DEMnode = self.SpheresModelPart.GetNode(Id)
				node.SetValue(KratosFemDem.INACTIVE_NODE, True)
				DEMnode.SetValue(KratosFemDem.INACTIVE_NODE, True)
				DEMnode.Set(KratosMultiphysics.TO_ERASE, True)

			# Reset the value to the next step
			node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

			# Remove inactive nodes
			self.SpheresModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)


#============================================================================================================================
	def UpdateDEMVariables(self):

		FEM_Nodes = self.FEM_Solution.main_model_part.Nodes
		DEM_Nodes = self.SpheresModelPart.Nodes

		for DEM_Node in DEM_Nodes:  # Loop over DEM nodes

			if (DEM_Node.GetValue(KratosFemDem.INACTIVE_NODE) == False):

				Id = DEM_Node.Id
				Corresponding_FEM_Node = self.FEM_Solution.main_model_part.GetNode(Id)

				Coordinates    = self.GetNodeCoordinates(Corresponding_FEM_Node)
				Velocity_x     = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
				Velocity_y     = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
				Velocity_z     = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
				Displacement_x = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
				Displacement_y = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
				Displacement_z = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

				# Update Coordinates
				DEM_Node.X = Coordinates[0]
				DEM_Node.Y = Coordinates[1]
				DEM_Node.Z = Coordinates[2]

				# Update Displacements
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, Displacement_x)
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, Displacement_y)
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, Displacement_z)

				# Update Velocities
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, Velocity_x)
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, Velocity_y)
				DEM_Node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, Velocity_z)
