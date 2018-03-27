import os
import KratosMultiphysics
import KratosMultiphysics.ALEApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import mesh_moving_analysis

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class MeshMovingTestCase(KratosUnittest.TestCase):

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # To avoid many prints
        echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        if (echo_level == 0):
            KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def runTest(self):
        mesh_moving_analysis.MeshMovingAnalysis(self.project_parameters).Run()
