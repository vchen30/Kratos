import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Structural_Test as Execute_Test

# This utiltiy will control the execution scope in case we need to acces files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class StructuralMechanichsTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Initialize GiD  I/O
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Creating the model part
            self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass
    
class SprismPanTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/pan_test"
    
class PendulusTLTest(StructuralMechanichsTestFactory):
    file_name = "pendulus_test/pendulus_TL_test"
    
class PendulusULTest(StructuralMechanichsTestFactory):
    file_name = "pendulus_test/pendulus_UL_test"

## From here new shell tests
class ShellQ4ThickLinearStaticUnstructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_static_unstruct_scordelis_lo_roof"
class ShellQ4ThinLinearStaticUnstructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_static_unstruct_scordelis_lo_roof"
class ShellT3ThickLinearStaticUnstructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_static_unstruct_scordelis_lo_roof"
class ShellT3ThinLinearStaticUnstructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_static_unstruct_scordelis_lo_roof"

class ShellQ4ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_static_unstruct_scordelis_lo_roof_orthotropic"
class ShellQ4ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_static_unstruct_scordelis_lo_roof_orthotropic"
class ShellT3ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_static_unstruct_scordelis_lo_roof_orthotropic"
class ShellT3ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_static_unstruct_scordelis_lo_roof_orthotropic"

class ShellQ4ThickLinearStaticUnstructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_static_unstruct_pinched_cylinder"
class ShellQ4ThinLinearStaticUnstructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_static_unstruct_pinched_cylinder"
class ShellT3ThickLinearStaticUnstructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_static_unstruct_pinched_cylinder"
class ShellT3ThinLinearStaticUnstructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_static_unstruct_pinched_cylinder"

class ShellQ4ThickLinearStaticUnstructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_static_unstruct_pinched_hemisphere"
class ShellQ4ThinLinearStaticUnstructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_static_unstruct_pinched_hemisphere"
class ShellT3ThickLinearStaticUnstructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_static_unstruct_pinched_hemisphere"
class ShellT3ThinLinearStaticUnstructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_static_unstruct_pinched_hemisphere"

class ShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"
class ShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"
class ShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"
class ShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"

class ShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"
class ShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"
class ShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"
class ShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellQ4ThickLinearDynamicUnstructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellQ4ThinLinearDynamicUnstructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThickLinearDynamicUnstructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThinLinearDynamicUnstructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_dynamic_unstruct_oscillating_plate"

class ShellQ4ThickLinearDynamicUnstructLumpedOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellQ4ThinLinearDynamicUnstructLumpedOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThickLinearDynamicUnstructLumpedOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThinLinearDynamicUnstructLumpedOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_dynamic_unstruct_oscillating_plate"

class ShellQ4ThickPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellQ4ThinPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThickPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThinPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_dynamic_unstruct_oscillating_plate"

class ShellQ4ThickPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellQ4ThinPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_thin_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThickPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thick_linear_dynamic_unstruct_oscillating_plate"
class ShellT3ThinPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_thin_linear_dynamic_unstruct_oscillating_plate"
