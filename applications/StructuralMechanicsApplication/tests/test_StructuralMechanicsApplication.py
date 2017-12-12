# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

# Import the tests o test_classes to create the suits
## SMALL 
# CL tests
from constitutive_law_test import TestConstitutiveLaw as TTestConstitutiveLaw
# Simple patch tests
from test_patch_test_small_strain import TestPatchTestSmallStrain as TTestPatchTestSmallStrain
from test_patch_test_large_strain import TestPatchTestLargeStrain as TTestPatchTestLargeStrain
from test_quadratic_elements import TestQuadraticElements as TTestQuadraticElements
from test_patch_test_shells import TestPatchTestShells as TTestPatchTestShells
from test_patch_test_truss import TestTruss3D2N as TTestTruss3D2N
from test_patch_test_cr_beam import TestCrBeam3D2N as TTestCrBeam3D2N
from test_patch_test_cr_beam import TestCrBeam2D2N as TTestCrBeam2D2N
from test_patch_test_shells_stress import TestPatchTestShellsStressRec as TTestPatchTestShellsStressRec
from test_patch_test_shells_orthotropic import TestPatchTestShellsOrthotropic as TTestPatchTestShellsOrthotropic
# Test loading conditions
from test_loading_conditions_point import TestLoadingConditionsPoint as TestLoadingConditionsPoint
from test_loading_conditions_line import TestLoadingConditionsLine as TestLoadingConditionsLine
from test_loading_conditions_surface import TestLoadingConditionsSurface as TestLoadingConditionsSurface
# Basic moving mesh test
from SmallTests import SimpleMeshMovingTest as TSimpleMeshMovingTest
# Dynamic basic tests
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests
# Patch test Small Displacements
from SmallTests import SDTwoDShearQuaPatchTest as TSDTwoDShearQuaPatchTest
from SmallTests import SDTwoDShearTriPatchTest as TSDTwoDShearTriPatchTest
from SmallTests import SDTwoDTensionQuaPatchTest as TSDTwoDTensionQuaPatchTest
from SmallTests import SDTwoDTensionTriPatchTest as TSDTwoDTensionTriPatchTest
from SmallTests import SDThreeDShearHexaPatchTest as TSDThreeDShearHexaPatchTest
from SmallTests import SDThreeDShearTetraPatchTest as TSDThreeDShearTetraPatchTest
from SmallTests import SDThreeDTensionHexaPatchTest as TSDThreeDTensionHexaPatchTest
from SmallTests import SDThreeDTensionTetraPatchTest as TSDThreeDTensionTetraPatchTest
# Patch test Total Lagrangian
from SmallTests import TLTwoDShearQuaPatchTest as TTLTwoDShearQuaPatchTest
from SmallTests import TLTwoDShearTriPatchTest as TTLTwoDShearTriPatchTest
from SmallTests import TLTwoDTensionQuaPatchTest as TTLTwoDTensionQuaPatchTest
from SmallTests import TLTwoDTensionTriPatchTest as TTLTwoDTensionTriPatchTest
from SmallTests import TLThreeDShearHexaPatchTest as TTLThreeDShearHexaPatchTest
from SmallTests import TLThreeDShearTetraPatchTest as TTLThreeDShearTetraPatchTest
from SmallTests import TLThreeDTensionHexaPatchTest as TTLThreeDTensionHexaPatchTest
from SmallTests import TLThreeDTensionTetraPatchTest as TTLThreeDTensionTetraPatchTest
# Patch test Updated Lagrangian
from SmallTests import ULTwoDShearQuaPatchTest as TULTwoDShearQuaPatchTest
from SmallTests import ULTwoDShearTriPatchTest as TULTwoDShearTriPatchTest
from SmallTests import ULTwoDTensionQuaPatchTest as TULTwoDTensionQuaPatchTest
from SmallTests import ULTwoDTensionTriPatchTest as TULTwoDTensionTriPatchTest
from SmallTests import ULThreeDShearHexaPatchTest as TULThreeDShearHexaPatchTest
from SmallTests import ULThreeDShearTetraPatchTest as TULThreeDShearTetraPatchTest
from SmallTests import ULThreeDTensionHexaPatchTest as TULThreeDTensionHexaPatchTest
from SmallTests import ULThreeDTensionTetraPatchTest as TULThreeDTensionTetraPatchTest
# SPRISM tests
from SmallTests import SprismMembranePatchTests as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests as TSprismBendingPatchTests
# Eigenvalues tests
from SmallTests import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests
from SmallTests import Eigen3D3NThinCircleTests as TEigen3D3NThinCircleTests
# Membrane tests
from SmallTests import Fofi4PointTentnoCableTests as TFofi4PointTentnoCableTests
from SmallTests import Fofi4PointTentCableTests as TFofi4PointTentCableTests
from SmallTests import MembraneQ4PointLoadTests as TMembraneQ4PointLoadTests
from SmallTests import MembraneQ4TrussPointLoadTests as TMembraneQ4TrussPointLoadTests
# 2Node Element tests
from SmallTests import Simple3D2NTrussTest as T3D2NTrussTest
from SmallTests import Simple3D2NTrussLinearTest as T3D2NTrussLinearTest
from SmallTests import Simple3D2NTrussDynamicTest as T3D2NTrussDynamicTest
from SmallTests import Simple3D2NBeamCrTest as T3D2NBeamCrTest
from SmallTests import Simple3D2NBeamCrLinearTest as T3D2NBeamCrLinearTest
from SmallTests import Simple3D2NBeamCrDynamicTest as T3D2NBeamCrDynamicTest
from SmallTests import Simple2D2NBeamCrTest as T2D2NBeamCrTest

# Multipoint constraint tests
from test_multipoint_contstraints import TestMultipointConstraints as TTestMultipointConstraints

# Nodal damping test
from test_nodal_damping import NodalDampingTests as TNodalDampingTests
# Spring damper element test
from test_spring_damper_element import SpringDamperElementTests as TSpringDamperElementTests
# Harmonic analysis tests
from test_harmonic_analysis import HarmonicAnalysisTests as THarmonicAnalysisTests

## NIGHTLY TESTS
# Shell tests
from NightlyTests import ShellT3IsotropicLinearStaticStructScordelisLoRoofTests as TShellT3IsotropicLinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellQ4ThickLinearStaticStructScordelisLoRoofTests as TShellQ4ThickLinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellQ4ThinLinearStaticStructScordelisLoRoofTests as TShellQ4ThinLinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellT3ThickLinearStaticStructScordelisLoRoofTests as TShellT3ThickLinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellT3ThinLinearStaticStructScordelisLoRoofTests as TShellT3ThinLinearStaticStructScordelisLoRoofTests

from NightlyTests import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from NightlyTests import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from NightlyTests import ShellQ4ThickOrthotropicLaminateLinearStaticTests as TShellQ4ThickOrthotropicLaminateLinearStaticTests

from NightlyTests import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from NightlyTests import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests
from NightlyTests import ShellT3ThinOrthotropicLaminateLinearStaticTests as TShellT3ThinOrthotropicLaminateLinearStaticTests

from NightlyTests import ShellT3ThickLinearStaticTests as TShellT3ThickLinearStaticTests
from NightlyTests import ShellT3ThickNonLinearStaticTests as TShellT3ThickNonLinearStaticTests
from NightlyTests import ShellT3ThickLinearDynamicTests as TShellT3ThickLinearDynamicTests
from NightlyTests import ShellT3ThickNonLinearDynamicTests as TShellT3ThickNonLinearDynamicTests
from NightlyTests import ShellT3ThickOrthotropicLaminateLinearStaticTests as TShellT3ThickOrthotropicLaminateLinearStaticTests

from NightlyTests import ShellQ4ThinLinearStaticTests as TShellQ4ThinLinearStaticTests
from NightlyTests import ShellQ4ThinNonLinearStaticTests as TShellQ4ThinNonLinearStaticTests
from NightlyTests import ShellQ4ThinLinearDynamicTests as TShellQ4ThinLinearDynamicTests
from NightlyTests import ShellQ4ThinNonLinearDynamicTests as TShellQ4ThinNonLinearDynamicTests
from NightlyTests import ShellQ4ThinOrthotropicLaminateLinearStaticTests as TShellQ4ThinOrthotropicLaminateLinearStaticTests


from NightlyTests import ShellQ4ThickLinearStaticStructScordelisLoRoofOrthotropicTests as TShellQ4ThickLinearStaticStructScordelisLoRoofOrthotropicTests
from NightlyTests import ShellQ4ThinLinearStaticStructScordelisLoRoofOrthotropicTests as TShellQ4ThinLinearStaticStructScordelisLoRoofOrthotropicTests
from NightlyTests import ShellT3ThickLinearStaticStructScordelisLoRoofOrthotropicTests as TShellT3ThickLinearStaticStructScordelisLoRoofOrthotropicTests
from NightlyTests import ShellT3ThinLinearStaticStructScordelisLoRoofOrthotropicTests as TShellT3ThinLinearStaticStructScordelisLoRoofOrthotropicTests

from NightlyTests import ShellQ4ThickLinearStaticStructPinchedCylinderTests as TShellQ4ThickLinearStaticStructPinchedCylinderTests
from NightlyTests import ShellQ4ThinLinearStaticStructPinchedCylinderTests as TShellQ4ThinLinearStaticStructPinchedCylinderTests
from NightlyTests import ShellT3ThickLinearStaticStructPinchedCylinderTests as TShellT3ThickLinearStaticStructPinchedCylinderTests
from NightlyTests import ShellT3ThinLinearStaticStructPinchedCylinderTests as TShellT3ThinLinearStaticStructPinchedCylinderTests

from NightlyTests import ShellQ4ThickLinearStaticStructPinchedHemisphereTests as TShellQ4ThickLinearStaticStructPinchedHemisphereTests
from NightlyTests import ShellQ4ThinLinearStaticStructPinchedHemisphereTests as TShellQ4ThinLinearStaticStructPinchedHemisphereTests
from NightlyTests import ShellT3ThickLinearStaticStructPinchedHemisphereTests as TShellT3ThickLinearStaticStructPinchedHemisphereTests
from NightlyTests import ShellT3ThinLinearStaticStructPinchedHemisphereTests as TShellT3ThinLinearStaticStructPinchedHemisphereTests

from NightlyTests import ShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughTests as TShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughTests
from NightlyTests import ShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughTests as TShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughTests
from NightlyTests import ShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughTests as TShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughTests
from NightlyTests import ShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughTests as TShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughTests

from NightlyTests import ShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests
from NightlyTests import ShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests
from NightlyTests import ShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests
from NightlyTests import ShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests

from NightlyTests import ShellQ4ThickLinearDynamicStructOscillatingPlateTests as TShellQ4ThickLinearDynamicStructOscillatingPlateTests
from NightlyTests import ShellQ4ThinLinearDynamicStructOscillatingPlateTests as TShellQ4ThinLinearDynamicStructOscillatingPlateTests
from NightlyTests import ShellT3ThickLinearDynamicStructOscillatingPlateTests as TShellT3ThickLinearDynamicStructOscillatingPlateTests
from NightlyTests import ShellT3ThinLinearDynamicStructOscillatingPlateTests as TShellT3ThinLinearDynamicStructOscillatingPlateTests


# CL tests
##from NightlyTests import IsotropicDamageSimoJuPSTest    as TIsotropicDamageSimoJuPSTest

## VALIDATION TESTS
# SPRISM tests
#from ValidationTests import SprismPanTests              as TSprismPanTests
from ValidationTests import PendulusTLTest              as TPendulusTLTest
from ValidationTests import PendulusULTest              as TPendulusULTest

from ValidationTests import ShellQ4ThickLinearStaticUnstructScordelisLoRoofTests as TShellQ4ThickLinearStaticUnstructScordelisLoRoofTests
from ValidationTests import ShellQ4ThinLinearStaticUnstructScordelisLoRoofTests as TShellQ4ThinLinearStaticUnstructScordelisLoRoofTests
from ValidationTests import ShellT3ThickLinearStaticUnstructScordelisLoRoofTests as TShellT3ThickLinearStaticUnstructScordelisLoRoofTests
from ValidationTests import ShellT3ThinLinearStaticUnstructScordelisLoRoofTests as TShellT3ThinLinearStaticUnstructScordelisLoRoofTests

from ValidationTests import ShellQ4ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests as TShellQ4ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests
from ValidationTests import ShellQ4ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests as TShellQ4ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests
from ValidationTests import ShellT3ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests as TShellT3ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests
from ValidationTests import ShellT3ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests as TShellT3ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests

from ValidationTests import ShellQ4ThickLinearStaticUnstructPinchedCylinderTests as TShellQ4ThickLinearStaticUnstructPinchedCylinderTests
from ValidationTests import ShellQ4ThinLinearStaticUnstructPinchedCylinderTests as TShellQ4ThinLinearStaticUnstructPinchedCylinderTests
from ValidationTests import ShellT3ThickLinearStaticUnstructPinchedCylinderTests as TShellT3ThickLinearStaticUnstructPinchedCylinderTests
from ValidationTests import ShellT3ThinLinearStaticUnstructPinchedCylinderTests as TShellT3ThinLinearStaticUnstructPinchedCylinderTests

from ValidationTests import ShellQ4ThickLinearStaticUnstructPinchedHemisphereTests as TShellQ4ThickLinearStaticUnstructPinchedHemisphereTests
from ValidationTests import ShellQ4ThinLinearStaticUnstructPinchedHemisphereTests as TShellQ4ThinLinearStaticUnstructPinchedHemisphereTests
from ValidationTests import ShellT3ThickLinearStaticUnstructPinchedHemisphereTests as TShellT3ThickLinearStaticUnstructPinchedHemisphereTests
from ValidationTests import ShellT3ThinLinearStaticUnstructPinchedHemisphereTests as TShellT3ThinLinearStaticUnstructPinchedHemisphereTests

from ValidationTests import ShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests
from ValidationTests import ShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests
from ValidationTests import ShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests
from ValidationTests import ShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests

from ValidationTests import ShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests
from ValidationTests import ShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests
from ValidationTests import ShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests
from ValidationTests import ShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests

from ValidationTests import ShellQ4ThickLinearDynamicUnstructOscillatingPlateTests as TShellQ4ThickLinearDynamicUnstructOscillatingPlateTests
from ValidationTests import ShellQ4ThinLinearDynamicUnstructOscillatingPlateTests as TShellQ4ThinLinearDynamicUnstructOscillatingPlateTests
from ValidationTests import ShellT3ThickLinearDynamicUnstructOscillatingPlateTests as TShellT3ThickLinearDynamicUnstructOscillatingPlateTests
from ValidationTests import ShellT3ThinLinearDynamicUnstructOscillatingPlateTests as TShellT3ThinLinearDynamicUnstructOscillatingPlateTests

from ValidationTests import ShellQ4ThickLinearDynamicUnstructLumpedOscillatingPlateTests as TShellQ4ThickLinearDynamicUnstructLumpedOscillatingPlateTests
from ValidationTests import ShellQ4ThinLinearDynamicUnstructLumpedOscillatingPlateTests as TShellQ4ThinLinearDynamicUnstructLumpedOscillatingPlateTests
from ValidationTests import ShellT3ThickLinearDynamicUnstructLumpedOscillatingPlateTests as TShellT3ThickLinearDynamicUnstructLumpedOscillatingPlateTests
from ValidationTests import ShellT3ThinLinearDynamicUnstructLumpedOscillatingPlateTests as TShellT3ThinLinearDynamicUnstructLumpedOscillatingPlateTests

from ValidationTests import ShellQ4ThickPendulusTests as TShellQ4ThickPendulusTests
from ValidationTests import ShellQ4ThinPendulusTests as TShellQ4ThinPendulusTests
from ValidationTests import ShellT3ThickPendulusTests as TShellT3ThickPendulusTests
from ValidationTests import ShellT3ThinPendulusTests as TShellT3ThinPendulusTests

from ValidationTests import ShellQ4ThickPendulusLumpedTests as TShellQ4ThickPendulusLumpedTests
from ValidationTests import ShellQ4ThinPendulusLumpedTests as TShellQ4ThinPendulusLumpedTests
from ValidationTests import ShellT3ThickPendulusLumpedTests as TShellT3ThickPendulusLumpedTests
from ValidationTests import ShellT3ThinPendulusLumpedTests as TShellT3ThinPendulusLumpedTests


def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    # Simple patch tests
    ## Solids
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestConstitutiveLaw]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestSmallStrain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestLargeStrain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestQuadraticElements]))
    ## Shells
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShells]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsStressRec]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsOrthotropic]))
    ## Trusses
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestTruss3D2N]))
    ## Beams
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam3D2N]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam2D2N]))
    # Test loading conditions
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingConditionsPoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingConditionsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingConditionsSurface]))
    # Basic moving mesh test
    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    # Dynamic basic tests
    smallSuite.addTest(TDynamicBossakTests('test_execution'))
    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
    # Patch test Small Displacements
    smallSuite.addTest(TSDTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Total Lagrangian
    smallSuite.addTest(TTLTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Updated Lagrangian
    smallSuite.addTest(TULTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionTetraPatchTest('test_execution'))
    # SPRISM tests
    smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
    smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
    # Membrane tests
    smallSuite.addTest(TFofi4PointTentnoCableTests('test_execution'))
    smallSuite.addTest(TFofi4PointTentCableTests('test_execution'))
    smallSuite.addTest(TMembraneQ4PointLoadTests('test_execution'))
    smallSuite.addTest(TMembraneQ4TrussPointLoadTests('test_execution'))
    # 2Node Element tests    
    smallSuite.addTest(T3D2NTrussDynamicTest('test_execution'))
    smallSuite.addTest(T3D2NTrussLinearTest('test_execution'))
    smallSuite.addTest(T3D2NTrussTest('test_execution'))
    smallSuite.addTest(T3D2NBeamCrTest('test_execution'))
    smallSuite.addTest(T3D2NBeamCrLinearTest('test_execution'))
    smallSuite.addTest(T3D2NBeamCrDynamicTest('test_execution'))
    smallSuite.addTest(T2D2NBeamCrTest('test_execution'))
    # Nodal damping test
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TNodalDampingTests]))

    if (missing_external_dependencies == False):
        if (hasattr(KratosMultiphysics.ExternalSolversApplication,
                    "FEASTSolver")):
            # Eigenvalues tests
            smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
            smallSuite.addTest(TEigen3D3NThinCircleTests('test_execution'))
            smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
            # Element damping test
            smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TSpringDamperElementTests]))
            # Harmonic analysis test
            smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([THarmonicAnalysisTests]))
        else:
            print(
                "FEASTSolver solver is not included in the compilation of the External Solvers Application"
            )


    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Shell tests
    nightSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    nightSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    nightSuite.addTest(TShellQ4ThickOrthotropicLaminateLinearStaticTests('test_execution'))

    nightSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    nightSuite.addTest(TShellT3ThinOrthotropicLaminateLinearStaticTests('test_execution'))

    nightSuite.addTest(TShellT3ThickLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellT3ThickOrthotropicLaminateLinearStaticTests('test_execution'))

    nightSuite.addTest(TShellQ4ThinLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinOrthotropicLaminateLinearStaticTests('test_execution'))


    nightSuite.addTest(TShellT3IsotropicLinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellQ4ThickLinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellT3ThinLinearStaticStructScordelisLoRoofTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickLinearStaticStructScordelisLoRoofOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearStaticStructScordelisLoRoofOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearStaticStructScordelisLoRoofOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellT3ThinLinearStaticStructScordelisLoRoofOrthotropicTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickLinearStaticStructPinchedCylinderTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearStaticStructPinchedCylinderTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearStaticStructPinchedCylinderTests('test_execution'))
    nightSuite.addTest(TShellT3ThinLinearStaticStructPinchedCylinderTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickLinearStaticStructPinchedHemisphereTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearStaticStructPinchedHemisphereTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearStaticStructPinchedHemisphereTests('test_execution'))
    nightSuite.addTest(TShellT3ThinLinearStaticStructPinchedHemisphereTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughTests('test_execution'))
    nightSuite.addTest(TShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellT3ThinNonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))

    nightSuite.addTest(TShellQ4ThickLinearDynamicStructOscillatingPlateTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearDynamicStructOscillatingPlateTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearDynamicStructOscillatingPlateTests('test_execution'))
    nightSuite.addTest(TShellT3ThinLinearDynamicStructOscillatingPlateTests('test_execution'))
	
    # CL tests
    ##nightSuite.addTest(TIsotropicDamageSimoJuPSTest('test_execution')) # FIXME: Needs get up to date

    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    # SPRISM tests
    ####validationSuite.addTest(TSprismPanTests('test_execution'))
    validationSuite.addTest(TPendulusTLTest('test_execution'))
    validationSuite.addTest(TPendulusULTest('test_execution'))
    validationSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    validationSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearStaticUnstructScordelisLoRoofTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearStaticUnstructScordelisLoRoofTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearStaticUnstructScordelisLoRoofTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearStaticUnstructScordelisLoRoofTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearStaticUnstructScordelisLoRoofOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearStaticUnstructScordelisLoRoofOrthotropicTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearStaticUnstructPinchedCylinderTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearStaticUnstructPinchedCylinderTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearStaticUnstructPinchedCylinderTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearStaticUnstructPinchedCylinderTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearStaticUnstructPinchedHemisphereTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearStaticUnstructPinchedHemisphereTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearStaticUnstructPinchedHemisphereTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearStaticUnstructPinchedHemisphereTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests('test_execution'))
    validationSuite.addTest(TShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughTests('test_execution'))
    validationSuite.addTest(TShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellT3ThickNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellT3ThinNonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearDynamicUnstructOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearDynamicUnstructOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearDynamicUnstructOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearDynamicUnstructOscillatingPlateTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickLinearDynamicUnstructLumpedOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearDynamicUnstructLumpedOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearDynamicUnstructLumpedOscillatingPlateTests('test_execution'))
    validationSuite.addTest(TShellT3ThinLinearDynamicUnstructLumpedOscillatingPlateTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickPendulusTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinPendulusTests('test_execution'))
    validationSuite.addTest(TShellT3ThickPendulusTests('test_execution'))
    validationSuite.addTest(TShellT3ThinPendulusTests('test_execution'))

    validationSuite.addTest(TShellQ4ThickPendulusLumpedTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinPendulusLumpedTests('test_execution'))
    validationSuite.addTest(TShellT3ThickPendulusLumpedTests('test_execution'))
    validationSuite.addTest(TShellT3ThinPendulusLumpedTests('test_execution'))
    
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
