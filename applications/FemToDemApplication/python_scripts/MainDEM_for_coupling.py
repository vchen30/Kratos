from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import main_script as MainDEM



class DEM_Solution(MainDEM.Solution):

    def Info(self):
        print("DEM part of the FEM-DEM application")