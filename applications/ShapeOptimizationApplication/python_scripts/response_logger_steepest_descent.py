# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Import logger base classes
from response_logger_base import ResponseLogger

# Import additional libraries
import csv
from custom_timer import Timer

# ==============================================================================
class ResponseLoggerSteepestDescent( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimizationSettings ):
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings

        self.onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()

        self.completeResponseLogFileName = self.__CreateCompleteResponseLogFilename( optimizationSettings )

        self.objectiveHistory = {}
        self.objectiveOutputReference = None
        self.absoluteChangeOfObjectiveHistory = {}
        self.relativeChangeOfObjectiveHistory = {}

        self.currentIteration = 0
        self.previousIteration = 0
        self.initialIteration = 0

    # --------------------------------------------------------------------------
    def __CreateCompleteResponseLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        responseLogFilename = optimizationSettings["output"]["response_log_filename"].GetString()+".csv"
        return os.path.join( resultsDirectory, responseLogFilename )

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.completeResponseLogFileName, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:<4s}".format("itr"))
            row.append("{:>20s}".format("f"))
            row.append("{:>12s}".format("df_abs[%]"))
            row.append("{:>12s}".format("df_rel[%]"))
            row.append("{:>13s}".format("step_size[-]"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def LogCurrentResponses( self, optimizationIteration ):
        self.currentIteration = optimizationIteration
        if self.__IsFirstLog():
            self.initialIteration = optimizationIteration
            self.__AddObjectiveValueToHistory()
            self.__DetermineReferenceValuesForOutput()
            self.__InitializeChangeOfObjectiveHistory()
        else:
            self.__AddObjectiveValueToHistory()
            self.__AddChangeOfObjectiveToHistory()
        self.__PrintInfoAboutResponseFunctionValues()
        self.__WriteDataToLogFile()
        self.previousIteration = optimizationIteration

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass # No finalization necessary here

    # --------------------------------------------------------------------------
    def GetValue( self, variableKey ):
        if variableKey=="RELATIVE_CHANGE_OF_OBJECTIVE_VALUE":
            if self.__IsFirstLog():
                raise RuntimeError("Relative change of objective function can not be computed since only one logged value is existing!")
            else:
                return self.relativeChangeOfObjectiveHistory[self.currentIteration]
        else:
            raise NameError("Value with the following variable key not defined in response_logger_penalized_projection.py: " + variableKey)

    # -------------------------------------------------------------------------
    def __IsFirstLog( self ):
        if len(self.objectiveHistory) == 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __AddObjectiveValueToHistory( self ):
        objectiveValue = self.communicator.getValue ( self.onlyObjective )
        self.objectiveHistory[self.currentIteration] = objectiveValue

    # --------------------------------------------------------------------------
    def __DetermineReferenceValuesForOutput( self ):
        self.objectiveOutputReference = self.objectiveHistory[self.initialIteration]

        if abs(self.objectiveOutputReference)<1e-12:
            print("\n> WARNING: Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
            self.objectiveOutputReference = 1.0
        else:
            self.objectiveOutputReference = self.objectiveHistory[self.initialIteration]

    # --------------------------------------------------------------------------
    def __InitializeChangeOfObjectiveHistory( self ):
        self.absoluteChangeOfObjectiveHistory[self.currentIteration] = 0.0
        self.relativeChangeOfObjectiveHistory[self.currentIteration] = 0.0

    # --------------------------------------------------------------------------
    def __AddChangeOfObjectiveToHistory( self ):
        objectiveValue = self.objectiveHistory[self.currentIteration]
        previousObjectiveValue = self.objectiveHistory[self.previousIteration]

        self.absoluteChangeOfObjectiveHistory[self.currentIteration] = 100*(objectiveValue-self.objectiveOutputReference) / abs(self.objectiveOutputReference)
        self.relativeChangeOfObjectiveHistory[self.currentIteration] = 100*(objectiveValue-previousObjectiveValue) / abs(self.objectiveOutputReference)

    # --------------------------------------------------------------------------
    def __PrintInfoAboutResponseFunctionValues( self ):
        objectiveValue = self.objectiveHistory[self.currentIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveHistory[self.currentIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveHistory[self.currentIteration]

        print("\n> Current value of objective function = ", round(objectiveValue,12))
        print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,4)," [%]")
        print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,4)," [%]")

    # --------------------------------------------------------------------------
    def __WriteDataToLogFile( self ):
        objectiveValue = self.objectiveHistory[self.currentIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveHistory[self.currentIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveHistory[self.currentIteration]

        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:<4s}".format(str(self.currentIteration)))
            row.append(str("{:>20f}".format(objectiveValue)))
            row.append(str("{:>12f}".format(absoluteChangeOfObjectiveValue)))
            row.append(str("{:>12f}".format(relativeChangeOfObjectiveValue)))
            row.append(str("{:>13f}".format(self.optimizationSettings["line_search"]["step_size"].GetDouble())))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
