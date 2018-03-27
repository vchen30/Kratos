# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# ==============================================================================
def CreateCommunicator( optimization_settings ):
    return Communicator( optimization_settings )

# ==============================================================================
class Communicator:

    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings ):
        self.optimization_settings = optimization_settings
        self.supported_objective_types = ["minimization", "maximization"]
        self.supported_constraint_types = ["=", "<", ">", "<=", ">="]
        self.supported_constraint_references = ["initial_value", "specified_value"]

        self.__initializeListOfRequests()
        self.__initializeListOfResponses()

    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def requestValueOf( self, response_id ):
        self.list_of_requests[response_id]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, response_id ):
        self.list_of_requests[response_id]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingValueOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportValue( self, response_id, value ):
        self.__storeValue( response_id, value)
        if self.__isResponseWaitingForInitialValueAsReference( response_id ):
            self.__setValueAsReference( response_id, value )

        standardized_value = self.__translateValueToStandardForm( response_id, value )
        self.__storeStandardizedValue( response_id, standardized_value )

    # --------------------------------------------------------------------------
    def reportGradient( self, response_id, gradient ):
        standardized_gradient = self.__translateGradientToStandardForm( response_id, gradient )
        self.__storeStandardizedGradient( response_id, standardized_gradient )

    # --------------------------------------------------------------------------
    def getValue( self, response_id ):
        if self.list_of_responses[response_id]["value"] is None:
            raise RuntimeError("The requested value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["value"]

    # --------------------------------------------------------------------------
    def getReferenceValue( self, response_id ):
        if self.list_of_responses[response_id]["reference_value"] is None:
            raise RuntimeError("The requested reference_value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["reference_value"]

    # --------------------------------------------------------------------------
    def getStandardizedValue( self, response_id ):
        if self.list_of_responses[response_id]["standardized_value"] is None:
            raise RuntimeError("The requested standardized_value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["standardized_value"]

    # --------------------------------------------------------------------------
    def getStandardizedGradient( self, response_id ):
        if self.list_of_responses[response_id]["standardized_gradient"] is None:
            raise RuntimeError("The requested standardized_gradient for ", response_id, " is None!")
        return self.list_of_responses[response_id]["standardized_gradient"]

    # --------------------------------------------------------------------------
    def __initializeListOfRequests( self ):
        self.list_of_requests = {}

        for objective_number in range( self.optimization_settings["objectives"].size()):
            objective_id = self.optimization_settings["objectives"][objective_number]["identifier"].GetString()
            self.list_of_requests[objective_id] = {"calculateValue": False, "calculateGradient": False}

        for constraint_number in range(self.optimization_settings["constraints"].size()):
            constraint_id = self.optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            self.list_of_requests[constraint_id] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def __initializeListOfResponses( self ):
        self.list_of_responses = {}
        self.__addObjectivesToListOfResponses()
        self.__addConstraintsToListOfResponses()

    # --------------------------------------------------------------------------
    def __addObjectivesToListOfResponses( self ):
        for objective_number in range(self.optimization_settings["objectives"].size()):
            objective = self.optimization_settings["objectives"][objective_number]
            objective_id =  objective["identifier"].GetString()

            if objective["type"].GetString() not in self.supported_objective_types:
                raise RuntimeError("Unsupported type defined for the following objective: " + objective_id)

            self.list_of_responses[objective_id] = { "type"                 : objective["type"].GetString(),
                                                     "value"                : None,
                                                     "standardized_value"   : None,
                                                     "standardized_gradient": None }

    # --------------------------------------------------------------------------
    def __addConstraintsToListOfResponses( self ):
        for constraint_number in range( self.optimization_settings["constraints"].size()):
            constraint = self.optimization_settings["constraints"][constraint_number]
            constraint_id = self.optimization_settings["constraints"][constraint_number]["identifier"].GetString()

            if constraint["type"].GetString() not in self.supported_constraint_types:
                raise RuntimeError("Unsupported type defined for the following constraint: " + constraint_id )

            if  constraint["reference"].GetString() == "specified_value":
                self.list_of_responses[constraint_id] = { "type"                 : constraint["type"].GetString(),
                                                          "value"                : None,
                                                          "standardized_value"   : None,
                                                          "standardized_gradient": None,
                                                          "reference_value"      : constraint["reference_value"].GetDouble() }
            elif constraint["reference"].GetString() == "initial_value":
                self.list_of_responses[constraint_id] = { "type"                 : constraint["type"].GetString(),
                                                          "value"                : None,
                                                          "standardized_value"   : None,
                                                          "standardized_gradient": None,
                                                          "reference_value"      : "waiting_for_initial_value" }
            else:
                raise RuntimeError("Unsupported reference defined for the following constraint: " + constraint_id )

    # --------------------------------------------------------------------------
    def __deleteAllRequests( self ):
        for response_id in self.list_of_requests:
            self.list_of_requests[response_id]["calculateValue"] = False
            self.list_of_requests[response_id]["calculateGradient"] = False

    # --------------------------------------------------------------------------
    def __deleteAllReportedValues( self ):
        for response_id in self.list_of_responses:
            self.list_of_responses[response_id]["value"] = None
            self.list_of_responses[response_id]["standardized_value"] = None
            self.list_of_responses[response_id]["standardized_gradient"] = None

    # --------------------------------------------------------------------------
    def __isResponseWaitingForInitialValueAsReference( self, response_id ):
        response = self.list_of_responses[response_id]
        if "reference_value" in response and response["reference_value"] == "waiting_for_initial_value":
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __setValueAsReference( self, response_id, value ):
        self.list_of_responses[response_id]["reference_value"] = value

    # --------------------------------------------------------------------------
    def __translateValueToStandardForm( self, response_id, value ):
        response = self.list_of_responses[response_id]
        response_type = response["type"]
        if response_type in self.supported_objective_types:
            if response_type == "maximization":
                return -value
            else:
                return value
        else:
            if response_type == ">" or response_type == ">=":
                return -(value-response["reference_value"])
            else:
                return (value-response["reference_value"])

    # --------------------------------------------------------------------------
    def __translateGradientToStandardForm( self, response_id, gradient ):
        response = self.list_of_responses[response_id]
        response_type = response["type"]
        if response_type == "maximization" or response_type == ">" or response_type == ">=":
            for local_gradient in gradient.values():
                local_gradient[0] = -local_gradient[0]
                local_gradient[1] = -local_gradient[1]
                local_gradient[2] = -local_gradient[2]
            return gradient
        else:
            return gradient

    # --------------------------------------------------------------------------
    def __storeValue( self, response_id, value ):
        self.list_of_responses[response_id]["value"] = value

    # --------------------------------------------------------------------------
    def __storeStandardizedValue( self, response_id, value ):
        self.list_of_responses[response_id]["standardized_value"] = value

    # --------------------------------------------------------------------------
    def __storeStandardizedGradient( self, response_id, gradient ):
        self.list_of_responses[response_id]["standardized_gradient"] = gradient

# ==============================================================================
