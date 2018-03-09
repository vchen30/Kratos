from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

class NonConformingTriangleRefinementUtility:
    def __init__(self, origin_model_part, destination_model_part):
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part
    
    def Initialize(self):
        pass
    
    def Execute(self):
        self.destination_model_part.Nodes = self.origin_model_part.Nodes
        nodes_hash = {}
        for element in self.origin_model_part.Elements:
            # print(element.__dir__())
            # for i in range(3):
            #     # print(nodes.__dir__())
            #     new_node_hash = (element.GetNode(i).Id, element.GetNode(i+1).Id)
            #     print(new_node_hash)
            #     if not new_node_hash in nodes_hash:
            #         self.destination_model_part.AddNode
            #         nodes_hash[new_node_hash] = 1
            
            # First edge
            new_node_hash = (element.GetNode(0).Id, element.GetNode(1).Id)
            if not new_node_hash in nodes_hash:
                # self.destination_model_part.AddNode()  # Añadir el nodo en el punto medio
                # nodes_hash[new_node_hash] = 1       # Tiene que ser el Id del nuevo nodo