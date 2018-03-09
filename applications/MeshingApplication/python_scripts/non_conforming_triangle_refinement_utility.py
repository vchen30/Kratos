from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

import pdb

class NonConformingTriangleRefinementUtility:
    def __init__(self, origin_model_part, destination_model_part):
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part
    
    def Initialize(self):
        pass
    
    def Execute(self):
        # pdb.set_trace()
        # Copy the nodes to the new model_part
        self.destination_model_part.Nodes = self.origin_model_part.Nodes
        # Initialize the hash to ensure that each node is created once
        self.nodes_hash = {}
        # Find the last node Id
        last_id = [6]

        # Loop the elements to add nodes and create the new nodes
        for element in self.origin_model_part.Elements:
            # Initialize the list of middle nodes id
            mid_nodes = [0, 1, 2]

            # First edge
            mid_nodes[0] = self.CreateMiddleNode(element.GetNode(0), element.GetNode(1), last_id)
            # Second edge
            mid_nodes[1] = self.CreateMiddleNode(element.GetNode(1), element.GetNode(2), last_id)
            # Third edge
            mid_nodes[2] = self.CreateMiddleNode(element.GetNode(2), element.GetNode(0), last_id)

    def CreateMiddleNode(self, node_a, node_b, last_id):
        new_node_hash = (min(node_a.Id, node_b.Id), max(node_a.Id, node_b.Id))
        # print('Do I need to create a new node?')
        # print('This is the nodes hash:')
        # print(self.nodes_hash)
        # print('And this would be the new node hash:')
        # print(new_node_hash)
        if not new_node_hash in self.nodes_hash:
            # print('The new node is not in the hash, so create a new node')
            new_x = 0.5 * ( node_a.X + node_b.X )
            new_y = 0.5 * ( node_a.Y + node_b.Y )
            new_id = last_id[0] + 1
            self.destination_model_part.CreateNewNode(new_id, new_x, new_y, 0)  # Añadir el nodo en el punto medio
            self.nodes_hash[new_node_hash] = new_id       # Tiene que ser el Id del nuevo nodo
            last_id[0] = new_id
            # print('The new node Id is : ', new_id[0])
            # print('')
        else:
            new_id = self.nodes_hash[new_node_hash]

        return new_id