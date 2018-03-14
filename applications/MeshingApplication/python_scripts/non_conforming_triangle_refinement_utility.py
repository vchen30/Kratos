from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

class NonConformingTriangleRefinementUtility:
    def __init__(self, model_part):
        self.model_part = model_part

    def Initialize(self):
        pass

    def Execute(self):
        # Initialize the hash to ensure that each node is created once
        self.nodes_hash = {}
        # Find the last node Id
        last_node_id = [self.GetMaxNodeId()]
        # Find the last condition Id
        last_cond_id = self.GetMaxConditionId()
        # Find the last element Id
        elem_id_list, last_elem_id = self.MarkInputElementsAndGetId()

        # Loop the origin elements to add nodes and create the new elements
        for this_id in elem_id_list:
            # Get the element
            element = self.model_part.GetElement(this_id)

            # Initialize the id list of the middle nodes
            mid_nodes = [0, 1, 2]

            # First edge
            mid_nodes[0] = self.CreateMiddleNode(element.GetNode(0), element.GetNode(1), last_node_id)
            # Second edge
            mid_nodes[1] = self.CreateMiddleNode(element.GetNode(1), element.GetNode(2), last_node_id)
            # Third edge
            mid_nodes[2] = self.CreateMiddleNode(element.GetNode(2), element.GetNode(0), last_node_id)

            # And now, create the new sub elements
            ## First sub element
            last_elem_id += 1
            self.model_part.CreateNewElement("Element2D3N", 
                last_elem_id,
                [element.GetNode(0).Id, mid_nodes[0], mid_nodes[2]],
                self.model_part.GetProperties()[0] )
            
            ## Second sub element
            last_elem_id += 1
            self.model_part.CreateNewElement("Element2D3N", 
                last_elem_id,
                [element.GetNode(1).Id, mid_nodes[1], mid_nodes[0]],
                self.model_part.GetProperties()[0] )
            
            ## Third sub element
            last_elem_id += 1
            self.model_part.CreateNewElement("Element2D3N", 
                last_elem_id,
                [element.GetNode(2).Id, mid_nodes[2], mid_nodes[1]],
                self.model_part.GetProperties()[0] )
            
            ## Fourth sub element
            last_elem_id += 1
            self.model_part.CreateNewElement("Element2D3N", 
                last_elem_id,
                [mid_nodes[0], mid_nodes[1], mid_nodes[2]],
                self.model_part.GetProperties()[0] )
            
            # Now we can delete the element
            self.model_part.RemoveElement(this_id)

        for condition in self.model_part.Conditions:
            # Create nodes inside the condition
            mid_node = self.CreateMiddleNode(condition.GetNode(0), condition.GetNode(1), last_node_id)

            # And now, create the sub conditions
            ## First sub condition
            last_cond_id += 1
            self.model_part.CreateNewCondition("Condition2D2N",
                last_cond_id,
                [condition.GetNode(0).Id, mid_node],
                self.model_part.GetProperties()[0] )

            ## Second sub condition
            last_cond_id += 1
            self.model_part.CreateNewCondition("Condition2D2N",
                last_cond_id,
                [mid_node, condition.GetNode(1).Id],
                self.model_part.GetProperties()[0] )
            
            # Now we can delete the condition
            self.model_part.RemoveCondition(condition.Id)


    def MarkInputElementsAndGetId(self):
        id_list = []
        max_id = 0
        for element in self.model_part.Elements:
            element.Set(KratosMultiphysics.TO_ERASE, True)
            id_list.append(element.Id)
            if element.Id > max_id:
                max_id = element.Id
        return id_list, max_id

    def GetMaxNodeId(self):
        max_id = 0
        for node in self.model_part.Nodes:
            if node.Id > max_id:
                max_id = node.Id
        return max_id
    
    def GetMaxConditionId(self):
        max_id = 0
        for condition in self.model_part.Conditions:
            if condition.Id > max_id:
                max_id = condition.Id
        return max_id

    def CreateMiddleNode(self, node_a, node_b, last_id):
        new_node_hash = (min(node_a.Id, node_b.Id), max(node_a.Id, node_b.Id))
        if not new_node_hash in self.nodes_hash:
            new_x = 0.5 * ( node_a.X + node_b.X )
            new_y = 0.5 * ( node_a.Y + node_b.Y )
            new_id = last_id[0] + 1
            self.model_part.CreateNewNode(new_id, new_x, new_y, 0)  # Añadir el nodo en el punto medio
            self.nodes_hash[new_node_hash] = new_id       # Tiene que ser el Id del nuevo nodo
            last_id[0] = new_id
        else:
            new_id = self.nodes_hash[new_node_hash]

        return new_id
