//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/compute_colors_utility.h"

namespace Kratos {

void ComputeColorsUtility::Initialize()
{        
    // Initialize and create the auxiliar maps
    const std::vector<std::string> sub_model_part_names = mrThisModelPart.GetSubModelPartNames();
    std::unordered_map<int,std::set<int>> aux_nodes_colors, aux_cond_colors, aux_elem_colors;
    
    std::vector<std::string> model_part_names;
    model_part_names.push_back(mrThisModelPart.Name());
    for (const auto & sub_model_part_name : sub_model_part_names)
        model_part_names.push_back(sub_model_part_name);
    
    // Initialize Colors
    int color = 0;
    for (SizeType i_sub_model_part = 0; i_sub_model_part < model_part_names.size(); ++i_sub_model_part) {
        mColors[i_sub_model_part].push_back(model_part_names[i_sub_model_part]);
        
        if (color > 0) {            
            ModelPart& r_sub_model_part = mrThisModelPart.GetSubModelPart(model_part_names[i_sub_model_part]);

            /* Nodes */
            NodesArrayType& nodes_array = r_sub_model_part.Nodes();
            for(SizeType i = 0; i < nodes_array.size(); ++i) 
                aux_nodes_colors[(nodes_array.begin() + i)->Id()].insert(color);
            
            /* Conditions */
            ConditionsArrayType& conditions_array = r_sub_model_part.Conditions();
            for(SizeType i = 0; i < conditions_array.size(); ++i) 
                aux_cond_colors[(conditions_array.begin() + i)->Id()].insert(color);
            
            /* Elements */
            ElementsArrayType& elements_array = r_sub_model_part.Elements();
            for(SizeType i = 0; i < elements_array.size(); ++i) 
                aux_elem_colors[(elements_array.begin() + i)->Id()].insert(color);
        }
        
        color += 1;
    }
    
    // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously 
    std::unordered_map<std::set<int>, int, KeyHasherRange<std::set<int>>, KeyComparorRange<std::set<int>> > combinations;
    
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const std::set<int>& value = aux_nodes_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }
    
    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const std::set<int>& value = aux_cond_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const std::set<int>& value = aux_elem_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }
    
    /* Combinations */
    for(auto & combination : combinations) {
        const std::set<int>& key = combination.first;
        for(int it : key) 
            mColors[color].push_back(mColors[it][0]);
        combinations[key] = color;
        color += 1;
    }
    
    // The final maps are created
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const int key = aux_nodes_color.first;
        const std::set<int>& value = aux_nodes_color.second;
        
        if (value.size() == 0)
            mNodesColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            mNodesColors[key] = *value.begin();
        else // There is a combination
            mNodesColors[key] = combinations[value];
    }
    
    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const int key = aux_cond_color.first;
        const std::set<int>& value = aux_cond_color.second;
        
        if (value.size() == 0)
            mConditionsColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            mConditionsColors[key] = *value.begin();
        else // There is a combination
            mConditionsColors[key] = combinations[value];
    }
    
    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const int key = aux_elem_color.first;
        const std::set<int>& value = aux_elem_color.second;
        
        if (value.size() == 0)
            mElementsColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            mElementsColors[key] = *value.begin();
        else // There is a combination
            mElementsColors[key] = combinations[value];
    }
}
    
/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos.
