//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Anna Rehr $
//   Date:                $Date: 2017-10-13 08:56:42 $
//   email:               anna.rehr@gmx.de
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

// #if !defined(EXTRACT_NODAL_CONTACT_PRESSURE )
// #define  EXTRACT_NODAL_CONTACT_PRESSURE

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/condition.h"
#include "includes/mortar_classes.h"
#include "utilities/mortar_utilities.h"

namespace Kratos
{

/*
This class extracts the nodal contact pressure for the nodes in contact.
For the slave body, the Lagrange Multiplier can be directly extracted,
for the master body a projection has to be executed wich projects the slave body LM to the master side.
*/

class ExtractNodalContactPressure
{
public:

    ExtractNodalContactPressure(ModelPart& mp): mThisModelPart (mp)
    {
    }

    void Execute()
    {
        std::cout<<"set bool for nodal contact"<<std::endl;

        // slave side: mark as in contact if a node has a Lagrange Multiplier unequal to zero
        ModelPart::NodesContainerType& rNodes = mThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            if(i_nodes->Has(AUGMENTED_NORMAL_CONTACT_PRESSURE)== true){
                std::cout<<"node "<<i_nodes->Id()<<" has contact pressure "<<i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE)<<std::endl;
                if(i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE)< 0)
                    i_nodes->SetValue(CONTACT_PRESSURE,i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE));
            }
        }   
        
        ModelPart::ConditionsContainerType conditions = mThisModelPart.Conditions();
        for(auto i_condition = conditions.begin(); i_condition!= conditions.end(); i_condition++)
        {
            unsigned int pair_number;
            boost::shared_ptr<ConditionMap>& all_conditions_maps = i_condition->GetValue( MAPPING_PAIRS );
            pair_number = all_conditions_maps->size();
            std::cout<<"condition"<<i_condition->Id()<<" has "<<pair_number<<" master elements."<<std::endl;
            
            for (auto i_master = all_conditions_maps->begin(); i_master != all_conditions_maps->end(); i_master++){
                if (i_master->second){
                //if (true){
                    std::cout<<"masterelement active:"<<i_master->second<<std::endl;
                    for(int i=0; i < i_master->first->GetGeometry().size(); i++){
                        
                        array_1d<double,3> i_nodes= i_master->first->GetGeometry()[i].Coordinates();
                        PointType master_point = PointType(i_nodes);
                        PointType slave_point = PointType(0);
                        array_1d<double,3> local_coord = array_1d<double,3>(0);
                        //array_1d<double, 3> normal_slave_element = Geometry::Normal()
                        // Projection of the master node to the slave node
                        std::cout<<"Normal Master: "<<i_master->first->GetGeometry()[i].GetValue(NORMAL)<<std::endl;
                        MortarUtilities::FastProjectDirection(i_condition->GetGeometry(),
                                                            master_point,
                                                            slave_point,
                                                            i_master->first->GetGeometry()[i].GetValue(NORMAL),
                                                            i_condition->GetGeometry().Normal(local_coord));
                    
                    }
                }
            }

        }
    }

private:
    ModelPart mThisModelPart;
};



}  // namespace Kratos.

// #endif // EXTRACT_NODAL_CONTACT_PRESSURE  defined 


