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
        
        // slave side: mark as in contact if a node has a Lagrange Multiplier unequal to zero
        ModelPart::NodesContainerType& rNodes = mThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            if(i_nodes->Has(AUGMENTED_NORMAL_CONTACT_PRESSURE)== true){
                //std::cout<<"node "<<i_nodes->Id()<<" has contact pressure "<<i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE)<<std::endl;
                if(i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE)!= 0)
                    i_nodes->SetValue(CONTACT_PRESSURE,i_nodes->GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE));
            }
        }  
        
        ModelPart::ConditionsContainerType conditions = mThisModelPart.Conditions();
        for(auto i_condition = conditions.begin(); i_condition!= conditions.end(); i_condition++)
        {
            // Slave side: Use ALM directly
            //for(auto j=0; j<i_condition->GetGeometry().size();j++)
            //    i_condition->GetGeometry()[j].SetValue(CONTACT_PRESSURE,i_condition->GetGeometry()[j].GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE));
            
            // Master side: Project ALM from the slave side to the master side
            boost::shared_ptr<ConditionMap>& all_conditions_maps = i_condition->GetValue( MAPPING_PAIRS );

            GeometryType::CoordinatesArrayType projected_1, projected_2;
            for (auto i_master = all_conditions_maps->begin(); i_master != all_conditions_maps->end(); i_master++){
                //if (i_master->second){
                if (true){
                    std::cout<<"masterelement active:"<<i_master->second<<std::endl;
                    for(int i=0; i < i_master->first->GetGeometry().size(); i++){
                        
                        array_1d<double,3> i_nodes= i_master->first->GetGeometry()[i].Coordinates();
                        PointType master_point = PointType(i_nodes);
                        PointType slave_point = PointType(0);
                        array_1d<double,3> local_coord = array_1d<double,3>(0);
                        // Projection of the master node to the slave node
                        //std::cout<<"Normal Master: "<<i_master->first->GetGeometry()[i].GetValue(NORMAL)<<std::endl;
                        MortarUtilities::FastProjectDirection(i_condition->GetGeometry(),
                                                            master_point,
                                                            slave_point,
                                                            i_master->first->GetGeometry()[i].GetValue(NORMAL),
                                                            i_condition->GetGeometry().Normal(local_coord));
                        //std::cout<<"point to project: "<<master_point<<" , projected point: "<<slave_point;

                        //get the local coordinates of the projected point
                        GeometryType::CoordinatesArrayType projected_point_local;
                        i_condition->GetGeometry().PointLocalCoordinates(projected_point_local,slave_point.Coordinates());
                        //std::cout<<" local coordinates: "<<projected_point_local<<std::endl;

                        

                        // if the projected point is inside the element: set contact pressure
                        if(projected_point_local[0]<=1.0 && projected_point_local[0]>= -1.0){
                            Vector shape_function_value;
                            i_condition->GetGeometry().ShapeFunctionsValues(shape_function_value,projected_point_local);
                            double contact_pressure=0;
                            for (auto j=0;j< i_condition->GetGeometry().size();j++)
                                contact_pressure += i_condition->GetGeometry()[j].GetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE)
                                                                                    *shape_function_value[j];
                            i_master->first->GetGeometry()[i].SetValue(CONTACT_PRESSURE, contact_pressure);
                            std::cout<<"node "<<i_master->first->GetGeometry()[i].Id()<<" has contact pressure "<<i_master->first->GetGeometry()[i].GetValue(CONTACT_PRESSURE)<<std::endl;
                        }

                    
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


