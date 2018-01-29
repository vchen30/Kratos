//
//   Project Name:        KratosFemToDemApplication $
//   Last modified by:    $Author:Alejandro Cornejo $
//   Date:                $Date:        Oct 2017    $
//   Revision:            $Revision:            0.0 $
//

#if !defined(KRATOS_STRESS_TO_NODES_PROCESS )
#define  KRATOS_STRESS_TO_NODES_PROCESS

#include <fstream>
#include <cmath>

#include "includes/model_part.h"
#include "boost/smart_ptr.hpp"
#include "processes/process.h"

#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/define.h"

namespace Kratos
{

    // ONLY FOR TRIANGLES IN 2D

    class StressToNodesProcess : public Process
    {

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    protected:
        
        struct NodeStresses
        {
            Vector EffectiveStressVector;
            int NElems;
            
            NodeStresses()
            {
                EffectiveStressVector = ZeroVector(3); // 2D: sigma=[Sxx,Syy,Sxy]
                NElems = 0;
            }
        };

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    public:
        
        typedef ModelPart::ElementsContainerType ElementsArrayType;

        // Constructor
        StressToNodesProcess(ModelPart& r_model_part) : mr_model_part(r_model_part)
        {
            mNNodes = mr_model_part.NumberOfNodes();
            mNElements = mr_model_part.NumberOfElements();
        }

        // Destructor
        virtual ~StressToNodesProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void Execute()
        {

            NodeStresses* NodeStressesVector = new NodeStresses[mNNodes];

            this->StressExtrapolationAndSmoothing(NodeStressesVector);

			//this->WriteAuxPost(NodeStressesVector);

            delete[] NodeStressesVector;
        }

        // --------------------------------------------------------------------
        void StressExtrapolationAndSmoothing(NodeStresses* pNodeStressesVector)
        {

            Vector GaussPointsStresses = ZeroVector(3);

            // Loop over elements to extrapolate the stress to the nodes
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
    
                bool condition_is_active = true;
                if ((*it)->IsDefined(ACTIVE))
                {
                    condition_is_active = (*it)->Is(ACTIVE);
                }
    
                if (condition_is_active)
                {
                   GaussPointsStresses = (*it)->GetValue(STRESS_VECTOR);
    
                    //Triangles2D3N
                    if((*it)->GetGeometry().PointsNumber() == 3)
                    {
                        for(int i = 0; i < 3; i++)
                        {   
                            pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[0] += GaussPointsStresses[0];
                            pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[1] += GaussPointsStresses[1];
                            pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[2] += GaussPointsStresses[2];
                            pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].NElems += 1;
                        }
                    }
                }
            }
    
            // Ponderate over the elements coincident on that node
            for(unsigned int i = 0; i < mNNodes; i++)
            {
                if (pNodeStressesVector[i].NElems == 0) pNodeStressesVector[i].NElems = 1; // node surrounded by inactive elems
    
                pNodeStressesVector[i].EffectiveStressVector = pNodeStressesVector[i].EffectiveStressVector / pNodeStressesVector[i].NElems;
            }

            // Loop over nodes to assign the value of NODAL_STRESS_VECTOR
            for(ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it)
            {
                int Id = (*it).Id();
                //(*it).SetValue(NODAL_STRESS_VECTOR, pNodeStressesVector[Id-1].EffectiveStressVector);
				Vector& nodal_stress = it->GetSolutionStepValue(NODAL_STRESS_VECTOR);
				nodal_stress = pNodeStressesVector[Id - 1].EffectiveStressVector;
            }
        }

        // --------------------------------------------------------------------
        void WriteAuxPost(const NodeStresses* NodeStressesVector)
        {
            std::string eletyp = "Triangle";
            // Writing Post Mesh File
            std::fstream meshfile;
            meshfile.open( "_AMR_parameters.post.msh",std::fstream::out);
        
            meshfile << "MESH \"AMR_Mesh\" dimension " << 2 << " ElemType " << eletyp << " Nnode " << (*(mr_model_part.Elements().ptr_begin()))->GetGeometry().PointsNumber() << std::endl;
            meshfile << "Coordinates" << std::endl;
            
            for(ModelPart::NodeIterator i = mr_model_part.NodesBegin(); i != mr_model_part.NodesEnd(); ++i)
                meshfile << i->Id() << " " << i->X0() << " " << i->Y0() << std::endl;
            
            meshfile << "End Coordinates" << std::endl << std::endl;
        
            meshfile << "Elements" << std::endl;
        
            if(eletyp == "Triangle")
            {
                for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
                {
                    meshfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                            << (*it)->GetGeometry().GetPoint(2).Id() << std::endl;
                }
            }
            else
            {
                for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
                {
                    meshfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                            << (*it)->GetGeometry().GetPoint(2).Id() << " " << (*it)->GetGeometry().GetPoint(3).Id() << std::endl;
                }
            }
        
            meshfile << "End Elements" << std::endl << std::endl;
        
            meshfile.close();
        
            // Writing Results File
            std::fstream resfile;
            resfile.open("_AMR_parameters.post.res",std::fstream::out);
        
            resfile << "GiD Post Results File 1.0" << std::endl;
        
            resfile << "GaussPoints \"AMR_PostProcess\" Elemtype " << eletyp << std::endl;
            resfile << "Number of Gauss Points: " << 1 << std::endl;
            resfile << "Natural Coordinates: Internal" << std::endl;
            resfile << "End GaussPoints" << std::endl << std::endl;
            
            resfile << "Result \"Nodal_Effective_Stresses\" \"Kratos_AMR\" " << 1 << " Vector OnNodes" << std::endl;
            resfile << "ComponentNames \"Sxx\", \"Syy\", \"Sxy\"" << std::endl;
            resfile << "Values" << std::endl;
            for(unsigned int i = 0; i < mNNodes; i++)
            {
            resfile << i+1 << " " << NodeStressesVector[i].EffectiveStressVector[0] << " " << NodeStressesVector[i].EffectiveStressVector[1] << " " << NodeStressesVector[i].EffectiveStressVector[2] << std::endl;
            }
            resfile << "End Values" << std::endl << std::endl;
            
            resfile.close();
        }






//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    protected:
        
        // Member Variables
        ModelPart& mr_model_part;
        unsigned int mNNodes,mNElements;

    }; // Class

} // namespace Kratos

#endif /* KRATOS_STRESS_TO_NODES_PROCESS defined */