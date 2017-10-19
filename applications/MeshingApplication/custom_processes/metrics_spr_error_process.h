// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Anna Rehr
//

#if !defined(KRATOS_SPR_ERROR_METRICS_PROCESS)
#define KRATOS_SPR_ERROR_METRICS_PROCESS

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "meshing_application.h"
#include "processes/find_nodal_neighbours_process.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "spaces/ublas_space.h"
namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    typedef Node <3>                                                                NodeType;
    
///@}
///@name  Enum's
///@{
    
    #if !defined(INTERPOLATION_METRIC)
    #define INTERPOLATION_METRIC
        enum Interpolation {Constant = 0, Linear = 1, Exponential = 2};
    #endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute the metrics of the model part with an Hessian approach

//template<unsigned int TDim, class TVarType>  
template<unsigned int TDim> 
class ComputeSPRErrorSolMetricProcess
    //: public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ComputeSPRErrorSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeSPRErrorSolMetricProcess);
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param ThisParameters: The input parameters
     */
    
    ComputeSPRErrorSolMetricProcess(
        ModelPart& rThisModelPart,
        //TVarType& rVariable,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mThisModelPart(rThisModelPart)
        //:mThisModelPart(rThisModelPart),
        //mVariable(rVariable)
    {               
        Parameters DefaultParameters = Parameters(R"(
        {
            "minimal_size"                        : 0.1,
            "maximal_size"                        : 10.0, 
            "error"                               : 0.1,
            "echo_level"                          : 0
        })" );
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
        mMinSize = ThisParameters["minimal_size"].GetDouble();
        mMaxSize = ThisParameters["maximal_size"].GetDouble();
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        
    }
    
    /// Destructor.
    virtual ~ComputeSPRErrorSolMetricProcess() {}
    
    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     */
    
    virtual double Execute()
    {
        return SuperconvergentPatchRecovery();

    }
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ComputeSPRErrorSolMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeSPRErrorSolMetricProcess";
    }

    /// Print object"s data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
    
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
    
private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{
    
    ModelPart& mThisModelPart;               // The model part to compute
    //TVarType mVariable;            // The variable to calculate the hessian
    double mMinSize;                         // The minimal size of the elements
    double mMaxSize;                         // The maximal size of the elements
    bool mEnforceCurrent;                    // With this we choose if we inforce the current nodal size (NODAL_H)
    double mInterpError;                     // The error of interpolation allowed
    double mMeshConstant;                    // The mesh constant to remesh (depends of the element type)
    double mAnisRatio;                       // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;                      // The boundary layer limit distance
    Interpolation mInterpolation;            // The interpolation type
    int mEchoLevel;
    

    double SuperconvergentPatchRecovery()
    {
        /************************************************************************
        --1-- calculate superconvergent stresses (at the nodes) --1--
        ************************************************************************/
        if(mEchoLevel>2){
            std::vector<std::string> submodels;
            submodels= mThisModelPart.GetSubModelPartNames();
            for (std::vector<std::string>::const_iterator i = submodels.begin();i!=submodels.end();i++) 
                std::cout << *i<<std::endl;}
        FindNodalNeighboursProcess findNeighbours(mThisModelPart);
        findNeighbours.Execute();
        //std::vector<Vector> stress_vector(1);
        //std::vector<array_1d<double,3>> coordinates_vector(1);
        //ariable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
        //iteration over all nodes -- construction of patches
        ModelPart::NodesContainerType& rNodes = mThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            int neighbour_size = i_nodes->GetValue(NEIGHBOUR_ELEMENTS).size();
            //std::cout << "Node: " << i_nodes->Id() << " has " << neighbour_size << " neighbouring elements: " << std::endl;

            Vector sigma_recovered(3,0);
            if(neighbour_size>2){ 
                CalculatePatch(i_nodes,i_nodes,neighbour_size,sigma_recovered);
                i_nodes->SetValue(RECOVERED_STRESS,sigma_recovered);
                if(mEchoLevel>2)
                    std::cout<<"recovered sigma"<<sigma_recovered<<std::endl;
            }
            else{
                for(WeakPointerVector< Node<3> >::iterator i_neighbour_nodes = i_nodes->GetValue(NEIGHBOUR_NODES).begin(); i_neighbour_nodes != i_nodes->GetValue(NEIGHBOUR_NODES).end(); i_neighbour_nodes++){
                    Vector sigma_recovered_i(3);
                    unsigned int count_i=0;
                    for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++){
                        if (i->Id() == i_neighbour_nodes->Id() && i->GetValue(NEIGHBOUR_ELEMENTS).size()>2){
                            CalculatePatch(i_nodes,i,neighbour_size,sigma_recovered_i);
                            count_i ++;
                        }
                    }
                    //average solution from different patches
                    if(count_i != 0)
                        sigma_recovered =sigma_recovered*(count_i-1)/count_i + sigma_recovered_i/count_i;
                }
                i_nodes->SetValue(RECOVERED_STRESS,sigma_recovered);
            }
       }
        /******************************************************************************
        --2-- calculate error estimation and new element size (for each element) --2--
        ******************************************************************************/
        //loop over all elements: 
        double error_overall=0;
        double energy_norm_overall=0;

        //compute the error estimate per element
        for(ModelPart::ElementsContainerType::iterator i_elements = mThisModelPart.Elements().begin() ; i_elements != mThisModelPart.Elements().end(); i_elements++) 
        {
            std::vector<double> error_integration_point;
            i_elements->GetValueOnIntegrationPoints(ERROR_INTEGRATION_POINT,error_integration_point,mThisModelPart.GetProcessInfo());
            double error_energy_norm=0;
            for(unsigned int i=0;i<error_integration_point.size();i++)
                error_energy_norm += error_integration_point[i];
            error_overall += error_energy_norm;
            error_energy_norm= sqrt(error_energy_norm);
            i_elements->SetValue(ELEMENT_ERROR,error_energy_norm);
            if(mEchoLevel>2)
                std::cout<<"element_error:"<<error_energy_norm<<std::endl;


            std::vector<double> strain_energy;
            i_elements->GetValueOnIntegrationPoints(STRAIN_ENERGY,strain_energy,mThisModelPart.GetProcessInfo());
            double energy_norm=0;
            for(unsigned int i=0;i<strain_energy.size();i++)
                energy_norm += 2*strain_energy[i];
            energy_norm_overall += energy_norm;
            energy_norm= sqrt(energy_norm);
            if(mEchoLevel>2)
                std::cout<<"energy norm:"<<energy_norm<<std::endl;
        }
        error_overall = sqrt(error_overall);
        energy_norm_overall = sqrt(energy_norm_overall);
        double error_percentage = error_overall/pow((error_overall*error_overall+energy_norm_overall*energy_norm_overall),0.5);
        if(mEchoLevel>1){
            std::cout<<"overall error norm :"<<error_overall<<std::endl;
            std::cout<<"overall energy norm :"<<energy_norm_overall<<std::endl;
            std::cout<<"error in percent: "<<error_percentage<<std::endl;}
        
        //compute new element size
        for(ModelPart::ElementsContainerType::iterator i_elements = mThisModelPart.Elements().begin() ; i_elements != mThisModelPart.Elements().end(); i_elements++) 
        {
            //compute the current element size h
            //i_elements->CalculateElementSize();
            ComputeElementSize(i_elements);

            //compute new element size
            double new_element_size;
            //old:
            new_element_size = i_elements->GetValue(ELEMENT_H)/i_elements->GetValue(ELEMENT_ERROR);
            //new_element_size *= sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/mThisModelPart.Elements().size())*0.1;
            new_element_size *= sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/1000)*0.15;
            
            //new: experiment
            //new_element_size = i_elements->GetValue(ELEMENT_H)*i_elements->GetValue(ELEMENT_H)/i_elements->GetValue(ELEMENT_ERROR);
            //new_element_size *= sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/mThisModelPart.Elements().size())*0.15;
            //new_element_size = sqrt(new_element_size);
            //std::cout<<"old element size: "<<i_elements->GetValue(ELEMENT_H)<<std::endl;
            i_elements->SetValue(ELEMENT_H,new_element_size);
            //std::cout<<"new element size: "<<i_elements->GetValue(ELEMENT_H)<<std::endl;
        }

        /******************************************************************************
        --3-- calculate metric (for each node) --3--
        ******************************************************************************/

        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            // get maximal element size from neighboring elements
            double h_min=0;
            for(WeakPointerVector< Element >::iterator i_neighbour_elements = i_nodes->GetValue(NEIGHBOUR_ELEMENTS).begin(); i_neighbour_elements != i_nodes->GetValue(NEIGHBOUR_ELEMENTS).end(); i_neighbour_elements++){
                if(h_min==0||h_min>i_neighbour_elements->GetValue(ELEMENT_H))
                    h_min = i_neighbour_elements->GetValue(ELEMENT_H);
                
            }
            //std::cout<<"h_min: "<<h_min<<std::endl;
            // set metric
            Matrix metric_matrix(2,2,0);
            metric_matrix(0,0)=1/(h_min*h_min);
            metric_matrix(1,1)=1/(h_min*h_min);
            // transform metric matrix to a vector
            const Vector metric = MetricsMathUtils<TDim>::TensorToVector(metric_matrix);
            i_nodes->SetValue(MMG_METRIC,metric);

            if(mEchoLevel>2)
                std::cout<<"node "<<i_nodes->Id()<<" has metric: "<<i_nodes->GetValue(MMG_METRIC)<<std::endl;
        }
        return error_overall/pow((error_overall*error_overall+energy_norm_overall*energy_norm_overall),0.5);
    }

    void CalculatePatch(
        ModelPart::NodesContainerType::iterator i_nodes,
        ModelPart::NodesContainerType::iterator i_patch_node,
        int neighbour_size,
        Vector& rsigma_recovered)
    {
        // determine if contact BC has to be regarded
        bool regard_contact = (i_nodes->GetValue(CONTACT_PRESSURE) != 0.0);
        //regard_contact = i_patch_node->Has(CONTACT_PRESSURE);
        /*if(regard_contact == false)
        {
            for( auto i_neighbour_nodes = i_patch_node->GetValue(NEIGHBOUR_NODES).begin(); i_neighbour_nodes != i_patch_node->GetValue(NEIGHBOUR_NODES).end(); i_neighbour_nodes++) {
                if (i_neighbour_nodes->Has(CONTACT_PRESSURE))
                {
                    regard_contact = true;
                    break;
                }
            }
        }*/
        if (regard_contact == false)
            CalculatePatchStandard(i_nodes, i_patch_node, neighbour_size, rsigma_recovered);
        else
            CalculatePatchContact(i_nodes, i_patch_node, neighbour_size, rsigma_recovered);

    }

    
    //calculates the recovered stress at a node in the case of a standard patch without contact BC
    // i_node: the node for which the recovered stress should be calculated
    // i_patch_node: the center node of the patch
    void CalculatePatchStandard(
        ModelPart::NodesContainerType::iterator i_nodes,
        ModelPart::NodesContainerType::iterator i_patch_node,
        int neighbour_size,
        Vector& rsigma_recovered)
    {
        std::vector<Vector> stress_vector(1);
        std::vector<array_1d<double,3>> coordinates_vector(1);
        Variable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
        Variable<Vector> variable_stress = CAUCHY_STRESS_VECTOR;
        Matrix A(3,3,0);
        Matrix b(3,3,0); 
        Matrix p_k(1,3,0);
        for( WeakPointerVector< Element >::iterator i_elements = i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).begin(); i_elements != i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).end(); i_elements++) {
            //std::cout << "\tElement: " << i_elements->Id() << std::endl;
            i_elements->GetValueOnIntegrationPoints(variable_stress,stress_vector,mThisModelPart.GetProcessInfo());
            i_elements->GetValueOnIntegrationPoints(variable_coordinates,coordinates_vector,mThisModelPart.GetProcessInfo());

            //std::cout << "\tstress: " << stress_vector[0] << std::endl;
            //std::cout << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;
            Matrix sigma(1,3);
            for(int j=0;j<3;j++)
                sigma(0,j)=stress_vector[0][j];
            p_k(0,0)=1;
            p_k(0,1)=coordinates_vector[0][0]-i_patch_node->X(); 
            p_k(0,2)=coordinates_vector[0][1]-i_patch_node->Y();   
            A+=prod(trans(p_k),p_k);
            b+=prod(trans(p_k),sigma);
        }
        Matrix invA(3,3);
        double det;
        MathUtils<double>::InvertMatrix(A,invA,det);
        //std::cout <<A<<std::endl;
        //std::cout <<invA<<std::endl;
        //std::cout << det<< std::endl;

        Matrix coeff(3,3);
        coeff = prod(invA,b);
        if(neighbour_size > 2)
            rsigma_recovered = MatrixRow(coeff,0);
        else{
            p_k(0,1)=i_nodes->X()-i_patch_node->X(); 
            p_k(0,2)=i_nodes->Y()-i_patch_node->Y();
            Matrix sigma(1,3);
            sigma = prod(p_k,coeff);
            rsigma_recovered = MatrixRow(sigma,0);
        }
        if(i_nodes->Id() == 1132)
            std::cout<<"strange, no contact, contact pressure: "<<i_nodes->GetValue(CONTACT_PRESSURE)<<std::endl;
    }

    //calculates the recovered stress at a node where contact BCs are regarded
    // i_node: the node for which the recovered stress should be calculated
    // i_patch_node: the center node of the patch

    void CalculatePatchContact(
        ModelPart::NodesContainerType::iterator i_nodes,
        ModelPart::NodesContainerType::iterator i_patch_node,
        int neighbour_size,
        Vector& rsigma_recovered)
    {
        if(mEchoLevel>1)
            std::cout<<"contact i regarded"<<std::endl;
        std::vector<Vector> stress_vector(1);
        std::vector<array_1d<double,3>> coordinates_vector(1);
        Variable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
        Variable<Vector> variable_stress = CAUCHY_STRESS_VECTOR;
        CompressedMatrix A(9,9,0);
        Matrix b(9,1,0); 
        Matrix p_k(3,9,0);
        Matrix N_k(1,3,0);
        Matrix T_k(1,3,0);
        double penalty_normal = 10000;
        double penalty_tangential = 10000;
        Matrix sigma(3,1);
        // computation A and b
        // PART 1: contributions from the neighboring elements
        for( WeakPointerVector< Element >::iterator i_elements = i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).begin(); i_elements != i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).end(); i_elements++) {
            //std::cout << "\tElement: " << i_elements->Id() << std::endl;
            i_elements->GetValueOnIntegrationPoints(variable_stress,stress_vector,mThisModelPart.GetProcessInfo());
            i_elements->GetValueOnIntegrationPoints(variable_coordinates,coordinates_vector,mThisModelPart.GetProcessInfo());

            //std::cout << "\tstress: " << stress_vector[0] << std::endl;
            //std::cout << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;
            for(int j=0;j<3;j++)
                sigma(j,0)=stress_vector[0][j];
            p_k(0,0)=1;
            p_k(0,1)=coordinates_vector[0][0]-i_patch_node->X(); 
            p_k(0,2)=coordinates_vector[0][1]-i_patch_node->Y();
            p_k(1,3)=1;
            p_k(1,4)=coordinates_vector[0][0]-i_patch_node->X(); 
            p_k(1,5)=coordinates_vector[0][1]-i_patch_node->Y();
            p_k(2,6)=1;
            p_k(2,7)=coordinates_vector[0][0]-i_patch_node->X(); 
            p_k(2,8)=coordinates_vector[0][1]-i_patch_node->Y();
               
            A+=prod(trans(p_k),p_k);
            b+=prod(trans(p_k),sigma);
        }
        // computing A and b
        Matrix A1(9,1,0), A2(1,9,0);
        
        p_k(0,1)= i_nodes->X()-i_patch_node->X();
        p_k(0,2)= i_nodes->Y()-i_patch_node->Y();
        p_k(1,4)= i_nodes->X()-i_patch_node->X();;
        p_k(1,5)= i_nodes->Y()-i_patch_node->Y();
        p_k(2,7)= i_nodes->X()-i_patch_node->X();;
        p_k(2,8)= i_nodes->Y()-i_patch_node->Y();
        N_k(0,0) = i_nodes->GetValue(NORMAL)[0]*i_nodes->GetValue(NORMAL)[0];
        N_k(0,1) = i_nodes->GetValue(NORMAL)[1]*i_nodes->GetValue(NORMAL)[1];
        N_k(0,2) = 2*i_nodes->GetValue(NORMAL)[0]*i_nodes->GetValue(NORMAL)[1];
        T_k(0,0) = i_nodes->GetValue(NORMAL)[0]*i_nodes->GetValue(NORMAL)[1];
        T_k(0,1) = -i_nodes->GetValue(NORMAL)[0]*i_nodes->GetValue(NORMAL)[1];
        T_k(0,2) = i_nodes->GetValue(NORMAL)[1]*i_nodes->GetValue(NORMAL)[1]-i_nodes->GetValue(NORMAL)[0]*i_nodes->GetValue(NORMAL)[0];

        A1 = prod(trans(p_k),trans(N_k));
        A2 = prod(N_k,p_k);
        A+= penalty_normal*prod(A1, A2);

        A1 = prod(trans(p_k),trans(T_k));
        A2 = prod(T_k,p_k);
        A+= penalty_tangential*prod(A1, A2);

        b= penalty_normal*prod(trans(p_k),trans(N_k))*i_nodes->GetValue(CONTACT_PRESSURE);
        /*
        //PART 2: contributions from contact nodes: regard all nodes from the patch which are in contact
        //patch center node:
        if (i_patch_node->Has(CONTACT_PRESSURE)){
            p_k(0,1)=0;
            p_k(0,2)=0;
            p_k(1,4)=0;
            p_k(1,5)=0;
            p_k(2,7)=0;
            p_k(2,8)=0;
            N_k(0,0) = i_patch_node->GetValue(NORMAL)[0]*i_patch_node->GetValue(NORMAL)[0];
            N_k(0,1) = i_patch_node->GetValue(NORMAL)[1]*i_patch_node->GetValue(NORMAL)[1];
            N_k(0,2) = 2*i_patch_node->GetValue(NORMAL)[0]*i_patch_node->GetValue(NORMAL)[1];
            T_k(0,0) = i_patch_node->GetValue(NORMAL)[0]*i_patch_node->GetValue(NORMAL)[1];
            T_k(0,1) = -i_patch_node->GetValue(NORMAL)[0]*i_patch_node->GetValue(NORMAL)[1];
            T_k(0,2) = i_patch_node->GetValue(NORMAL)[1]*i_patch_node->GetValue(NORMAL)[1]-i_patch_node->GetValue(NORMAL)[0]*i_patch_node->GetValue(NORMAL)[0];

            A1 = prod(trans(p_k),trans(N_k));
            A2 = prod(N_k,p_k);
            A+= penalty_normal*prod(A1, A2);

            A1 = prod(trans(p_k),trans(T_k));
            A2 = prod(T_k,p_k);
            A+= penalty_tangential*prod(A1, A2);
            //A+= penalty_normal*prod(prod(trans(p_k),trans(N_k)),prod(N_k,p_k));
            //A+= penalty_tangential*prod(prod(prod(trans(p_k),trans(T_k)),T_k),p_k);

            b-= penalty_normal*prod(trans(p_k),trans(N_k))*i_patch_node->GetValue(CONTACT_PRESSURE);
        }

        //neighboring nodes:
        
        for( auto i_neighbour_nodes = i_patch_node->GetValue(NEIGHBOUR_NODES).begin(); i_neighbour_nodes != i_patch_node->GetValue(NEIGHBOUR_NODES).end(); i_neighbour_nodes++) {
            if (i_neighbour_nodes->Has(CONTACT_PRESSURE)){
                p_k(0,1)= i_neighbour_nodes->X()-i_patch_node->X();
                p_k(0,2)= i_neighbour_nodes->Y()-i_patch_node->Y();
                p_k(1,4)= i_neighbour_nodes->X()-i_patch_node->X();;
                p_k(1,5)= i_neighbour_nodes->Y()-i_patch_node->Y();
                p_k(2,7)= i_neighbour_nodes->X()-i_patch_node->X();;
                p_k(2,8)= i_neighbour_nodes->Y()-i_patch_node->Y();
                N_k(0,0) = i_neighbour_nodes->GetValue(NORMAL)[0]*i_neighbour_nodes->GetValue(NORMAL)[0];
                N_k(0,1) = i_neighbour_nodes->GetValue(NORMAL)[1]*i_neighbour_nodes->GetValue(NORMAL)[1];
                N_k(0,2) = 2*i_neighbour_nodes->GetValue(NORMAL)[0]*i_neighbour_nodes->GetValue(NORMAL)[1];
                T_k(0,0) = i_neighbour_nodes->GetValue(NORMAL)[0]*i_neighbour_nodes->GetValue(NORMAL)[1];
                T_k(0,1) = -i_neighbour_nodes->GetValue(NORMAL)[0]*i_neighbour_nodes->GetValue(NORMAL)[1];
                T_k(0,2) = i_neighbour_nodes->GetValue(NORMAL)[1]*i_neighbour_nodes->GetValue(NORMAL)[1]-i_neighbour_nodes->GetValue(NORMAL)[0]*i_neighbour_nodes->GetValue(NORMAL)[0];

                A1 = prod(trans(p_k),trans(N_k));
                A2 = prod(N_k,p_k);
                A+= penalty_normal*prod(A1, A2);

                A1 = prod(trans(p_k),trans(T_k));
                A2 = prod(T_k,p_k);
                A+= penalty_tangential*prod(A1, A2);

                b+= penalty_normal*prod(trans(p_k),trans(N_k))*i_neighbour_node->GetValue(CONTACT_PRESSURE);
            }
        }*/

        // computing coefficients a: A*a=b
        //UblasSpace<double,CompressedMatrix,Vector> U1 = UblasSpace<double, Matrix,Vector>();
        //UblasSpace<double, Matrix,Vector> U2 = UblasSpace<double, Matrix,Vector>();
        SkylineLUFactorizationSolver< UblasSpace<double,CompressedMatrix,Vector>, UblasSpace<double,Matrix,Vector>> solver = SkylineLUFactorizationSolver< UblasSpace<double,CompressedMatrix,Vector>, UblasSpace<double,Matrix,Vector>>();
        //std::cout<<A<<std::endl;
        /*CompressedMatrix compA(9,9,81);// = A.sparseView();
        for (unsigned i = 0; i < compA.size1 (); ++ i){
            for (unsigned j = 0; j < compA.size2 (); ++ j)
                compA (i, j) = 9 * i + j;
        }*/
        Vector coeff(9);
        Vector b_vector = MatrixColumn(b,0);
        solver.Solve(A,coeff,b_vector);
                        
        p_k(0,1)= i_nodes->X()-i_patch_node->X();
        p_k(0,2)= i_nodes->Y()-i_patch_node->Y();
        p_k(1,4)= i_nodes->X()-i_patch_node->X();;
        p_k(1,5)= i_nodes->Y()-i_patch_node->Y();
        p_k(2,7)= i_nodes->X()-i_patch_node->X();;
        p_k(2,8)= i_nodes->Y()-i_patch_node->Y();
        Matrix coeff_matrix(9,1);
        for (unsigned int i=0; i<9; i++)
            coeff_matrix(i,0)=coeff(i);
        sigma = prod(p_k,coeff_matrix);

        rsigma_recovered = MatrixColumn(sigma,0);
        if(mEchoLevel>1)
            std::cout<<"recovered pressure: "<<prod(N_k,sigma)<<", LM: "<<i_nodes->GetValue(CONTACT_PRESSURE)<<std::endl;
        if(i_nodes->Id() == 1132)
            std::cout<<"strange, contact, contact pressure: "<<i_nodes->GetValue(CONTACT_PRESSURE)<<std::endl;
        
    }

    // set the element size. The element size is defined as the diameter of the smallest ball containing the element
    void ComputeElementSize(ModelPart::ElementsContainerType::iterator pElement){

        // triangular elements
        if (pElement->GetGeometry().size()==3){
        pElement->SetValue(ELEMENT_H,2*pElement->GetGeometry().Circumradius());
        }
        
    }
    /// Assignment operator.
    ComputeSPRErrorSolMetricProcess& operator=(ComputeSPRErrorSolMetricProcess const& rOther);

    /// Copy constructor.
    //ComputeSPRErrorSolMetricProcess(ComputeSPRErrorSolMetricProcess const& rOther);

};// class ComputeSPRErrorSolMetricProcess

};// namespace Kratos.
#endif /* KRATOS_SPR_ERROR_METRICS_PROCESS defined */
