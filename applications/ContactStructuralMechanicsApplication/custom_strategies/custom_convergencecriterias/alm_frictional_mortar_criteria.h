// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "utilities/color_utilities.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
///@{
    
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**  
 * @class ALMFrictionalMortarConvergenceCriteria 
 * @ingroup ContactStructuralMechanicsApplication 
 * @brief Custom convergence criteria for the mortar condition for frictional case
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class ALMFrictionalMortarConvergenceCriteria 
    : public  BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ALMFrictionalMortarConvergenceCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaBaseType;
    
    typedef BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >          BaseType;

    typedef TSparseSpace                                                 SparseSpaceType;

    typedef typename BaseType::TDataType                                       TDataType;

    typedef typename BaseType::DofsArrayType                               DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType                                 NodesArrayType;
    
    typedef TableStreamUtility::Pointer                          TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ALMFrictionalMortarConvergenceCriteria(        
        TablePrinterPointerType pTable = nullptr,
        const bool PrintingOutput = false,
        const bool GiDIODebug = false
        ) : BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >(GiDIODebug),
        mpTable(pTable),
        mPrintingOutput(PrintingOutput),
        mTableIsInitialized(false)
    {
    }

    ///Copy constructor 
    ALMFrictionalMortarConvergenceCriteria( ALMFrictionalMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
      ,mpTable(rOther.mpTable)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    /// Destructor
    ~ALMFrictionalMortarConvergenceCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{
    
    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        BaseType::PreCriteria(rModelPart, rDofSet, A, Dx, b);
        
        return true;
    }
    
    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // We call the base class
        BaseType::PostCriteria(rModelPart, rDofSet, A, Dx, b);
        
        // Defining the convergence
        unsigned int is_converged_active = 0;
        unsigned int is_converged_slip = 0;
        
//         const double epsilon = rModelPart.GetProcessInfo()[INITIAL_PENALTY]; 
        const double scale_factor = rModelPart.GetProcessInfo()[SCALE_FACTOR];
        const double tangent_factor = rModelPart.GetProcessInfo()[TANGENT_FACTOR];
        
        const array_1d<double,3> zero_vector = ZeroVector(3);
        
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();

        #pragma omp parallel for reduction(+:is_converged_active, is_converged_slip)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            
            const double epsilon = it_node->GetValue(INITIAL_PENALTY); 
            
            const array_1d<double,3>& lagrange_multiplier = it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            const array_1d<double,3>& nodal_normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double normal_lagrange_multiplier = inner_prod(nodal_normal, lagrange_multiplier);
            
            const double augmented_normal_pressure = scale_factor * normal_lagrange_multiplier + epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);     
            
            it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)
            
            if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                if (it_node->Is(ACTIVE) == false ) {
                    it_node->Set(ACTIVE, true);
                    is_converged_active += 1;
                }
                
                // Computing the augmented tangent pressure
                const array_1d<double,3> tangent_lagrange_multiplier = lagrange_multiplier - normal_lagrange_multiplier * nodal_normal;
                const double lambda_tangent = norm_2(tangent_lagrange_multiplier); 
                
                // The friction coefficient
                const double mu = it_node->GetValue(FRICTION_COEFFICIENT);
                
                // Finally we compute the augmented tangent pressure
                const double gt = it_node->FastGetSolutionStepValue(WEIGHTED_SLIP);
                const double augmented_tangent_pressure = std::abs(scale_factor * lambda_tangent + tangent_factor * epsilon * gt) + mu * augmented_normal_pressure;
                
                it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_tangent_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)
                
                if (augmented_tangent_pressure <= 0.0) { // TODO: Check if it is minor equal or just minor
                    if (it_node->Is(SLIP) == true ) {
                        it_node->Set(SLIP, false);
                        is_converged_slip += 1;
                    }
                } else {
                    if (it_node->Is(SLIP) == false) {
                        it_node->Set(SLIP, true);
                        is_converged_slip += 1;
                    }
                }   
            } else {
                if (it_node->Is(ACTIVE) == true ) {
                    it_node->Set(ACTIVE, false);
                    is_converged_active += 1;
                }
            }
        }
        
        // We save to the process info if the active set has converged
        const bool active_set_converged = (is_converged_active + is_converged_slip) == 0 ? true : false;
        rModelPart.GetProcessInfo()[ACTIVE_SET_CONVERGED] = active_set_converged;
        
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (mpTable != nullptr) {
                auto& table = mpTable->GetTable();
                if (is_converged_active == 0) {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FGRN("       Achieved"));
                    else
                        table << "Achieved";
                } else {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FRED("   Not achieved"));
                    else
                        table << "Not achieved";
                }
                if (is_converged_slip == 0) {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FGRN("       Achieved"));
                    else
                        table << "Achieved";
                } else {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FRED("   Not achieved"));
                    else
                        table << "Not achieved";
                }
            } else {
                if (is_converged_active == 0) {
                    if (mPrintingOutput == false)
                        std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    else
                        std::cout << "\tActive set convergence is achieved" << std::endl;
                } else {
                    if (mPrintingOutput == false)
                        std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        std::cout << "\tActive set convergence is not achieved" << std::endl;
                }
                
                if (is_converged_slip == 0) {
                    if (mPrintingOutput == false)
                        std::cout << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    else
                        std::cout << "\tSlip/stick set convergence is achieved" << std::endl;
                } else {
                    if (mPrintingOutput == false)
                        std::cout << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        std::cout << "\tSlip/stick set  convergence is not achieved" << std::endl;
                }
            }
        }
        
        return active_set_converged;
    }
    
    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        ConvergenceCriteriaBaseType::mConvergenceCriteriaIsInitialized = true;
        
        if (mpTable != nullptr && mTableIsInitialized == false) {
            auto& table = mpTable->GetTable();
            table.AddColumn("ACTIVE SET CONV", 15);
            table.AddColumn("SLIP/STICK CONV", 15);
            mTableIsInitialized = true;
        }
    }
    
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    
    /**
     * @brief This method resets the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    
    void ResetWeightedGap(ModelPart& rModelPart) override
    {       
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, nodes_array);
        VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_SLIP, 0.0, nodes_array);
    }
    
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
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{
    
    TablePrinterPointerType mpTable; /// Pointer to the fancy table 
    bool mPrintingOutput;            /// If the colors and bold are printed
    bool mTableIsInitialized;        /// If the table is already initialized
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class ALMFrictionalMortarConvergenceCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H  defined */

