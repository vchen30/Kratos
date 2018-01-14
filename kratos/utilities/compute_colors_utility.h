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
#if !defined(KRATOS_COMPUTE_COLORS_UTILITY_H_INCLUDED)
#define KRATOS_COMPUTE_COLORS_UTILITY_H_INCLUDED

// System includes
#include <unordered_map>
#include <set>
#include <unordered_set>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/key_hash.h"

namespace Kratos {
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // Index definition
    typedef std::size_t                                         SizeType;
    
    // Maps definition
    typedef std::unordered_map<int,int>                     ColorMapType; 
    typedef std::unordered_map<int,std::vector<std::string>> NameMapType;
    
    // Containers definition
    typedef ModelPart::NodesContainerType                 NodesArrayType;
    typedef ModelPart::ElementsContainerType           ElementsArrayType;
    typedef ModelPart::ConditionsContainerType       ConditionsArrayType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief ComputeColorsUtility 
 * This utility can be used to create a "database" to organize the submodelparts and identify where belongs each component
 */
class KRATOS_API(KRATOS_CORE) ComputeColorsUtility {
   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeColorsUtility
    KRATOS_CLASS_POINTER_DEFINITION(ComputeColorsUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    /**
     * This is the default constructor
     */
    
    ComputeColorsUtility(ModelPart& rThisModelPart) :mrThisModelPart(rThisModelPart)
    {
    }

    /// Destructor.
    virtual ~ComputeColorsUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This functions gets the "colors", parts of a model part to process
     */
    
    void Initialize();
    
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
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{
    
    ModelPart& mrThisModelPart;                                    // The model part to compute  
    
    ColorMapType mNodesColors, mElementsColors, mConditionsColors; // This class stores the colors of the components
    
    NameMapType mColors;                                           // Where the sub model parts IDs are stored
    
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};  // Class ComputeColorsUtility
}
#endif /* KRATOS_COMPUTE_COLORS_UTILITY_H_INCLUDED defined */
