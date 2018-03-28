//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"
//

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_CONDITION_H_INCLUDED)
#define KRATOS_NEAREST_NEIGHBOR_MAPPER_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{

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

class KRATOS_API(MAPPING_APPLICATION) NearestNeighborMapperCondition : public Condition
{
  public:
    ///@name Type Definitions
    ///@{

    typedef Condition BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapperCondition
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapperCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
   * Constructor.
   */
    NearestNeighborMapperCondition(IndexType NewId = 0);

    /**
   * Constructor using an array of nodes
   */
    NearestNeighborMapperCondition(IndexType NewId, const NodesArrayType &ThisNodes);

    /**
   * Constructor using Geometry
   */
    NearestNeighborMapperCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
   * Constructor using Properties
   */
    NearestNeighborMapperCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
   * Copy Constructor
   */
    NearestNeighborMapperCondition(NearestNeighborMapperCondition const &rOther);

    /**
   * Destructor
   */
    virtual ~NearestNeighborMapperCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NearestNeighborMapperCondition &operator=(NearestNeighborMapperCondition const &rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
   * CONDITIONS inherited from this class have to implement next
   * Create and Clone methods: MANDATORY
   */

    /**
   * creates a new condition pointer
   * @param NewId: the ID of the new condition
   * @param ThisNodes: the nodes of the new condition
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
   * creates a new condition pointer
   * @param NewId: the ID of the new condition
   * @param pGeom: the geometry to be employed
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;

    /**
   * creates a new condition pointer and clones the previous condition data
   * @param NewId: the ID of the new condition
   * @param ThisNodes: the nodes of the new condition
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

    /**
    * this determines the condition equation ID vector for all conditional
    * DOFs
    * @param rResult: the condition equation ID vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);


    /**
   * CONDITIONS inherited from this class have to implement next
   * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
   * they can be managed internally with a private method to do the same calculations
   * only once: MANDATORY
   */


    /**
   * this is called during the assembling process in order
   * to calculate the condition left hand side matrix only
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
    virtual void CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo);



    /**
   * This method provides the place to perform checks on the completeness of the input
   * and the compatibility with the problem options as well as the contitutive laws selected
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   * this method is: MANDATORY
   */
    virtual int Check(const ProcessInfo &rCurrentProcessInfo);

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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer &rSerializer) const;
    virtual void load(Serializer &rSerializer);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class NearestNeighborMapperCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_CONDITION_H_INCLUDED  defined
