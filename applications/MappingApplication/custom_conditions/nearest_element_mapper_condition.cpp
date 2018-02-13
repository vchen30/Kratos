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

// System includes

// External includes

// Include Base h
#include "custom_conditions/nearest_element_mapper_condition.h"

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

/**
 * Constructor.
 */
NearestElementMapperCondition::NearestElementMapperCondition(IndexType NewId)
    : BaseMapperCondition(NewId) // @{KRATOS_INIT_MEMBER_LIST}
{
}

/**
 * Constructor using an array of nodes
 */
NearestElementMapperCondition::NearestElementMapperCondition(IndexType NewId, const NodesArrayType &ThisNodes)
    : BaseMapperCondition(NewId, ThisNodes) // @{KRATOS_INIT_MEMBER_LIST}
{
}

/**
 * Constructor using Geometry
 */
NearestElementMapperCondition::NearestElementMapperCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseMapperCondition(NewId, pGeometry) // @{KRATOS_INIT_MEMBER_LIST}
{
}

/**
 * Constructor using Properties
 */
NearestElementMapperCondition::NearestElementMapperCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseMapperCondition(NewId, pGeometry, pProperties) // @{KRATOS_INIT_MEMBER_LIST}
{
}

/**
 * Copy Constructor
 */
NearestElementMapperCondition::NearestElementMapperCondition(NearestElementMapperCondition const &rOther)
    : BaseMapperCondition(rOther) // @{KRATOS_CC_INIT_MEMBER_LIST}
{
}

/**
 * Destructor
 */
NearestElementMapperCondition::~NearestElementMapperCondition()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
NearestElementMapperCondition &NearestElementMapperCondition::operator=(NearestElementMapperCondition const &rOther)
{
    BaseType::operator=(rOther);
    Flags::operator=(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

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
Condition::Pointer NearestElementMapperCondition::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Condition::Pointer(new NearestElementMapperCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestElementMapperCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Condition::Pointer(new NearestElementMapperCondition(NewId, pGeom, pProperties));
    KRATOS_CATCH("");
}

/**
 * creates a new condition pointer and clones the previous condition data
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestElementMapperCondition::Clone(IndexType NewId, NodesArrayType const &ThisNodes) const
{
    KRATOS_TRY
    return Condition::Pointer(new NearestElementMapperCondition(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
    KRATOS_CATCH("");
}

/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &CurrentProcessInfo)
{
    const SizeType num_points = GetGeometry().PointsNumber();

    if (rResult.size() != 2*num_points)
        rResult.resize(2*num_points, false);

    for (SizeType i=0; i < num_points; ++i)
    {
        rResult[i] = GetGeometry().GetPoint(i).GetValue(MAPPING_MATRIX_EQUATION_ID); // ID on Destination
        // rResult[i+1] = GetGeometry().GetPoint(i).GetValue(MAPPER_NEIGHBOR_INFORMATION)[0]; // ID on Origin. This is written by the Communicator
    }

    rResult[0] = GetGeometry().GetPoint(0).GetValue(MAPPING_MATRIX_EQUATION_ID); // ID on Destination

    // rResult[1] = GetGeometry().GetPoint(0).GetValue(MAPPER_NEIGHBOR_INFORMATION)[0]; // ID on Origin. This is written by the Communicator
}

/**
 * determines the condition equation list of DOFs
 * @param ConditionDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::GetDofList(DofsVectorType &rConditionDofList, ProcessInfo &CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rConditionDofList.size() != number_of_nodes)
        rConditionDofList.resize(number_of_nodes);

    // @{ KRATOS_CONDITION_LIST_DOFS }
}

/**
 * CONDITIONS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all condition contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix only
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector only
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateFirstDerivativesContributions(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{

    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateFirstDerivativesLHS(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateFirstDerivativesRHS(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * CONDITION inherited from this class must implement this methods
 * if they need to add dynamic condition contributions
 * note: second derivatives means the accelerations if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateSecondDerivativesContributions,
 * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
 */

/**
 * this is called during the assembling process in order
 * to calculate the second derivative contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateSecondDerivativesContributions(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{

    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateSecondDerivativesLHS(
    MatrixType &rLeftHandSideMatrix,
    ProcessInfo &rCurrentProcessInfo)
{

    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateSecondDerivativesRHS(
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{

    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition mass matrix
 * @param rMassMatrix: the condition mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 0)
        rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition damping matrix
 * @param rDampingMatrix: the condition damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementMapperCondition::CalculateDampingMatrix(MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 0)
        rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int NearestElementMapperCondition::Check(const ProcessInfo &rCurrentProcessInfo)
{

    KRATOS_TRY

    if (this->Id() < 1)
    {
        KRATOS_THROW_ERROR(std::logic_error, "NearestElementMapperCondition found with Id 0 or negative", "")
    }

    if (this->GetGeometry().Area() <= 0)
    {
        std::cout << "error on NearestElementMapperCondition -> " << this->Id() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0", "")
    }

    return 0;

    KRATOS_CATCH("");
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

std::string NearestElementMapperCondition::Info() const
{
    std::stringstream buffer;
    buffer << "NearestElementMapperCondition #" << Id();
    return buffer.str();
}

/// Print information about this object.

void NearestElementMapperCondition::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "NearestElementMapperCondition #" << Id();
}

/// Print object's data.

void NearestElementMapperCondition::PrintData(std::ostream &rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

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

void NearestElementMapperCondition::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseMapperCondition);

    // List
    // To be completed with the class member list
}

void NearestElementMapperCondition::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseMapperCondition);

    // List
    // To be completed with the class member list
}

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
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream, NearestElementMapperCondition &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream, const NearestElementMapperCondition &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
