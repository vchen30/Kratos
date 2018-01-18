// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Klaus B. Sautter
//

// System includes
#if !defined(KRATOS_LINE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_CONDITION_H_INCLUDED


// Project includes
#include "includes/define.h"
#include "custom_conditions/base_load_condition.h"
#include "includes/variables.h"

namespace Kratos
{

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  LineLoadCondition
    : public BaseLoadCondition
{
public:

    /// Counted pointer of LineLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineLoadCondition );


    /// Default constructor.
    LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    ~LineLoadCondition() override;

    
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;
        
    Condition::Pointer Create( 
        IndexType NewId, 
        NodesArrayType const& ThisNodes, 
        PropertiesType::Pointer pProperties 
        ) const override;


    /*     Adds work equivalent nodal moments  */
    void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
        const Vector& rForceInput, VectorType& rRightHandSideVector) const;

    void AddLineLoad(VectorType& rRightHandSideVector);

    void AddLinePressure(VectorType& rRightHandSideVector,
         MatrixType& rLeftHandSideMatrix,
         bool CalculateStiffnessMatrixFlag,
         bool CalculateResidualVectorFlag);

    void UpdateMemberLength();

    void UpdateNodalPosition();

protected:
    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        ) override;

    void CalculateAndSubKp(
        Matrix& K,
        const Matrix& DN_De,
        const Vector& N,
        const double Pressure,
        const double IntegrationWeight
        );

    void CalculateAndAddPressureForce(
        VectorType& rRightHandSideVector,
        const Vector& N,
        const array_1d<double, 3>& Normal,
        const double Pressure,
        const double IntegrationWeight 
        );

    // A protected default constructor necessary for serialization
    LineLoadCondition() {};


private:

    double mL;
    int mTimeStep = 0;
    Vector mNodeA,mNodeB;

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
        rSerializer.save("MemberLength", this->mL);
        rSerializer.save("TimeStep", this->mTimeStep);
        rSerializer.save("PositionNodeA", this->mNodeA);
        rSerializer.save("PositionNodeB", this->mNodeB);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
        rSerializer.load("MemberLength", this->mL);
        rSerializer.load("TimeStep", this->mTimeStep);
        rSerializer.load("PositionNodeA", this->mNodeA);
        rSerializer.load("PositionNodeB", this->mNodeB);
    }



}; // Class LineLoadCondition


}  // namespace Kratos.

#endif // KRATOS_LINE_LOAD_CONDITION_H_INCLUDED  defined 


