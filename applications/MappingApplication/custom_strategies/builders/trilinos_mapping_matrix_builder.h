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

#if !defined(KRATOS_TRILINOS_MAPPING_MATRIX_BUILDER_H_INCLUDED)
#define KRATOS_TRILINOS_MAPPING_MATRIX_BUILDER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"

namespace Kratos
{
///@addtogroup MappingApplication
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

/// Short class definition.
/** Detail class definition.
  */
template <class TSparseSpace,
          class TDenseSpace, // = DenseSpace<double>,
          >
class TrilinosMappingMatrixBuilder : public MappingMatrixBuilder<TSparseSpace, TDenseSpace>
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosMappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMappingMatrixBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosMappingMatrixBuilder() {}

    /// Destructor.
    virtual ~TrilinosMappingMatrixBuilder() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        std::stringstream buffer;
        buffer << "TrilinosMappingMatrixBuilder";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "TrilinosMappingMatrixBuilder"; }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosMappingMatrixBuilder &operator=(TrilinosMappingMatrixBuilder const &rOther) {}

    /// Copy constructor.
    TrilinosMappingMatrixBuilder(TrilinosMappingMatrixBuilder const &rOther) {}

    ///@}

}; // Class TrilinosMappingMatrixBuilder

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                TrilinosMappingMatrixBuilder &rThis) {}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const TrilinosMappingMatrixBuilder &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TRILINOS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
