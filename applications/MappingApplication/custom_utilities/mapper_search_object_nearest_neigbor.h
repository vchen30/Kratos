//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

#error Please change the author fullname and if your license is different please state it above
#error Please change the KRATOS_FILENAME_H_INCLUDED and InterfaceObjectNearestNeighbor and remove these error pragmas

#if !defined(KRATOS_FILENAME_H_INCLUDED )
#define  KRATOS_FILENAME_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
  class NearestNeighborSearchObject : public MapperSearchObject
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceObjectNearestNeighbor
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectNearestNeighbor);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      InterfaceObjectNearestNeighbor();

      /// Destructor.
      virtual ~InterfaceObjectNearestNeighbor();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

        void ProcessSearchResult(InterfaceBaseObject::iterator OtherSearchObject, double Distance)
        {
            mSearchSuccessful = true;

            if (Distance < mClosestNeighborDistance)
            {
                mClosestNeighborDistance = Distance;
                mNeighborMappingMatrixEquationId = OtherSearchObject->pGetBaseNode()->GetValue(MAPPING_MATRIX_EQUATION_ID);
            }
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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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

        int mNeighborMappingMatrixEquationId;
        double mClosestNeighborDistance;


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
      InterfaceObjectNearestNeighbor& operator=(InterfaceObjectNearestNeighbor const& rOther);

      /// Copy constructor.
      InterfaceObjectNearestNeighbor(InterfaceObjectNearestNeighbor const& rOther);


      ///@}

    }; // Class InterfaceObjectNearestNeighbor

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceObjectNearestNeighbor& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceObjectNearestNeighbor& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
