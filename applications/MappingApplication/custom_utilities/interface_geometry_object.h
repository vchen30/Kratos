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

#if !defined(KRATOS_INTERFACE_GEOMETRY_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_GEOMETRY_OBJECT_INCLUDED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_object.h"


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

/// GeometricalObject-based objects (Element or Condition) on the Interface for Searching
/** This class Is the "wrapper" for Elements/Conditions on the interface. It uses the fact that both
* Elements and Conditions are deriving from "GeometricalObject". The search is caarried out using the
* center of the geometry.
* It saves a pointer to the original geometry, not to the Condition/Element itself. This is e.g. why the Id is not accessible.
* It selects the best result by the closest projection distance of the successful projections.
* In case no projection is successful, it uses an approximation (closest node of the geometry with the
* smallest center distance to the point for which a neighbor is to be found)
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceGeometryObject : public InterfaceObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceGeometryObject
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceGeometryObject);

    ///@}
    ///@name Life Cycle
    ///@{

    // A default constructor necessary for serialization
    InterfaceGeometryObject() : InterfaceObject()
    {
    }

    InterfaceGeometryObject(GeometricalObject& rGeometricalObject) :
        InterfaceObject(rGeometricalObject.GetGeometry().Center()), mpGeometricalObject(&rGeometricalObject)
    { }

    /// Destructor.
    virtual ~InterfaceGeometryObject() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    GeometricalObject* pGetBaseGeometricalObject() override
    {
        return mpGeometricalObject;
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceGeometryObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceGeometryObject";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    GeometricalObject* mpGeometricalObject;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }

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
    InterfaceGeometryObject& operator=(InterfaceGeometryObject const& rOther);

    //   /// Copy constructor.
    //   InterfaceGeometryObject(InterfaceGeometryObject const& rOther){}


    ///@}

}; // Class InterfaceGeometryObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceGeometryObject& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceGeometryObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_GEOMETRY_OBJECT_INCLUDED_H_INCLUDED  defined
