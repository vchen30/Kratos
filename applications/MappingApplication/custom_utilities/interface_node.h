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

#if !defined(KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_base_object.h"


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

/// Node on the Interface for Searching
/** This class Is the "wrapper" for nodes on the interface. It selects the best result by the closest distance to the
* point of which neighbor have to be found
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceNode : public InterfaceBaseObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceNode
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceNode);

    ///@}
    ///@name Life Cycle
    ///@{

    // A default constructor necessary for serialization
    InterfaceNode() : InterfaceBaseObject()
    {
    }

    InterfaceNode(Node<3>& rNode, const int EchoLevel) :
        InterfaceBaseObject(rNode.Coordinates()), mpNode(&rNode)
    {
    }

    /// Destructor.
    virtual ~InterfaceNode() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Node<3>* pGetBaseNode() override
    {
        return mpNode;
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
        buffer << "InterfaceNode" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceNode";
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

    Node<3>* mpNode;

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
    InterfaceNode& operator=(InterfaceNode const& rOther);

    //   /// Copy constructor.
    //   InterfaceNode(InterfaceNode const& rOther){}


    ///@}

}; // Class InterfaceNode

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceNode& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceNode& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED  defined
