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

#if !defined(KRATOS_INTERFACE_BASE_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_BASE_OBJECT_INCLUDED_H_INCLUDED

// System includes
#include<set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "mapper_utilities.h"
#include "../mapping_application_variables.h" // TODO remove "../"


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

/// Base Class for Searching Objects
/** This class provides a set of functions that is used for identifying the nearest object.
* It is needed such that the bin-search can be used with both nodes and elements/conditions
* The bin search is implemented to work with this kind of object
* It implements the function "EvaluateResult", which is used by the local search to determine which
* of the objects in teh vicinity of a point is the best search result. This function has to be
* implemented by all subclasses. It cannot be made pure virtual, because for remote searching,
* objects of this type have to be created
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceBaseObject : public Point
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::set<int> CandidateRankContainer;

    /// Pointer definition of InterfaceBaseObject
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceBaseObject);

    ///@}
    ///@name  Enum's
    ///@{

    enum PairingStatus
    {
        NoNeighbor = 0,
        Approximation = 1,
        NeighborFound = 2,
    };

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceBaseObject(const array_1d<double, 3>& rCoordinates) : Point(rCoordinates)
    { }

    /// Destructor.
    virtual ~InterfaceBaseObject() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual Node<3>* pGetBaseNode()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return nullptr;
    }

    virtual Geometry<Node<3>>* pGetBaseGeometry()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return nullptr;
    }

    /**
    Returns the custom search radius for the search
    */
    virtual double GetSearchRadius()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    /**
    This function is called from inside the bounding box
    The BB calls the Coordinates() function and checks whether the
    object is inside it or not. If yes it sets its rank to this object
    The BB is containing a partitioned ModelPart
    */
    void SetCandidatePartitionRank(int CandidateRank)
    {
        mCandidateRanks.insert(CandidateRank);
    }

    void ProcessSearchResult(InterfaceBaseObject::iterator)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    bool GetSuccessfulSearch()
    {
        return mSearchSuccessful;
    }

    int GetSourceRank()
    {
        return mSourceRank;
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
        std::stringstream buffer;
        buffer << "InterfaceBaseObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceBaseObject";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    int mEchoLevel = 0;

    CandidateRankContainer mCandidateRanks;
    bool mSearchSuccessful;
    int mSourceRank = 0;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
    }
    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
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
    InterfaceBaseObject& operator=(InterfaceBaseObject const& rOther);

    //   /// Copy constructor.
    //   InterfaceBaseObject(InterfaceBaseObject const& rOther){}


    ///@}

}; // Class InterfaceBaseObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceBaseObject& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceBaseObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_BASE_OBJECT_INCLUDED_H_INCLUDED  defined
