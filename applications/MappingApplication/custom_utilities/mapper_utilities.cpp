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


// System includes


// External includes


// Project includes
#include "mapper_utilities.h"


namespace Kratos
{
namespace MapperUtilities
{

    typedef std::unordered_map<std::string, InterfaceSearchObject::Pointer> InterfaceObjectMap; // TODO redefinition necessary?

    void MapperUtilities::RegisterInterfaceSearchObject(const std::string& InterfaceSearchObjectName,
                                                        InterfaceSearchObject::Pointer pInterfaceSearchObjectPrototype)
    {
        GetRegisteredInterfaceSearchObjectsList().insert(make_pair(InterfaceSearchObjectName, pInterfaceSearchObjectPrototype));
    }

    InterfaceObjectMap& MapperUtilities::GetRegisteredInterfaceSearchObjectsList()
    {
        static InterfaceObjectMap registered_interface_search_objects;

        return registered_interface_search_objects;
    }

    InterfaceSearchObject::Pointer MapperUtilities::GetInterfaceSearchObject(const std::string& rInterfaceSearchObjectName)
    {
        const auto& r_interface_search_object_list = MapperFactory::GetRegisteredMappersList();

        if (r_interface_search_object_list.find(rInterfaceSearchObjectName) != r_interface_search_object_list.end())
        {
            return r_interface_search_object_list.at(rInterfaceSearchObjectName)->Clone();
        }
        else
        {
            // TODO throw proper error message

            // std::stringstream err_msg;
            // err_msg << "The requested Mapper \"" << mapper_name <<"\" is not registered! The following mappers are registered:" << std::endl;
            // for (auto const& registered_mapper : mapper_list)
            // {
            //     err_msg << registered_mapper.first << ", ";
            // }
            // KRATOS_ERROR << err_msg.str() << std::endl;
        }
    }

}  // namespace MapperUtilities.
}  // namespace Kratos.