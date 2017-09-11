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
#include "mapping_application.h"
#include "mapping_application_variables.h"


#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"


#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

namespace Kratos
{

KratosMappingApplication::KratosMappingApplication() :
    mNearestNeighborMapperCondition2D1N( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1 ) ) ) ),
    mNearestNeighborMapperCondition3D1N( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1 ) ) ) ),
    
    mNearestElementMapperCondition2D2NLine( 0, Element::GeometryType::Pointer(new Line2D2 <Node<3> >(Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mNearestElementMapperCondition3D2NLine( 0, Element::GeometryType::Pointer(new Line3D2 <Node<3> >(Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mNearestElementMapperCondition2D3NTri( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mNearestElementMapperCondition3D3NTri( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mNearestElementMapperCondition2D4NQuad( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mNearestElementMapperCondition3D4NQuad( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mNearestElementMapperCondition3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mNearestElementMapperCondition3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),

    mMortarMapperCondition2D2NLine( 0, Element::GeometryType::Pointer(new Line2D2 <Node<3> >(Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMortarMapperCondition3D2NLine( 0, Element::GeometryType::Pointer(new Line3D2 <Node<3> >(Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMortarMapperCondition2D3NTri( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMortarMapperCondition3D3NTri( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMortarMapperCondition2D4NQuad( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mMortarMapperCondition3D4NQuad( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mMortarMapperCondition3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mMortarMapperCondition3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
{ }

void KratosMappingApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner;
    banner << "    KRATOS ______  ___                      _____  "                          << std::endl;
    banner << "           ___   |/  /_____ ___________________(_)_____________ _  "          << std::endl;
    banner << "           __  /|_/ /_  __ `/__  __ \\__  __ \\_  /__  __ \\_  __ `/  "       << std::endl;
    banner << "           _  /  / / / /_/ /__  /_/ /_  /_/ /  / _  / / /  /_/ /  "           << std::endl;
    banner << "           /_/  /_/  \\__,_/ _  .___/_  .___//_/  /_/ /_/_\\__, /  "          << std::endl;
    banner << "                            /_/     /_/                 /____/ Application"   << std::endl;

    banner << "Initializing KratosMappingApplication... " << std::endl;

    int rank = 0;

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized)   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0) std::cout << banner.str();

    // Needed to exchange Information abt the found neighbors (i.e. only for debugging)
    KRATOS_REGISTER_VARIABLE( MAPPING_MATRIX_EQUATION_ID )
    KRATOS_REGISTER_VARIABLE( NEIGHBOR_RANK )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEIGHBOR_COORDINATES )
}
}  // namespace Kratos.