//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Feb 10, 2017$
//   Revision:            $Revision: 1.0 $
//
// 


// System includes


// External includes


// Project includes
#include "finite_cell_structural_application.h"
#include "finite_cell_application/custom_geometries/finite_cell_geometry.h"

namespace Kratos
{
    
    KratosFiniteCellStructuralApplication::KratosFiniteCellStructuralApplication()
    :
      mKinematicLinearFiniteCell2D3N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D3 <Node<3> > >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D6N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D6 <Node<3> > >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D9N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D9 <Node<3> > >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D10N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D10 <Node<3> > >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D20N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D20 <Node<3> > >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D27N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D27 <Node<3> > >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mLineForceWithFunction2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForceWithFunction2D3N( 0, Element::GeometryType::Pointer( new Line2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mLineForceWithFunction3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForceWithFunction3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    {}

    void KratosFiniteCellStructuralApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosFiniteCellStructuralApplication... " << std::endl;

        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell2D3N", mKinematicLinearFiniteCell2D3N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell2D6N", mKinematicLinearFiniteCell2D6N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell2D4N", mKinematicLinearFiniteCell2D4N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell2D8N", mKinematicLinearFiniteCell2D8N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell2D9N", mKinematicLinearFiniteCell2D9N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell3D4N", mKinematicLinearFiniteCell3D4N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell3D10N", mKinematicLinearFiniteCell3D10N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell3D8N", mKinematicLinearFiniteCell3D8N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell3D20N", mKinematicLinearFiniteCell3D20N )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearFiniteCell3D27N", mKinematicLinearFiniteCell3D27N )

        KRATOS_REGISTER_CONDITION( "LineForceWithFunction2D2N", mLineForceWithFunction2D2N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction2D3N", mLineForceWithFunction2D3N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction3D2N", mLineForceWithFunction3D2N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction3D3N", mLineForceWithFunction3D3N )

    }

} // namespace Kratos

