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
#include "geometries/point_3d.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "finite_cell_application/custom_geometries/finite_cell_geometry.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE( double, FORCE_MAGNITUDE )
    KRATOS_CREATE_VARIABLE( double, PHYSICAL_STRESS_OFFSET_PARAMETER )
    KRATOS_CREATE_VARIABLE( double, STRESS_STABILIZATION )
    KRATOS_CREATE_VARIABLE( int, GHOST_PENALTY_STABILIZATION_ORDER )
    KRATOS_CREATE_VARIABLE( double, GHOST_PENALTY_STABILIZATION_FACTOR )

    KratosFiniteCellStructuralApplication::KratosFiniteCellStructuralApplication()
    : KratosApplication()
    , mKinematicLinearFiniteCell2D3N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D3 <Node<3> > >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D6N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D6 <Node<3> > >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell2D9N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D9 <Node<3> > >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D10N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D10 <Node<3> > >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D20N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D20 <Node<3> > >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mKinematicLinearFiniteCell3D27N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D27 <Node<3> > >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell2D3N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D3 <Node<3> > >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell2D6N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D6 <Node<3> > >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell2D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell2D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell2D9N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D9 <Node<3> > >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell3D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell3D10N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D10 <Node<3> > >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell3D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell3D20N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D20 <Node<3> > >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mTotalLagrangianFiniteCell3D27N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D27 <Node<3> > >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell2D3N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D3 <Node<3> > >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell2D6N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D6 <Node<3> > >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell2D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell2D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell2D9N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D9 <Node<3> > >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell3D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell3D10N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D10 <Node<3> > >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell3D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell3D20N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D20 <Node<3> > >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mExtrapolatedKinematicLinearFiniteCell3D27N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D27 <Node<3> > >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell2D3N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D3 <Node<3> > >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell2D6N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Triangle2D6 <Node<3> > >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell2D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell2D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell2D9N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Quadrilateral2D9 <Node<3> > >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell3D4N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D4 <Node<3> > >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell3D10N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Tetrahedra3D10 <Node<3> > >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell3D8N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D8 <Node<3> > >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell3D20N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D20 <Node<3> > >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mExtrapolatedConstantStressKinematicLinearFiniteCell3D27N( 0, Element::GeometryType::Pointer( new FiniteCellGeometry< Hexahedra3D27 <Node<3> > >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mLineForceWithFunction2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForceWithFunction2D3N( 0, Element::GeometryType::Pointer( new Line2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mLineForceWithFunction3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForceWithFunction3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFaceForceWithFunction3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFaceForceWithFunction3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mFaceForceWithFunction3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mFaceForceWithFunction3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mFaceForceWithFunction3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mImmersedPointForce3D( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    {}

    void KratosFiniteCellStructuralApplication::Register()
    {
        std::cout << "Initializing KratosFiniteCellStructuralApplication... " << std::endl;

        KRATOS_REGISTER_VARIABLE( FORCE_MAGNITUDE )
        KRATOS_REGISTER_VARIABLE( PHYSICAL_STRESS_OFFSET_PARAMETER )
        KRATOS_REGISTER_VARIABLE( STRESS_STABILIZATION )
        KRATOS_REGISTER_VARIABLE( GHOST_PENALTY_STABILIZATION_ORDER )
        KRATOS_REGISTER_VARIABLE( GHOST_PENALTY_STABILIZATION_FACTOR )

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

        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell2D3N", mTotalLagrangianFiniteCell2D3N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell2D6N", mTotalLagrangianFiniteCell2D6N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell2D4N", mTotalLagrangianFiniteCell2D4N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell2D8N", mTotalLagrangianFiniteCell2D8N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell2D9N", mTotalLagrangianFiniteCell2D9N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell3D4N", mTotalLagrangianFiniteCell3D4N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell3D10N", mTotalLagrangianFiniteCell3D10N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell3D8N", mTotalLagrangianFiniteCell3D8N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell3D20N", mTotalLagrangianFiniteCell3D20N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFiniteCell3D27N", mTotalLagrangianFiniteCell3D27N )


        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell2D3N", mExtrapolatedKinematicLinearFiniteCell2D3N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell2D6N", mExtrapolatedKinematicLinearFiniteCell2D6N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell2D4N", mExtrapolatedKinematicLinearFiniteCell2D4N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell2D8N", mExtrapolatedKinematicLinearFiniteCell2D8N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell2D9N", mExtrapolatedKinematicLinearFiniteCell2D9N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell3D4N", mExtrapolatedKinematicLinearFiniteCell3D4N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell3D10N", mExtrapolatedKinematicLinearFiniteCell3D10N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell3D8N", mExtrapolatedKinematicLinearFiniteCell3D8N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell3D20N", mExtrapolatedKinematicLinearFiniteCell3D20N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedKinematicLinearFiniteCell3D27N", mExtrapolatedKinematicLinearFiniteCell3D27N )

        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell2D3N", mExtrapolatedConstantStressKinematicLinearFiniteCell2D3N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell2D6N", mExtrapolatedConstantStressKinematicLinearFiniteCell2D6N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell2D4N", mExtrapolatedConstantStressKinematicLinearFiniteCell2D4N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell2D8N", mExtrapolatedConstantStressKinematicLinearFiniteCell2D8N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell2D9N", mExtrapolatedConstantStressKinematicLinearFiniteCell2D9N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell3D4N", mExtrapolatedConstantStressKinematicLinearFiniteCell3D4N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell3D10N", mExtrapolatedConstantStressKinematicLinearFiniteCell3D10N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell3D8N", mExtrapolatedConstantStressKinematicLinearFiniteCell3D8N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell3D20N", mExtrapolatedConstantStressKinematicLinearFiniteCell3D20N )
        KRATOS_REGISTER_ELEMENT( "ExtrapolatedConstantStressKinematicLinearFiniteCell3D27N", mExtrapolatedConstantStressKinematicLinearFiniteCell3D27N )

        KRATOS_REGISTER_CONDITION( "LineForceWithFunction2D2N", mLineForceWithFunction2D2N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction2D3N", mLineForceWithFunction2D3N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction3D2N", mLineForceWithFunction3D2N )
        KRATOS_REGISTER_CONDITION( "LineForceWithFunction3D3N", mLineForceWithFunction3D3N )

        KRATOS_REGISTER_CONDITION( "FaceForceWithFunction3D3N", mFaceForceWithFunction3D3N )
        KRATOS_REGISTER_CONDITION( "FaceForceWithFunction3D6N", mFaceForceWithFunction3D6N )
        KRATOS_REGISTER_CONDITION( "FaceForceWithFunction3D4N", mFaceForceWithFunction3D4N )
        KRATOS_REGISTER_CONDITION( "FaceForceWithFunction3D8N", mFaceForceWithFunction3D8N )
        KRATOS_REGISTER_CONDITION( "FaceForceWithFunction3D9N", mFaceForceWithFunction3D9N )

        KRATOS_REGISTER_CONDITION( "ImmersedPointForce3D", mImmersedPointForce3D )
    }

} // namespace Kratos
