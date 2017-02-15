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
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
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
    {}

    void KratosFiniteCellStructuralApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosFiniteCellStructuralApplication... " << std::endl;
    }

} // namespace Kratos

