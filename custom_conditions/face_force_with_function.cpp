// see finite_cell_application/LICENSE.txt
/* **************************************************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: Aug 26, 2019 $
*   Revision:            $Revision: 1.0 $
*
* ***************************************************************************************/


// System includes


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/legacy_structural_app_vars.h"
#include "custom_conditions/face_force_with_function.h"
#include "brep_application/custom_algebra/function/function.h"
#include "brep_application/brep_application.h"

namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
FaceForceWithFunction::FaceForceWithFunction()
{
}

// Constructor
FaceForceWithFunction::FaceForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
}

// Constructor
FaceForceWithFunction::FaceForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer FaceForceWithFunction::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FaceForceWithFunction( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer FaceForceWithFunction::Create( IndexType NewId,
                                        GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FaceForceWithFunction( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
FaceForceWithFunction::~FaceForceWithFunction()
{
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::EquationIdVector( EquationIdVectorType& rResult,
                                    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    DofsVectorType ConditionalDofList;
    GetDofList(ConditionalDofList, rCurrentProcessInfo);

    if (rResult.size() != ConditionalDofList.size())
        rResult.resize(ConditionalDofList.size(), false);

    for(unsigned int i = 0; i < ConditionalDofList.size(); ++i)
    {
        rResult[i] = ConditionalDofList[i]->EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::GetDofList( DofsVectorType& ElementalDofList,
                              const ProcessInfo& rCurrentProcessInfo ) const
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if(dim == 3)
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    //resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    GeometryData::IntegrationMethod ThisIntegrationMethod;
    if(GetProperties().Has( INTEGRATION_ORDER ) == true)
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Does not support for more integration points", *this)
    }
    else
        ThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(ThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(ThisIntegrationMethod);

    // DN_DeContainer is the array of shape function gradients at each integration points
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod);

    // Ncontainer is the array of shape function values at each integration points
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(ThisIntegrationMethod);

//    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
//        KRATOS_WATCH(integration_points[PointNumber])
//    KRATOS_WATCH(Ncontainer)

    //loop over integration points
    Vector Load( dim );
    Vector t( dim );
    PointType GlobalCoords;

    FunctionR3Rn::Pointer pLoadFuntion;

    #pragma omp critical
    {
        boost::python::object pyObject = GetProperties()[LOAD_FUNCTION];
        try
        {
            pLoadFuntion = boost::python::extract<FunctionR3Rn::Pointer>(pyObject);
        }
        catch (boost::python::error_already_set&)
        {
            std::cout << "Extracting load function has some problem" << std::endl;
            PyErr_Print();
        }
    }

    if(pLoadFuntion == NULL)
        KRATOS_THROW_ERROR(std::logic_error, "LOAD_FUNCTION is not provided for FaceForceWithFunction", Id())

    Vector t1(3), t2(3), v3(3);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        noalias(GlobalCoords) = ZeroVector(3);
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
            noalias(GlobalCoords) += Ncontainer(PointNumber, i) * GetGeometry()[i].GetInitialPosition();

        noalias( Load ) = pLoadFuntion->GetValue(GlobalCoords);
//        KRATOS_WATCH(GlobalCoords)
//        KRATOS_WATCH(Load)

        double IntegrationWeight = integration_points[PointNumber].Weight();
//        KRATOS_WATCH(IntegrationWeight)

        noalias(t1) = ZeroVector( 3 );//first tangential vector
        noalias(t2) = ZeroVector( 3 );//second tangential vector
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t1[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
            t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
            t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 0 );
            t2[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 1 );
            t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 1 );
            t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 1 );
        }

        //calculating normal
        v3[0] = t1[1] * t2[2] - t1[2] * t2[1];

        v3[1] = t1[2] * t2[0] - t1[0] * t2[2];

        v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

        // calculating area
        double dA = norm_2(v3);
//        KRATOS_WATCH(v3)
//        KRATOS_WATCH(dA)

        // RIGHT HAND SIDE VECTOR
        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
            for ( unsigned int i = 0; i < dim; ++i )
                rRightHandSideVector( prim * dim + i ) +=
                    Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dA;
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //clean the internal data of the geometry
    GetGeometry().Clean();
    #endif

//    KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    rLeftHandSideMatrix = ZeroMatrix( dim * number_of_nodes, dim * number_of_nodes );
    CalculateRightHandSide( rRightHandSideVector, rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void FaceForceWithFunction::CalculateDampingMatrix( MatrixType& rDampingMatrix,
                              const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rDampingMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int FaceForceWithFunction::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
