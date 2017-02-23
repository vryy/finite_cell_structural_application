// see finite_cell_application/LICENSE.txt
/* **************************************************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: Feb 20, 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***************************************************************************************/


// System includes


// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_conditions/line_force_with_function.h"
#include "finite_cell_application/finite_cell_application.h"
#include "structural_application/structural_application.h"

namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
LineForceWithFunction::LineForceWithFunction()
{
}

// Constructor
LineForceWithFunction::LineForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
}

// Constructor
LineForceWithFunction::LineForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineForceWithFunction::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LineForceWithFunction( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
LineForceWithFunction::~LineForceWithFunction()
{
}

//***********************************************************************************
//***********************************************************************************
void LineForceWithFunction::EquationIdVector( EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo )
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
void LineForceWithFunction::GetDofList( DofsVectorType& ElementalDofList,
                              ProcessInfo& rCurrentProcessInfo )
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
void LineForceWithFunction::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
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

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(ThisIntegrationMethod);

    // DN_DeContainer is the array of shape function gradients at each integration points
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod);

    // Ncontainer is the array of shape function values at each integration points
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(ThisIntegrationMethod);

    //loop over integration points
    Vector Load( dim );
    PointType GlobalCoords;

    typedef Function<PointType, Vector> FunctionR3RnType;

    FunctionR3RnType::Pointer pLoadFuntion;
    #pragma omp critical
    {
        boost::python::object pyObject = GetProperties()[LOAD_FUNCTION];
        pLoadFuntion = boost::python::extract<FunctionR3RnType::Pointer>(pyObject);

        if(pLoadFuntion == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "LOAD_FUNCTION is not provided for LineForceWithFunction", Id())

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            noalias(GlobalCoords) = ZeroVector(3);
            for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
                noalias(GlobalCoords) += Ncontainer(PointNumber, i) * GetGeometry()[i].GetInitialPosition();

            noalias( Load ) = pLoadFuntion->GetValue(GlobalCoords);
    //        KRATOS_WATCH(Load)

            double IntegrationWeight = integration_points[PointNumber].Weight();

            if(dim == 2) IntegrationWeight *= GetProperties()[THICKNESS];

            Vector t = ZeroVector( dim );//tangential vector
            for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
            {
                t[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
                t[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
                if(dim == 3)
                    t[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 0 );
            }

            // calculating length
            double dL = norm_2(t);
    //        KRATOS_WATCH(dL)

            // RIGHT HAND SIDE VECTOR
            for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
                for ( unsigned int i = 0; i < dim; ++i )
                    rRightHandSideVector( prim * dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dL;
        }

    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineForceWithFunction::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    rLeftHandSideMatrix = ZeroMatrix( dim * number_of_nodes, dim * number_of_nodes );
    CalculateRightHandSide( rRightHandSideVector, rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************
void LineForceWithFunction::CalculateMassMatrix( MatrixType& rMassMatrix,
                              ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineForceWithFunction::CalculateDampingMatrix( MatrixType& rDampingMatrix,
                              ProcessInfo& rCurrentProcessInfo )
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
int LineForceWithFunction::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
