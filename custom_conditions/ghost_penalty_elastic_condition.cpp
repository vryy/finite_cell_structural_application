//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/ghost_penalty_elastic_condition.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "finite_cell_application/custom_utilities/ghost_penalty_utility.h"
#include "structural_application/structural_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
GhostPenaltyElasticCondition::GhostPenaltyElasticCondition()
{
}

GhostPenaltyElasticCondition::GhostPenaltyElasticCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: GhostPenaltyCondition( NewId, pGeometry )
{
}

GhostPenaltyElasticCondition::GhostPenaltyElasticCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: GhostPenaltyCondition( NewId, pGeometry, pProperties )
{
}

GhostPenaltyElasticCondition::GhostPenaltyElasticCondition( IndexType NewId, GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement, Element::Pointer pMasterElement, PropertiesType::Pointer pProperties )
: GhostPenaltyCondition(NewId, pGeometry, pSlaveElement, pMasterElement, pProperties)
{
}

/**
 * Destructor. Never to be called manually
 */
GhostPenaltyElasticCondition::~GhostPenaltyElasticCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer GhostPenaltyElasticCondition::Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties)
{
    return Condition::Pointer( new GhostPenaltyElasticCondition(NewId, pGeometry, pSlaveElement, pMasterElement, pProperties) );
}

void GhostPenaltyElasticCondition::Initialize()
{
    KRATOS_TRY

    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void GhostPenaltyElasticCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector,
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag,
                  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void GhostPenaltyElasticCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************    /
/**
 * This function calculates all system contributions due to the ghost penalty problem
 */
void GhostPenaltyElasticCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    unsigned int Dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = GetGeometry().size() * Dim;

    // resize the LHS accordingly
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    // compute the shape function gradient on the edge on the normal direction
    Matrix dN1dn(GetGeometry().size(), integration_points.size());
    GhostPenaltyUtility::ComputeShapeFunctionNormalGradient(dN1dn, pSlave()->GetGeometry(), GetGeometry(), integration_points);
    KRATOS_WATCH(dN1dn)

    Matrix dN2dn(GetGeometry().size(), integration_points.size());
    GhostPenaltyUtility::ComputeShapeFunctionNormalGradient(dN2dn, pMaster()->GetGeometry(), GetGeometry(), integration_points);
    KRATOS_WATCH(dN2dn)

    //initializing the Jacobian in the reference configuration
    Matrix DeltaPosition(GetGeometry().size(), 3);
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );
    double DetJ;

    // compute the length of the geometry
    double Length = 0.0;
    for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        DetJ = MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber])));
        Length += DetJ * integration_points[PointNumber].Weight();
    }
    KRATOS_WATCH(Length)

    const double& Penalty = this->GetValue(INITIAL_PENALTY);
    double Aux = Penalty * pow(Length, 2);

    if ( CalculateResidualVectorFlag == true )
    {
        Vector jump(Dim);

        // quadrature loop
        for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            // compute the jump term
            noalias(jump) = ZeroVector(Dim);
            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                jump(0) += dN1dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_X);
                jump(1) += dN1dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_Y);
                if (Dim == 3)
                    jump(2) += dN1dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_Z);
            }

            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                jump(0) -= dN2dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_X);
                jump(1) -= dN2dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_Y);
                if (Dim == 3)
                    jump(2) -= dN2dn(prim, PointNumber) * GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT_Z);
            }
            KRATOS_WATCH(jump)

            // compute determinant
            DetJ = MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber])));
//            KRATOS_WATCH(DetJ)

            // compute the right hand side vector
            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    rRightHandSideVector(prim*Dim + dim) += Aux * dN1dn(prim, PointNumber) * jump(dim) * integration_points[PointNumber].Weight() * DetJ;
                }
            }

            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    rRightHandSideVector(prim*Dim + dim) -= Aux * dN2dn(prim, PointNumber) * jump(dim) * integration_points[PointNumber].Weight() * DetJ;
                }
            }
        }

        KRATOS_WATCH(rRightHandSideVector)
    }

    if ( CalculateStiffnessMatrixFlag == true )
    {
        Vector djump(GetGeometry().size());

        // quadrature loop
        for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            // compute derivatives of the jumps
            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                djump(prim) = dN1dn(prim, PointNumber) - dN2dn(prim, PointNumber);
            }
            KRATOS_WATCH(djump)

            // compute determinant
            DetJ = MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber])));

            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    for (std::size_t sec = 0; sec < GetGeometry().size(); ++sec)
                    {
                        rLeftHandSideMatrix(prim*Dim + dim, sec*Dim + dim) -= Aux * dN1dn(prim, PointNumber) * djump(sec) * integration_points[PointNumber].Weight() * DetJ;
                    }
                }
            }

            // compute the left hand side matrix
            for (std::size_t prim = 0; prim < GetGeometry().size(); ++prim)
            {
                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    for (std::size_t sec = 0; sec < GetGeometry().size(); ++sec)
                    {
                        rLeftHandSideMatrix(prim*Dim + dim, sec*Dim + dim) += Aux * dN2dn(prim, PointNumber) * djump(sec) * integration_points[PointNumber].Weight() * DetJ;
                    }
                }
            }
        }

        KRATOS_WATCH(rLeftHandSideMatrix)
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void GhostPenaltyElasticCondition::EquationIdVector( EquationIdVectorType& rResult,
                                      ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = ( GetGeometry().WorkingSpaceDimension() );
    unsigned int mat_size = GetGeometry().size() * dim;

    if ( rResult.size() != mat_size )
        rResult.resize( mat_size, false );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        int index = i * dim;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if(dim == 3)
            rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list for the current partners.
 * All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per Node.
 * All DOF are given Master first, Slave second
 */
//************************************************************************************
//************************************************************************************
void GhostPenaltyElasticCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    ConditionalDofList.resize( 0 );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if(dim == 3)
            ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

} // Namespace Kratos

