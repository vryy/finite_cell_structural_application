//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/ghost_penalty_stress_condition.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "finite_cell_application/custom_utilities/ghost_penalty_utility.h"
#include "structural_application/structural_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
GhostPenaltyStressCondition::GhostPenaltyStressCondition()
{
    mIsInitialized = false;
}

GhostPenaltyStressCondition::GhostPenaltyStressCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: GhostPenaltyCondition( NewId, pGeometry )
{
    mIsInitialized = false;
}

GhostPenaltyStressCondition::GhostPenaltyStressCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: GhostPenaltyCondition( NewId, pGeometry, pProperties )
{
    mIsInitialized = false;
}

GhostPenaltyStressCondition::GhostPenaltyStressCondition( IndexType NewId, GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement, Element::Pointer pMasterElement, PropertiesType::Pointer pProperties )
: GhostPenaltyCondition(NewId, pGeometry, pSlaveElement, pMasterElement, pProperties)
{
    mIsInitialized = false;
}

/**
 * Destructor. Never to be called manually
 */
GhostPenaltyStressCondition::~GhostPenaltyStressCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer GhostPenaltyStressCondition::Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties)
{
    return Condition::Pointer( new GhostPenaltyStressCondition(NewId, pGeometry, pSlaveElement, pMasterElement, pProperties) );
}

//************************************************************************************
//************************************************************************************

void GhostPenaltyStressCondition::Initialize()
{
    GhostPenaltyCondition::Initialize();

    //number of integration points used, mThisIntegrationMethod refers to the
    //integration method defined in the constructor
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Initialization of the constitutive law vector and
    // declaration, definition and initialization of the material
    // lwas at each integration point
    if ( mSlaveConstitutiveLaws.size() != integration_points.size() )
    {
        mSlaveConstitutiveLaws.resize( integration_points.size() );

        #ifdef ENABLE_BEZIER_GEOMETRY
        pSlave()->GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        for ( unsigned int i = 0; i < mSlaveConstitutiveLaws.size(); ++i )
        {
            // REMARKS!!! It is noted that the constitutive law on the interface may be different to the constitutive law of the bulk, if their Properties are different.
            mSlaveConstitutiveLaws[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mSlaveConstitutiveLaws[i]->SetValue( PARENT_ELEMENT_ID, pSlave()->Id(), *(ProcessInfo*)0);
            mSlaveConstitutiveLaws[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
            mSlaveConstitutiveLaws[i]->InitializeMaterial( GetProperties(), pSlave()->GetGeometry(), row( pSlave()->GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

            //check constitutive law
            mSlaveConstitutiveLaws[i]->Check( GetProperties(), pSlave()->GetGeometry(), *(ProcessInfo*)0 );
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        pSlave()->GetGeometry().Clean();
        #endif
    }

    if ( mMasterConstitutiveLaws.size() != integration_points.size() )
    {
        mMasterConstitutiveLaws.resize( integration_points.size() );

        #ifdef ENABLE_BEZIER_GEOMETRY
        pMaster()->GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        for ( unsigned int i = 0; i < mMasterConstitutiveLaws.size(); ++i )
        {
            // REMARKS!!! It is noted that the constitutive law on the interface may be different to the constitutive law of the bulk, if their Properties are different.
            mMasterConstitutiveLaws[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mMasterConstitutiveLaws[i]->SetValue( PARENT_ELEMENT_ID, pMaster()->Id(), *(ProcessInfo*)0);
            mMasterConstitutiveLaws[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
            mMasterConstitutiveLaws[i]->InitializeMaterial( GetProperties(), pMaster()->GetGeometry(), row( pMaster()->GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

            //check constitutive law
            mMasterConstitutiveLaws[i]->Check( GetProperties(), pMaster()->GetGeometry(), *(ProcessInfo*)0 );
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        pMaster()->GetGeometry().Clean();
        #endif
    }
}

//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void GhostPenaltyStressCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void GhostPenaltyStressCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void GhostPenaltyStressCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    unsigned int Dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = GetGeometry().size() * Dim;
    unsigned int strain_size = Dim * (Dim + 1) / 2;

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

    if ( CalculateResidualVectorFlag || CalculateStiffnessMatrixFlag )
    {
        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        // reading integration points and local gradients on the edge/face
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );
        double DetJ;

        // compute the length/area of the geometry
        double DomainSize = 0.0;
        for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            DetJ = MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber])));
            DomainSize += DetJ * integration_points[PointNumber].Weight();
        }
//        KRATOS_WATCH(DomainSize)

        // compute the penalty coefficients
        const double& Penalty = GetProperties()[INITIAL_PENALTY];
        double Aux;
        if (Dim == 2)
            Aux = Penalty * pow(DomainSize, 2);
        else if (Dim == 3)
            Aux = Penalty * DomainSize;

        Matrix B1_Operator( strain_size, pSlave()->GetGeometry().size()*Dim );
        Matrix TanC_1( strain_size, strain_size );
        Vector StrainVector_1( strain_size );
        Vector StressVector_1( strain_size );
        Matrix CurrentDisp_1( pSlave()->GetGeometry().size(), Dim );
        Vector Nvalues_1( pSlave()->GetGeometry().size() );
        Matrix DN1_DX( pSlave()->GetGeometry().size(), Dim );
        IntegrationPointType integration_point_1;
        Vector RightHandSideVector_1( pSlave()->GetGeometry().size()*Dim );
        Matrix LeftHandSideMatrix_1( pSlave()->GetGeometry().size()*Dim, pSlave()->GetGeometry().size()*Dim );

        noalias( RightHandSideVector_1 ) = ZeroVector( pSlave()->GetGeometry().size()*Dim );
        noalias( LeftHandSideMatrix_1 ) = ZeroMatrix( pSlave()->GetGeometry().size()*Dim, pSlave()->GetGeometry().size()*Dim );

        Matrix B2_Operator( strain_size, pMaster()->GetGeometry().size()*Dim );
        Matrix TanC_2( strain_size, strain_size );
        Vector StrainVector_2( strain_size );
        Vector StressVector_2( strain_size );
        Matrix CurrentDisp_2( pMaster()->GetGeometry().size(), Dim );
        Vector Nvalues_2( pMaster()->GetGeometry().size() );
        Matrix DN2_DX( pMaster()->GetGeometry().size(), Dim );
        IntegrationPointType integration_point_2;
        Vector RightHandSideVector_2( pMaster()->GetGeometry().size()*Dim );
        Matrix LeftHandSideMatrix_2( pMaster()->GetGeometry().size()*Dim, pMaster()->GetGeometry().size()*Dim );

        noalias( RightHandSideVector_2 ) = ZeroVector( pMaster()->GetGeometry().size()*Dim );
        noalias( LeftHandSideMatrix_2 ) = ZeroMatrix( pMaster()->GetGeometry().size()*Dim, pMaster()->GetGeometry().size()*Dim );

        // extract current displacements
        for (unsigned int node = 0; node < pSlave()->GetGeometry().size(); ++node)
            noalias(row(CurrentDisp_1, node)) =
                pSlave()->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

        for (unsigned int node = 0; node < pMaster()->GetGeometry().size(); ++node)
            noalias(row(CurrentDisp_2, node)) =
                pMaster()->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

        int side_1 = GhostPenaltyUtility::FindSide(pSlave()->GetGeometry(), GetGeometry());
        int side_2 = GhostPenaltyUtility::FindSide(pMaster()->GetGeometry(), GetGeometry());
//        KRATOS_WATCH(side_1)
//        KRATOS_WATCH(side_2)
        // Remark: I don't want to check size because when setting up the ghost penalty, the size is automatically correct.

        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index_1;
        GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index_1,
                GetGeometry(), pSlave()->GetGeometry(), GhostPenaltyUtility::Faces(pSlave()->GetGeometry(), side_1));

        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index_2;
        GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index_2,
                GetGeometry(), pMaster()->GetGeometry(), GhostPenaltyUtility::Faces(pMaster()->GetGeometry(), side_2));

        // quadrature loop
        for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            // compute the integration points of the edge on the slave and the master
            noalias( integration_point_1 ) = GhostPenaltyUtility::ComputeIntegrationPoint(pSlave()->GetGeometry(), side_1, integration_points[PointNumber]);
            noalias( integration_point_2 ) = GhostPenaltyUtility::ComputeIntegrationPoint(pMaster()->GetGeometry(), side_2, integration_points[PointNumber]);

            /* compute the jump term */

            // compute the shape function values
            Nvalues_1 = pSlave()->GetGeometry().ShapeFunctionsValues(Nvalues_1, integration_point_1);
            Nvalues_2 = pMaster()->GetGeometry().ShapeFunctionsValues(Nvalues_2, integration_point_2);

            // compute the shape function gradient
            GhostPenalty_Helper::ComputeShapeFunctionGradient(DN1_DX, pSlave()->GetGeometry(), integration_point_1);
            GhostPenalty_Helper::ComputeShapeFunctionGradient(DN2_DX, pMaster()->GetGeometry(), integration_point_2);

            // compute the B Operator
            CalculateBoperator(B1_Operator, DN1_DX);
            CalculateBoperator(B2_Operator, DN2_DX);

            // compute the strain
            CalculateStrain(B1_Operator, CurrentDisp_1, StrainVector_1);
            CalculateStrain(B2_Operator, CurrentDisp_2, StrainVector_2);

            // compute the stress
            mSlaveConstitutiveLaws[PointNumber]->CalculateMaterialResponse(
                StrainVector_1,
                ZeroMatrix( 1 ),
                StressVector_1,
                TanC_1,
                rCurrentProcessInfo,
                GetProperties(),
                pSlave()->GetGeometry(),
                Nvalues_1,
                true,
                ( int )CalculateStiffnessMatrixFlag,
                true
            );

            mMasterConstitutiveLaws[PointNumber]->CalculateMaterialResponse(
                StrainVector_2,
                ZeroMatrix( 1 ),
                StressVector_2,
                TanC_2,
                rCurrentProcessInfo,
                GetProperties(),
                pMaster()->GetGeometry(),
                Nvalues_2,
                true,
                ( int )CalculateStiffnessMatrixFlag,
                true
            );

            if ( CalculateResidualVectorFlag == true )
            {
                noalias( RightHandSideVector_1 ) += prod( trans(B1_Operator), StressVector_1 );
                noalias( RightHandSideVector_2 ) += prod( trans(B2_Operator), StressVector_2 );
            }

            if ( CalculateStiffnessMatrixFlag == true )
            {
                noalias( LeftHandSideMatrix_1 ) += prod( trans( B1_Operator ), Matrix( prod( TanC_1, B1_Operator ) ) );
                noalias( LeftHandSideMatrix_2 ) += prod( trans( B2_Operator ), Matrix( prod( TanC_2, B2_Operator ) ) );
            }
        }

        if ( CalculateResidualVectorFlag == true )
        {
            // assemble the jump term to the nodes on edge/face
            for (std::map<std::size_t, std::size_t>::iterator it = map_edge_node_index_to_element_node_index_1.begin();
                    it != map_edge_node_index_to_element_node_index_1.end(); ++it)
            {
                const std::size_t& edge_node_index = it->first;
                const std::size_t& element_node_index = it->second;

                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    rRightHandSideVector( edge_node_index*Dim + dim ) += RightHandSideVector_1( element_node_index*Dim + dim );
                }
            }

            for (std::map<std::size_t, std::size_t>::iterator it = map_edge_node_index_to_element_node_index_2.begin();
                    it != map_edge_node_index_to_element_node_index_2.end(); ++it)
            {
                const std::size_t& edge_node_index = it->first;
                const std::size_t& element_node_index = it->second;

                for (std::size_t dim = 0; dim < Dim; ++dim)
                {
                    rRightHandSideVector( edge_node_index*Dim + dim ) -= RightHandSideVector_2( element_node_index*Dim + dim );
                }
            }
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            for (std::map<std::size_t, std::size_t>::iterator it1 = map_edge_node_index_to_element_node_index_1.begin();
                    it1 != map_edge_node_index_to_element_node_index_1.end(); ++it1)
            {
                const std::size_t& edge_node_index1 = it1->first;
                const std::size_t& element_node_index1 = it1->second;

                for (std::size_t dim1 = 0; dim1 < Dim; ++dim1)
                {
                    for (std::map<std::size_t, std::size_t>::iterator it2 = map_edge_node_index_to_element_node_index_1.begin();
                            it2 != map_edge_node_index_to_element_node_index_1.end(); ++it2)
                    {
                        const std::size_t& edge_node_index2 = it2->first;
                        const std::size_t& element_node_index2 = it2->second;

                        for (std::size_t dim2 = 0; dim2 < Dim; ++dim2)
                        {
                            rLeftHandSideMatrix( edge_node_index1*Dim + dim1, edge_node_index2*Dim + dim2 )
                                -= LeftHandSideMatrix_1( element_node_index1*Dim + dim1, element_node_index2*Dim + dim2 );
                        }
                    }
                }
            }

            for (std::map<std::size_t, std::size_t>::iterator it1 = map_edge_node_index_to_element_node_index_2.begin();
                    it1 != map_edge_node_index_to_element_node_index_2.end(); ++it1)
            {
                const std::size_t& edge_node_index1 = it1->first;
                const std::size_t& element_node_index1 = it1->second;

                for (std::size_t dim1 = 0; dim1 < Dim; ++dim1)
                {
                    for (std::map<std::size_t, std::size_t>::iterator it2 = map_edge_node_index_to_element_node_index_2.begin();
                            it2 != map_edge_node_index_to_element_node_index_2.end(); ++it2)
                    {
                        const std::size_t& edge_node_index2 = it2->first;
                        const std::size_t& element_node_index2 = it2->second;

                        for (std::size_t dim2 = 0; dim2 < Dim; ++dim2)
                        {
                            rLeftHandSideMatrix( edge_node_index1*Dim + dim1, edge_node_index2*Dim + dim2 )
                                += LeftHandSideMatrix_2( element_node_index1*Dim + dim1, element_node_index2*Dim + dim2 );
                        }
                    }
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Clean();
        #endif
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
void GhostPenaltyStressCondition::EquationIdVector( EquationIdVectorType& rResult,
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
void GhostPenaltyStressCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
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

//************************************************************************************
//************************************************************************************
void GhostPenaltyStressCondition::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX ) const
{
    KRATOS_TRY

    const unsigned int number_of_nodes = DN_DX.size1();
    unsigned int dim = DN_DX.size2();
    unsigned int strain_size = dim * (dim + 1) / 2;

    if ( B_Operator.size1() != strain_size || B_Operator.size2() != number_of_nodes * dim )
        B_Operator.resize( strain_size, number_of_nodes * dim, false );
    noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

    if ( dim == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*2 ) = DN_DX( i, 0 );
            B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
        }
    }
    else if ( dim == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*3 ) = DN_DX( i, 0 );
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
            B_Operator( 3, i*3 ) = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3 ) = DN_DX( i, 2 );
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void GhostPenaltyStressCondition::CalculateStrain( const Matrix& B,
        const Matrix& Displacements, Vector& StrainVector ) const
{
    KRATOS_TRY

    unsigned int Dim = Displacements.size2();
    unsigned int strain_size = StrainVector.size();
    noalias( StrainVector ) = ZeroVector( strain_size );

    for ( unsigned int node = 0; node < Displacements.size1(); ++node )
    {
        for ( unsigned int item = 0; item < strain_size; ++item )
            for ( unsigned int dim = 0; dim < Dim; ++dim )
                StrainVector[item] += B( item, Dim * node + dim ) * Displacements( node, dim );
    }

    KRATOS_CATCH( "" )
}

} // Namespace Kratos

