/*
see finite_cell_structural_application/LICENSE.txt
 */
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 24 Mar 2017 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
// #include "includes/define.h"
#include "extrapolated_kinematic_linear.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "finite_cell_application/finite_cell_application_variables.h"
#include "structural_application/structural_application_variables.h"
#include "finite_cell_structural_application/finite_cell_structural_application.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

//TODO: there is a potential bug at the CalculateRightHandSide, which is used for calculating the reaction. In principal, it should not tell the material to update themself, however, CalculateRightHandSide indirectly call CalculateMaterialResponse. THis should be fixed, by introducing another abstract layer to update the material in the input parameters for CalculateAll

namespace Kratos
{
    ExtrapolatedKinematicLinear::ExtrapolatedKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
    {
        mIsInitialized = false;
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    /**
     * A simple kinematic linear 3D element for the solution
     * of the momentum balance in structural mechanics.
     * This element is used for students training at the Ruhr University Bochum.
     * Therefore it may includes comments that are obvious for the
     * experienced user.
     */
    ExtrapolatedKinematicLinear::ExtrapolatedKinematicLinear( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
    {
        mIsInitialized = false;
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
    }

    Element::Pointer ExtrapolatedKinematicLinear::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new ExtrapolatedKinematicLinear( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer ExtrapolatedKinematicLinear::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new ExtrapolatedKinematicLinear( NewId, pGeom, pProperties ) );
    }

    ExtrapolatedKinematicLinear::~ExtrapolatedKinematicLinear()
    {
    }

    /**
     * Initialization of the element, called at the begin of each simulation.
     * Membervariables and the Material law are initialized here
     */
    void ExtrapolatedKinematicLinear::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY//EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

        //dimension of the problem
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        // integration rule
        if(this->Has( INTEGRATION_ORDER ))
        {
            if(this->GetValue(INTEGRATION_ORDER) == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "ExtrapolatedKinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
        }
        else if(GetProperties().Has( INTEGRATION_ORDER ))
        {
            if(GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "ExtrapolatedKinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
        }
        else
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

        if ( mIsInitialized )
        {
            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize( GetGeometry().size(), dim, false );

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                for ( unsigned int i = 0; i < dim; ++i )
                    mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
//            std::cout << "Element " << Id() << " mInitialDisp is reinitialized to " << mInitialDisp << std::endl;
            #endif

            return;
        }

        //number of integration points used, mThisIntegrationMethod refers to the
        //integration method defined in the constructor
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //calculating the domain size
        mTotalDomainInitialSize = 0.00;
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating the total domain size
            mTotalDomainInitialSize += MathUtils<double>::Det(J0[PointNumber]) * IntegrationWeight;
        }
        if ( mTotalDomainInitialSize < 0.0 )
        {
            std::stringstream ss;
            ss << "error on element -> " << this->Id() << std::endl;
            ss << ". Domain size can not be less than 0, mTotalDomainInitialSize = " << mTotalDomainInitialSize;
            KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" );
        }
        this->SetValue(GEOMETRICAL_DOMAIN_SIZE, mTotalDomainInitialSize);

        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize( GetGeometry().size(), dim, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            for ( unsigned int i = 0; i < dim; ++i )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
//        std::cout << "Element " << Id() << " mInitialDisp is initialized to " << mInitialDisp << std::endl;
        #endif

        mCurrentStrainVector.resize(integration_points.size());
        mOldStrainVector.resize(integration_points.size());
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            mCurrentStrainVector[PointNumber].resize(strain_size);
            noalias(mCurrentStrainVector[PointNumber]) = ZeroVector(strain_size);
            mOldStrainVector[PointNumber].resize(strain_size);
            noalias(mOldStrainVector[PointNumber]) = ZeroVector(strain_size);
        }

        InitializeMaterial();

        mIsInitialized = true;

        KRATOS_CATCH( "" )
    }

    /**
     * Initialization of the Material law at each integration point
     */
    void ExtrapolatedKinematicLinear::InitializeMaterial()
    {
    }

    void ExtrapolatedKinematicLinear::ResetConstitutiveLaw()
    {
    }

    /**
     * THIS is the main method here the integration in space (loop over the integration points) is done,
     * the algorithmic tangent and the (inner and outer) load vector is computed
     * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     * @param CalculateStiffnessMatrixFlag true: algorithmic tangent has to be computed
     * @param CalculateResidualVectorFlag true: load vector has to be computed
     */
    void ExtrapolatedKinematicLinear::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       bool CalculateStiffnessMatrixFlag,
                                       bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

//        std::cout << "start computing extrapolated element " << Id() << std::endl;

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();;
        unsigned int strain_size = dim * ( dim + 1 ) / 2;

        //this is the size of the elements stiffness matrix/force vector
        unsigned int mat_size = GetGeometry().size() * dim;

        //Initialize local variables
        Matrix B( strain_size, mat_size );
        Matrix TanC( strain_size, strain_size );
        Matrix TanCe( strain_size, strain_size );
        Vector StrainVector( strain_size );
        Vector StressVector( strain_size );
        Vector PhysicalStrainVector( strain_size );
        Vector PhysicalStressVector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );
        Matrix InvJ0(dim, dim);
        double DetJ0;

        //resize the LHS=StiffnessMatrix if its size is not correct
        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        ////// make a safety check
        if(integration_points.size() == 0)
        {
            std::cout << "!!!!WARNING!!!!The integration_points at extrapolated element " << Id() << " is null. Something must be wrong." << std::endl;
        }

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        // TODO check if the body forces, i.e. gravity is computed correctly in this extrapolation scheme

        //auxiliary terms
        const Vector& BodyForce = GetProperties()[BODY_FORCE];

        // the consitutive_law that this element relies on
        ConstitutiveLaw::Pointer pConstitutiveLaw = GetValue(CONSTITUTIVE_LAW);

        // extract the required value at the master consitutive_law
        TanCe = pConstitutiveLaw->GetValue(ELASTIC_TANGENT, TanCe);
        TanC = pConstitutiveLaw->GetValue(ALGORITHMIC_TANGENT, TanC);
        array_1d<double, 3> PhysicalIntPoint = GetValue(PHYSICAL_INTEGRATION_POINT_LOCAL);


        /****
        DEBUGGING
        ****/
        int parent_element_id, parent_integration_point_index;
        parent_element_id = pConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, parent_element_id);
        parent_integration_point_index = pConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, parent_integration_point_index);
//        KRATOS_WATCH(parent_element_id)
//        KRATOS_WATCH(parent_integration_point_index)

//        std::vector<std::size_t> equation_ids;
//        this->EquationIdVector(equation_ids, rCurrentProcessInfo);

//        if(parent_element_id == 37)
//        {
//            std::cout << "integration_points of extrapolated kinematic_element " << Id()
//                      << ", parent_element " << parent_element_id << ":" << std::endl;
//            for(std::size_t i = 0; i < integration_points.size(); ++i)
//            {
//                std::cout << "(" << integration_points[i].X() << ", "
//                                 << integration_points[i].Y() << ", "
//                                 << integration_points[i].Z() << "), "
//                                 << integration_points[i].Weight() << std::endl;
//            }

//            std::cout << "nodes of extrapolated kinematic_element " << Id() << ":";
//            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
//                std::cout << " " << GetGeometry()[i].Id();
//            std::cout << std::endl;

//            std::cout << "equation ids of extrapolated kinematic_element " << Id() << ":" << std::endl;
//            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
//            {
//                std::cout << "\t" << GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId()
//                          << " " << GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId()
//                          << " " << GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId()
//                          << std::endl;
//            }

//        }
        /*********************/

        // compute B operator at physical point
        Matrix DN_Dep;
        DN_Dep = GetGeometry().ShapeFunctionsLocalGradients(DN_Dep, PhysicalIntPoint);

        Matrix Jp(dim, dim);
        Matrix InvJp(dim, dim);
        double DetJp;
        Jp = GetGeometry().Jacobian(Jp, PhysicalIntPoint, DeltaPosition);

        Matrix Bp( strain_size, mat_size );
        Matrix DN_DXp( number_of_nodes, dim );
        MathUtils<double>::InvertMatrix( Jp, InvJp, DetJp );
        noalias( DN_DXp ) = prod( DN_Dep, InvJp );
        CalculateBoperator( Bp, DN_DXp );

        double beta = 1.0;
        if(GetProperties().Has(PHYSICAL_STRESS_OFFSET_PARAMETER))
        {
            beta = GetProperties()[PHYSICAL_STRESS_OFFSET_PARAMETER];
        }
//        KRATOS_WATCH(beta)

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            // update stress
            PhysicalStressVector = pConstitutiveLaw->GetValue(STRESSES, PhysicalStressVector);
            PhysicalStrainVector = pConstitutiveLaw->GetValue(STRAIN, PhysicalStrainVector);
//            noalias(StressVector) = PhysicalStressVector + prod(TanCe, StrainVector - PhysicalStrainVector);
            noalias(StressVector) = prod(TanCe, StrainVector) + beta * (PhysicalStressVector - prod(TanCe, PhysicalStrainVector));
//            KRATOS_WATCH(PhysicalStrainVector)
//            KRATOS_WATCH(PhysicalStressVector)

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight();

            //modify integration weight in case of 2D
            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculate stiffness matrix
                noalias( rLeftHandSideMatrix ) += ( IntToReferenceWeight * DetJ0 ) * (
                               prod( trans( B ), Matrix( prod( TanCe, B ) ) )
                      + beta * prod( trans( B ), Matrix( prod( TanC, Bp ) ) )
                      - beta * prod( trans( B ), Matrix( prod( TanCe, Bp ) ) )
                    );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //contribution of external forces
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight, DetJ0);

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, row( Ncontainer, PointNumber ), IntToReferenceWeight, DetJ0 );

                //contribution of internal forces
                AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, IntToReferenceWeight, DetJ0 );
            }
        }//loop over integration points
//        if(parent_element_id == 37)
//        {
//            std::cout << "At extrapolated element " << Id() << " with parent_element_id " << parent_element_id << ":" << std::endl;
//            KRATOS_WATCH(rLeftHandSideMatrix)
//            KRATOS_WATCH(rRightHandSideVector)
//            std::cout << "----------------------------------------" << std::endl;
//        }

//        for(std::size_t i = 0; i < equation_ids.size(); ++i)
//        {
//            if(equation_ids[i] == 21807)
//            {
//                std::cout << "extrapolated element " << Id() << " with parent_element_id " << parent_element_id
//                          << " contributes lhs = " << rLeftHandSideMatrix(i, i)
//                          << ", rhs = " << rRightHandSideVector(i)
//                          << " to row " << equation_ids[i]
//                          << std::endl;
//            }
//        }

        // modify the right hand side to account for prescribed displacement
        // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
        // // However, I have to temporarily disable it to keep the consistency.
        if(CalculateStiffnessMatrixFlag && CalculateResidualVectorFlag)
        {
            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) && dim == 3)
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

//        for(std::size_t i = 0; i < GetGeometry().size(); ++i)
//        {
//            if(GetGeometry()[i].Id() == 690)
//            {
//                std::cout << "extrapolated_kinematic_linear " << Id() << " contributes "
//                          << subslice(rLeftHandSideMatrix, 3*i, 3*i+1, 3*i+2, 3*i, 3*i+1, 3*i+2)
//                          << " to the stiffness matrix"
//                          << std::endl;
//            }
//        }


        KRATOS_CATCH( "" )
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = true
     * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    /**
     * THIS method is called from the scheme at the start of each solution step
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
    }

    void ExtrapolatedKinematicLinear::InitializeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
    }

    void ExtrapolatedKinematicLinear::FinalizeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
    }

    /**
     * THIS method is called from the scheme after each solution step, here the time step
     * start and end point variables can be transferred n --> n+1
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        for(std::size_t i = 0; i < mOldStrainVector.size(); ++i)
            noalias(mOldStrainVector[i]) = mCurrentStrainVector[i];
    }

    void ExtrapolatedKinematicLinear::MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", __FUNCTION__)
    }

    void ExtrapolatedKinematicLinear::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int mat_size = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != mat_size )
            rMassMatrix.resize( mat_size, mat_size, false );

        rMassMatrix = ZeroMatrix( mat_size, mat_size );

        double TotalMass = 0.0;
        if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
        {
            TotalMass = mTotalDomainInitialSize * GetValue(DENSITY);
        }
        else
            TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

        if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

        for ( unsigned int i = 0; i < NumberOfNodes; ++i )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; ++j )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;
            }
        }

        KRATOS_CATCH( "" )
    }

    void ExtrapolatedKinematicLinear::DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", __FUNCTION__)
    }

    void ExtrapolatedKinematicLinear::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dim;

        if ( rDampMatrix.size1() != mat_size )
            rDampMatrix.resize( mat_size, mat_size, false );

        noalias( rDampMatrix ) = ZeroMatrix( mat_size, mat_size );

        Matrix StiffnessMatrix = ZeroMatrix( mat_size, mat_size );

        Vector RHS_Vector = ZeroVector( mat_size );

        //rayleigh damping
        CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

        double alpha = 0.001, beta = 0.001;

        if(GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
        {
            alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        }

        if(GetProperties().Has(RAYLEIGH_DAMPING_BETA))
        {
            beta = GetProperties()[RAYLEIGH_DAMPING_BETA];
        }

        CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

        noalias( rDampMatrix ) = alpha * rDampMatrix + beta * StiffnessMatrix;

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void ExtrapolatedKinematicLinear::GetValuesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )    values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }
    }


    //************************************************************************************
    //************************************************************************************
    void ExtrapolatedKinematicLinear::GetFirstDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )   values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
        }
    }

    //************************************************************************************
    //************************************************************************************
    void ExtrapolatedKinematicLinear::GetSecondDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size ) values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }

    /**
     * returns the used integration method
     */
    ExtrapolatedKinematicLinear::IntegrationMethod ExtrapolatedKinematicLinear::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a Vector containing the EquationIds of the differnt Dofs is created
     * @param rResult Vector of the EquationIds
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::EquationIdVector( EquationIdVectorType& rResult,
            const ProcessInfo& CurrentProcessInfo ) const
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

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a Container containing the pointers rto the DOFs of this element is created
     * @param ElementalDofList Container with of the DOFs associated with the nodes
     *                           of this element
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::GetDofList( DofsVectorType& ElementalDofList,
                const ProcessInfo& CurrentProcessInfo ) const
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            if(dim == 3)
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    /**
     * Adds the Body Forces to the load vector
     * @param R RHS Vector
     * @param N_DISP shape function values at the current integration points
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    inline void ExtrapolatedKinematicLinear::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector gravity( dim );

        double density = 0.0;
        if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        {
            noalias( gravity ) = GetValue(GRAVITY);
            density = GetValue(DENSITY);
        }
        else
        {
            noalias( gravity ) = GetProperties()[GRAVITY];
            density = GetProperties()[DENSITY];
        }

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    inline void ExtrapolatedKinematicLinear::CalculateAndAdd_ExtForceContribution(
            const Vector& N,
            const ProcessInfo& CurrentProcessInfo,
            const Vector& BodyForce,
            VectorType& rRightHandSideVector,
            double weight,
            double detJ
            )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; ++j )
                rRightHandSideVector[index + j] += weight * detJ * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Internal Forces to the load vector
     * @param R RHS Vector
     * @param B_Operator B-Operator at the current integration point
     * @param StressVector current stress vector
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
//    void ExtrapolatedKinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
//    {
//        KRATOS_TRY

//        unsigned int dim = GetGeometry().WorkingSpaceDimension();
//        unsigned int strain_size = dim * (dim + 1) / 2;
//
//        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
//        {
//            for ( unsigned int i = 0; i < dim; ++i )
//            {
//                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
//                {
//                    R( prim * dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight );
//                }
//            }
//        }
//        //         noalias(R) -= detJ * Weight * prod(trans(B_Operator), StressVector);
//
//        KRATOS_CATCH( "" )
//    }

    void ExtrapolatedKinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        Vector InternalForces(3);

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                InternalForces(i) = 0.0;
                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
                {
                    InternalForces(i) += B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight;
                }

                R( prim * dim + i ) -= InternalForces(i);
            }
//            GetGeometry()[prim].GetSolutionStepValue( REACTION ) += InternalForces;
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Contribution of the current quadrature point to the load vector
     * @param K LHS Matrix
     * @param tan_C 6*6 algorithmic tangent of the materia law (derivation of stresses
     *               regarding strains
     * @param B_Operator B-Operator at the current integration point
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    void ExtrapolatedKinematicLinear::CalculateStiffnesMatrix( Matrix& K, const Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2 ;

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                for ( unsigned int sec = 0; sec < GetGeometry().size(); ++sec )
                {
                    for ( unsigned int j = 0; j < dim; ++j )
                    {
                        for ( unsigned int alpha = 0; alpha < strain_size; ++alpha )
                            for ( unsigned int beta = 0; beta < strain_size; ++beta )
                                K( prim*dim + i, sec*dim + j ) += B_Operator( alpha, dim * prim + i )
                                    * tan_C( alpha, beta ) * B_Operator( beta, dim * sec + j ) * detJ * Weight;
                    }
                }
            }
        }

//        noalias(K) -= prod( trans(B_Operator), (Weight * detJ) * Matrix(prod(C, B_Operator)) );

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the stress vector and the algorithmic tangent at the current quadrature points
     * regarding the current strain vector
     * @param StressVector output of the method, Cauchy stress vector
     * @param tanC_U algorithmic tangen (derivative Cauchy stress vector/StrainVector) [strain_size*strain_size]
     * @param StrainVector output: current strain vector
     * @param B_Operator current B-operator
     * @param PointNumber number of the current integration point
     * @param CurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateStressAndTangentialStiffness( Vector& StressVector, Matrix& tanC_U,
            Vector& StrainVector, const Matrix& B_Operator, int PointNumber,
            const ProcessInfo& CurrentProcessInfo )
    {
    }

    /**
     * Computes the strain vector
     */
    void ExtrapolatedKinematicLinear::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
    {
        KRATOS_TRY
        unsigned int Dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = Dim * (Dim + 1) / 2;
        noalias( StrainVector ) = ZeroVector( strain_size );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            for ( unsigned int item = 0; item < strain_size; ++item )
                for ( unsigned int dim = 0; dim < Dim; ++dim )
                    StrainVector[item] += B( item, Dim * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the B-Operator at the current quadrature point
     * @param B_Operator current B-operator
     * @param DN_DX shape function values at the current integration point
     */
    void ExtrapolatedKinematicLinear::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

        if(dim == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*2 ) = DN_DX( i, 0 );
                B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
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

    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == PK2_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        return;
    }

    /**
     * Calculate Matrix Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     *
     * Calculate Vector Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if( rVariable == PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
        {
            ConstitutiveLaw::Pointer pConstitutiveLaw = GetValue(CONSTITUTIVE_LAW);

             if ( rValues.size() != 1 )
                rValues.resize( 1 );

            rValues[0] = pConstitutiveLaw->GetValue(THREED_STRESSES, rValues[0]);

            return;
        }

        if( rVariable == THREED_STRESSES )
        {
            ConstitutiveLaw::Pointer pConstitutiveLaw = GetValue(CONSTITUTIVE_LAW);

            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

             if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            for (std::size_t i = 0; i < integration_points.size(); ++i)
                rValues[i] = pConstitutiveLaw->GetValue(THREED_STRESSES, rValues[i]);

            return;
        }
    }

    /**
     * Calculate double Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if( rVariable == DISPLACEMENT )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = ZeroVector(3);
                for(std::size_t i = 0; i < GetGeometry().size(); ++i)
                {
                    const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                    noalias(rValues[point]) += Ncontainer(point, i) * displacement;
                }
            }

            return;
        }

        if( rVariable == INTEGRATION_POINT_GLOBAL )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
            }

            return;
        }

        if( rVariable == INTEGRATION_POINT_LOCAL )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = integration_points[point];
            }

            return;
        }

        if( rVariable == PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )
        {
            if ( rValues.size() != 1 )
                rValues.resize( 1 );

            Vector Ncontainer;
            noalias(Ncontainer) = GetGeometry().ShapeFunctionsValues( Ncontainer, GetValue(PHYSICAL_INTEGRATION_POINT_LOCAL) );

            noalias(rValues[0]) = ZeroVector(3);
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
            {
                const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                noalias(rValues[0]) += Ncontainer(i) * displacement;
            }

            return;
        }

        if( rVariable == PHYSICAL_INTEGRATION_POINT_LOCAL )
        {
            if ( rValues.size() != 1 )
                rValues.resize( 1 );

            rValues[0] = GetValue(PHYSICAL_INTEGRATION_POINT_LOCAL);

            return;
        }

        if( rVariable == PHYSICAL_INTEGRATION_POINT_GLOBAL )
        {
            if ( rValues.size() != 1 )
                rValues.resize( 1 );

            rValues[0] = GetGeometry().GlobalCoordinates(rValues[0], GetValue(PHYSICAL_INTEGRATION_POINT_LOCAL));

            return;
        }
    }

    /**
     * Calculate double Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if( rVariable == JACOBIAN_0 )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points
                    = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            // initializing the Jacobian in the reference configuration
            GeometryType::JacobiansType J0;
            Matrix DeltaPosition(GetGeometry().size(), 3);

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
            }

            J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

            // compute the Jacobian
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = MathUtils<double>::Det(J0[i]);
            }
        }
        else if( rVariable == INTEGRATION_WEIGHT )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points
                    = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = integration_points[i].Weight();
            }
        }
        else if( rVariable == SUBCELL_DOMAIN_SIZE )
        {
            if ( rValues.size() != 1 )
                rValues.resize( 1 );

            rValues[0] = GetValue(SUBCELL_DOMAIN_SIZE);
        }
        else if( rVariable == VON_MISES_STRESS )
        {
            ConstitutiveLaw::Pointer pConstitutiveLaw = GetValue(CONSTITUTIVE_LAW);

            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

             if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            Vector stresses(6);
            for (std::size_t i = 0; i < integration_points.size(); ++i)
            {
                stresses = pConstitutiveLaw->GetValue(THREED_STRESSES, stresses);
                double p = (stresses[0] + stresses[1] + stresses[2]) / 3;
                double sxx = stresses[0] - p;
                double syy = stresses[1] - p;
                double szz = stresses[2] - p;
                double sxy = stresses[3];
                double syz = stresses[4];
                double sxz = stresses[5];
                rValues[i] = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
            }
        }
    }

    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void ExtrapolatedKinematicLinear::CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    /**
     * Set a Matrix Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector of the values on the quadrature points
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable,
            const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    /**
     * Set a Vector Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector of the values on the quadrature points
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable,
            const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    /**
     * Set a Double Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValue value on the quadrature points
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
            const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }
    /**
     * Set an Int Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValue value on the quadrature points
     * @param rCurrentProcessInfo
     */
    void ExtrapolatedKinematicLinear::SetValuesOnIntegrationPoints( const Variable<int>& rVariable,
            const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void ExtrapolatedKinematicLinear::SetValuesOnIntegrationPoints( const Kratos::Variable< ConstitutiveLaw::Pointer >& rVariable,
            const std::vector< ConstitutiveLaw::Pointer >& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
    }

    int ExtrapolatedKinematicLinear::Check( const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

//	    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        if ( this->Id() < 1 )
        {
            KRATOS_THROW_ERROR( std::logic_error, "Element found with Id 0 or negative, Id() =", Id() );
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property", this->GetProperties().Id() );
        }

        //Verify that the body force is defined
        if ( this->GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property", this->GetProperties().Id() )
        }

        return 0;

        KRATOS_CATCH( "" );
    }

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW

