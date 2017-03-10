// see finite_cell_structural_application/LICENSE.txt
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Feb 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes 

// External includes 

// Project includes 
#include "custom_conditions/immersed_point_force.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
ImmersedPointForce::ImmersedPointForce()
{
}

ImmersedPointForce::ImmersedPointForce( IndexType NewId, 
                              GeometryType::Pointer pGeometry)
: Condition( NewId, pGeometry )
{
}

ImmersedPointForce::ImmersedPointForce( IndexType NewId, 
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Condition( NewId, pGeometry, pProperties )
{
}

ImmersedPointForce::ImmersedPointForce( IndexType NewId,
                                    GeometryType::Pointer pGeometry,
                                    const array_1d<double, 3>& rForce,  
                                    Element::Pointer pMasterElement,
                                    Point<3>& rMasterLocalPoint )
: Condition( NewId, pGeometry )
, mPointForce(rForce), mMasterLocalPoint(rMasterLocalPoint), mpMasterElement(pMasterElement)
{
}

ImmersedPointForce::ImmersedPointForce( IndexType NewId,
                                    GeometryType::Pointer pGeometry,
                                    PropertiesType::Pointer pProperties,
                                    const array_1d<double, 3>& rForce,  
                                    Element::Pointer pMasterElement,
                                    Point<3>& rMasterLocalPoint )
: Condition( NewId, pGeometry, pProperties )
, mPointForce(rForce), mMasterLocalPoint(rMasterLocalPoint), mpMasterElement(pMasterElement)
{
}

/**
 * Destructor. Never to be called manually
 */
ImmersedPointForce::~ImmersedPointForce()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer ImmersedPointForce::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ImmersedPointForce(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Condition::Pointer ImmersedPointForce::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ImmersedPointForce(NewId, pGeom, pProperties));
}

void ImmersedPointForce::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

//************************************************************************************ 
//************************************************************************************
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void ImmersedPointForce::CalculateRightHandSide( VectorType& rRightHandSideVector, 
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
void ImmersedPointForce::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
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
 * This function calculates all system contributions due to the contact problem
 * with regard to the current master and slave partners.
 * All Conditions are assumed to be defined in 2D/3D space and having 2/3 DOFs per node 
 */
void ImmersedPointForce::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                  VectorType& rRightHandSideVector, 
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    unsigned int Dim = mpMasterElement->GetGeometry().WorkingSpaceDimension();
    unsigned int MasterNN = mpMasterElement->GetGeometry().size();
    unsigned int MatSize = MasterNN*Dim;

    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if(rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize, false);
        noalias(rRightHandSideVector) = ZeroVector(MatSize);
    }

    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
    {
        if(rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize)
            rLeftHandSideMatrix.resize(MatSize, MatSize, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); 
    }

    if( ( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) ) )
    {
        return;
    }

    if( !CalculateStiffnessMatrixFlag && !CalculateResidualVectorFlag )
        return;

    Vector ShapeFunctionValuesOnMaster;
    ShapeFunctionValuesOnMaster = mpMasterElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnMaster, mMasterLocalPoint);

    for( unsigned int node = 0; node < MasterNN; ++node )
    {
        for( unsigned int dim = 0; dim < Dim; ++dim )
            rRightHandSideVector[Dim*node + dim] += mPointForce[dim] * ShapeFunctionValuesOnMaster[node];
    }
//KRATOS_WATCH(mpMasterElement->Id())
//KRATOS_WATCH(rRightHandSideVector)
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.    
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void ImmersedPointForce::EquationIdVector( EquationIdVectorType& rResult, 
                                      ProcessInfo& CurrentProcessInfo)
{
    //determining size of DOF list
    unsigned int Dim = mpMasterElement->GetGeometry().WorkingSpaceDimension();
    unsigned int MasterNN = mpMasterElement->GetGeometry().size();
    unsigned int ndofs = Dim*(MasterNN);
    unsigned int index;

    if(rResult.size() != ndofs)
        rResult.resize(ndofs, false);

    for( unsigned int node = 0; node < MasterNN; ++node )
    {
        index = node*Dim;
        rResult[index  ] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index+1] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
        if(Dim == 3)
            rResult[index+2] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
    }
//KRATOS_WATCH(mpMasterElement->Id())
//std::cout << "EquationId:";
//for(int i = 0; i < ndofs; ++i)
//    std::cout << " " << rResult[i];
//std::cout << std::endl;
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
void ImmersedPointForce::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    ConditionalDofList.resize(0);

    //determining size of DOF list
    unsigned int Dim = mpMasterElement->GetGeometry().WorkingSpaceDimension();
    unsigned int MasterNN = mpMasterElement->GetGeometry().size();

    for( unsigned int node = 0; node < MasterNN; ++node )
    {
        ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X) );
        ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y) );
        if(Dim == 3)
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z) );
    }
}

} // Namespace Kratos

