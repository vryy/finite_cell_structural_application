/*
see finite_cell_structural_application/LICENSE.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Mar 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_EXTRAPOLATED_CONSTANT_STRESS_KINEMATIC_LINEAR_INCLUDED )
#define  KRATOS_EXTRAPOLATED_CONSTANT_STRESS_KINEMATIC_LINEAR_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
// #include "includes/constitutive_law_new.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

Using for moment-fitted subcell scheme, constant stress extrapolation from physical point
stress at Gauss point = stress at physical point + alpha*Ce(strain at Gauss point - strain at physical point)
alpha = 0: stress at Gauss point = stress at physical point (constant stress over sub-cell)
alpha = 1: for linear elastic, it is always: stress at Gauss point = Ce*strain at Gauss point
 */

class ExtrapolatedConstantStressKinematicLinear : public Element
{

public:
    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    KRATOS_CLASS_POINTER_DEFINITION( ExtrapolatedConstantStressKinematicLinear );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExtrapolatedConstantStressKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry );
    ExtrapolatedConstantStressKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~ExtrapolatedConstantStressKinematicLinear();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    IntegrationMethod GetIntegrationMethod() const final;

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const final;

    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const final;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

    void ResetConstitutiveLaw() final;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) final;

    void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const final;

    void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const final;

    void MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) final;

    void DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) final;

    void FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo ) final;

    void InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo ) final;

    void InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) final;

    void FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void SetValuesOnIntegrationPoints( const Variable<double>& rVariable, const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void SetValuesOnIntegrationPoints( const Variable<int>& rVariable, const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable, const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void SetValuesOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, const std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void GetValuesVector( Vector& values, int Step = 0 ) const final;

    void GetFirstDerivativesVector( Vector& values, int Step = 0 ) const final;

    void GetSecondDerivativesVector( Vector& values, int Step = 0 ) const final;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) const final;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    bool mIsInitialized;
    Matrix mInitialDisp;
    IntegrationMethod mThisIntegrationMethod;
    double mTotalDomainInitialSize;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeMaterial();

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    ExtrapolatedConstantStressKinematicLinear() {}

    void save( Serializer& rSerializer ) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  Element );
        rSerializer.save( "mIsInitialized", mIsInitialized );
        rSerializer.save( "mInitialDisp", mInitialDisp );
//        rSerializer.save( "mThisIntegrationMethod", mThisIntegrationMethod );
        rSerializer.save( "mTotalDomainInitialSize", mTotalDomainInitialSize );
    }

    void load( Serializer& rSerializer ) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  Element );
        rSerializer.load( "mIsInitialized", mIsInitialized );
        rSerializer.load( "mInitialDisp", mInitialDisp );
//        rSerializer.load( "mThisIntegrationMethod", mThisIntegrationMethod );
        rSerializer.load( "mTotalDomainInitialSize", mTotalDomainInitialSize );
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{
    /** K += weight*Btrans*D*B */
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    void CalculateBodyForces( Vector& BodyForce, const ProcessInfo& CurrentProcessInfo );

    void InitializeVariables();

    void InitializeMaterial( std::vector<std::vector<Matrix> >& C );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE FORCEVECTORS DISPLACEMENT

    void AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ );

    void CalculateAndAdd_ExtForceContribution( const Vector& N, const ProcessInfo& CurrentProcessInfo,
                                               const Vector& BodyForce, VectorType& rRightHandSideVector,
                                               double weight, double detJ);

    void AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ );

    void CalculateStiffnesMatrix( Matrix& K, const Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ );

    void CalculateStressAndTangentialStiffness( Vector& StressVector, Matrix& tanC_U,
                                                Vector& StrainVector, const Matrix& B_Operator,
                                                int PointNumber, const ProcessInfo& CurrentProcessInfo );

    void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector );

//     Matrix CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, int PointNumber );

    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX );


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //ExtrapolatedConstantStressKinematicLinear& operator=(const ExtrapolatedConstantStressKinematicLinear& rOther);

    /// Copy constructor.
    //ExtrapolatedConstantStressKinematicLinear(const ExtrapolatedConstantStressKinematicLinear& rOther);


    ///@}

}; // Class ExtrapolatedConstantStressKinematicLinear

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   ExtrapolatedConstantStressKinematicLinear& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                   const ExtrapolatedConstantStressKinematicLinear& rThis)
           {
                   rThis.PrintInfo(rOStream);
                   rOStream << std::endl;
                   rThis.PrintData(rOStream);

                   return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_EXTRAPOLATED_CONSTANT_STRESS_KINEMATIC_LINEAR_INCLUDED defined


