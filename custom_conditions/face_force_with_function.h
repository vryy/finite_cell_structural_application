// see finite_cell_application/LICENSE.txt
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: Aug 26, 2019 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_FACE_FORCE_WITH_FUNCTION_H_INCLUDED )
#define  KRATOS_FACE_FORCE_WITH_FUNCTION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{

class FaceForceWithFunction : public Condition
{
public:

    // Counted pointer of FaceForceWithFunction
    KRATOS_CLASS_POINTER_DEFINITION( FaceForceWithFunction );

    typedef Condition BaseType;

    typedef BaseType::GeometryType::PointType NodeType;

    typedef NodeType::PointType PointType;

    // Constructor void
    FaceForceWithFunction();

    // Constructor using an array of nodes
    FaceForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    FaceForceWithFunction( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~FaceForceWithFunction();


    // Name Operations

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties ) const;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo );

    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo );

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo );

//    void GetValuesVector(
//        Vector& values,
//        int Step = 0 );

//    void GetFirstDerivativesVector(
//        Vector& values,
//        int Step = 0 );

//    void GetSecondDerivativesVector(
//        Vector& values,
//        int Step = 0 );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


protected:


private:
    ///@name Static Member Variables

    /// privat variables


    // privat name Operations

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class FaceForceWithFunction.

} // namespace Kratos.

#endif // KRATOS_FACE_FORCE_WITH_FUNCTION_H_INCLUDED defined
