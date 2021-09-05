// see finite_cell_structural_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Feb 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_IMMERSED_POINT_FORCE_H_INCLUDED )
#define  KRATOS_IMMERSED_POINT_FORCE_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"

namespace Kratos
{

/**
 */
class ImmersedPointForce : public Condition
{
    public:
        // Counted pointer of ImmersedPointForce
        KRATOS_CLASS_POINTER_DEFINITION(ImmersedPointForce);

        /**
         * Default constructor.
         */
        ImmersedPointForce();
        ImmersedPointForce( IndexType NewId, GeometryType::Pointer pGeometry);
        ImmersedPointForce( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
        ImmersedPointForce( IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        const double& rWeight,
                        const array_1d<double, 3>& rForce,
                        Element::Pointer pMasterElement,
                        Point<3>& rMasterLocalPoint );
        ImmersedPointForce( IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
                        const double& rWeight,
                        const array_1d<double, 3>& rForce,
                        Element::Pointer pMasterElement,
                        Point<3>& rMasterLocalPoint );

        /**
         * Destructor.
         */
        virtual ~ImmersedPointForce();

        /**
         * Operations.
         */

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const final;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo);

        void CalculateDampingMatrix( MatrixType& rDampMatrix,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo) const final;

        void GetDofList( DofsVectorType& ConditionalDofList,
                         const ProcessInfo& CurrentProcessInfo) const final;

        void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

        void SetMagnitude(const double& P);

        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();

        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        //virtual void PrintData(std::ostream& rOStream) const;

    protected:


    private:

        friend class Serializer;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

        array_1d<double, 3> mPointForceDirection;
        double mWeight; // Weight * DetJ
        Point<3> mMasterLocalPoint;
        Element::Pointer mpMasterElement;

}; // Class ImmersedPointForce

}  // namespace Kratos.


#endif // KRATOS_IMMERSED_POINT_FORCE_H_INCLUDED defined

