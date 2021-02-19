//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_GHOST_PENALTY_STRESS_CONDITION_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_STRESS_CONDITION_H_INCLUDED


// External includes

// Project includes
#include "finite_cell_application/custom_conditions/ghost_penalty_condition.h"

namespace Kratos
{

/**
 * Ghost penalty implementation for linear elasticity
 * This class add a weak variational term to enforce the continuity of the stress along the boundary
 * The term is \int_edge gamma*[|d\epsilon|]*[|\sigma|]
 */
class GhostPenaltyStressCondition : public GhostPenaltyCondition
{
    public:
        // Counted pointer of GhostPenaltyStressCondition
        KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyStressCondition);

        typedef GeometryType::IntegrationPointType IntegrationPointType;

        /**
         * Default constructor.
         */
        GhostPenaltyStressCondition();

        GhostPenaltyStressCondition( IndexType NewId, GeometryType::Pointer pGeometry);

        GhostPenaltyStressCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        GhostPenaltyStressCondition( IndexType NewId, GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement, Element::Pointer pMasterElement, PropertiesType::Pointer pProperties );

        /**
         * Destructor.
         */
        virtual ~GhostPenaltyStressCondition();

        Condition::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties) const final;

        /**
         * Operations.
         */

        void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                                const ProcessInfo& rCurrentProcessInfo) const final;

        void GetDofList( DofsVectorType& ConditionalDofList,
                                const ProcessInfo& CurrentProcessInfo) const final;

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

    private:

        bool mIsInitialized;
        std::vector<ConstitutiveLaw::Pointer> mSlaveConstitutiveLaws;
        std::vector<ConstitutiveLaw::Pointer> mMasterConstitutiveLaws;

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

        void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX ) const;

        void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector ) const;

}; // Class GhostPenaltyStressCondition

}  // namespace Kratos.


#endif // KRATOS_GHOST_PENALTY_STRESS_CONDITION_H_INCLUDED defined

