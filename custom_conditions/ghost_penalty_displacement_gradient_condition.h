//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_GHOST_PENALTY_DISPLACEMENT_GRADIENT_CONDITION_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_DISPLACEMENT_GRADIENT_CONDITION_H_INCLUDED


// External includes

// Project includes
#include "finite_cell_application/custom_conditions/ghost_penalty_condition.h"

namespace Kratos
{

/**
 * Ghost penalty implementation for linear elasticity
 * This class add a weak variational term to enforce the continuity of the gradient of displacement along the boundary
 */
class GhostPenaltyDisplacementGradientCondition : public GhostPenaltyCondition
{
    public:
        // Counted pointer of GhostPenaltyDisplacementGradientCondition
        KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyDisplacementGradientCondition);

        /**
         * Default constructor.
         */
        GhostPenaltyDisplacementGradientCondition();

        GhostPenaltyDisplacementGradientCondition( IndexType NewId, GeometryType::Pointer pGeometry);

        GhostPenaltyDisplacementGradientCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        GhostPenaltyDisplacementGradientCondition( IndexType NewId, GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement, Element::Pointer pMasterElement, PropertiesType::Pointer pProperties );

        /**
         * Destructor.
         */
        virtual ~GhostPenaltyDisplacementGradientCondition();

        Condition::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties) const final;

        /**
         * Operations.
         */

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
         */
        std::string Info() const final
        {
            std::stringstream buffer;
            buffer << "GhostPenaltyDisplacementGradientCondition #" << Id();
            return buffer.str();
        }

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

}; // Class GhostPenaltyDisplacementGradientCondition

}  // namespace Kratos.


#endif // KRATOS_GHOST_PENALTY_DISPLACEMENT_GRADIENT_CONDITION_H_INCLUDED defined

