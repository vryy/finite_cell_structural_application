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
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
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

        virtual Condition::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties);

        /**
         * Operations.
         */

        virtual void Initialize();

        /**
         * Calculates the local system contributions for this contact element
         */
        virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     ProcessInfo& rCurrentProcessInfo);

        virtual void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo);

        virtual void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo);

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

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

        void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX ) const;

        void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector ) const;

}; // Class GhostPenaltyStressCondition

}  // namespace Kratos.


#endif // KRATOS_GHOST_PENALTY_STRESS_CONDITION_H_INCLUDED defined

