//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_GHOST_PENALTY_ELASTIC_CONDITION_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_ELASTIC_CONDITION_H_INCLUDED


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
 */
class GhostPenaltyElasticCondition : public GhostPenaltyCondition
{
    public:
        // Counted pointer of GhostPenaltyElasticCondition
        KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyElasticCondition);

        /**
         * Default constructor.
         */
        GhostPenaltyElasticCondition();

        GhostPenaltyElasticCondition( IndexType NewId, GeometryType::Pointer pGeometry);

        GhostPenaltyElasticCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        GhostPenaltyElasticCondition( IndexType NewId, GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement, Element::Pointer pMasterElement, PropertiesType::Pointer pProperties );

        /**
         * Destructor.
         */
        virtual ~GhostPenaltyElasticCondition();

        virtual Condition::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties);

        /**
         * Operations.
         */

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

        virtual void Initialize();

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

        IntegrationMethod mThisIntegrationMethod;

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

}; // Class GhostPenaltyElasticCondition

}  // namespace Kratos.


#endif // KRATOS_GHOST_PENALTY_CONDITION_H_INCLUDED defined

