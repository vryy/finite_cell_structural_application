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
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
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
                        const array_1d<double, 3>& rForce,
                        Element::Pointer pMasterElement,
                        Point<3>& rMasterLocalPoint );
        ImmersedPointForce( IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
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

        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

        virtual Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                   VectorType& rRightHandSideVector, 
                                   ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                                     ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector( EquationIdVectorType& rResult, 
                               ProcessInfo& rCurrentProcessInfo);

        void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo);

        void Initialize();
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

        array_1d<double, 3> mPointForce;
        Point<3> mMasterLocalPoint;
        Element::Pointer mpMasterElement;

}; // Class ImmersedPointForce 

}  // namespace Kratos.
  

#endif // KRATOS_IMMERSED_POINT_FORCE_H_INCLUDED defined 

