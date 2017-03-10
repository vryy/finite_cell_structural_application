//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Feb 10, 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FINITE_CELL_STRUCTURAL_APPLICATION_H_INCLUDED)
#define KRATOS_FINITE_CELL_STRUCTURAL_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "structural_application/custom_elements/kinematic_linear.h"
#include "custom_conditions/line_force_with_function.h"
#include "custom_conditions/immersed_point_force.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KratosFiniteCellStructuralApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{
        
        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosFiniteCellStructuralApplication);

        ///@}
        ///@name Life Cycle
        ///@{ 

        /// Default constructor.
        KratosFiniteCellStructuralApplication();

        /// Destructor.
        virtual ~KratosFiniteCellStructuralApplication(){}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register();

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
        virtual std::string Info() const
        {
            return "Finite Cell method for structural simulation";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << "in KratosFiniteCellStructuralApplication:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
        }


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


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


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

        const KinematicLinear mKinematicLinearFiniteCell2D3N;
        const KinematicLinear mKinematicLinearFiniteCell2D6N;
        const KinematicLinear mKinematicLinearFiniteCell2D4N;
        const KinematicLinear mKinematicLinearFiniteCell2D8N;
        const KinematicLinear mKinematicLinearFiniteCell2D9N;
        const KinematicLinear mKinematicLinearFiniteCell3D4N;
        const KinematicLinear mKinematicLinearFiniteCell3D10N;
        const KinematicLinear mKinematicLinearFiniteCell3D8N;
        const KinematicLinear mKinematicLinearFiniteCell3D20N;
        const KinematicLinear mKinematicLinearFiniteCell3D27N;

        const LineForceWithFunction mLineForceWithFunction2D2N;
        const LineForceWithFunction mLineForceWithFunction2D3N;
        const LineForceWithFunction mLineForceWithFunction3D2N;
        const LineForceWithFunction mLineForceWithFunction3D3N;

        const ImmersedPointForce mImmersedPointForce3D;

        ///@}
        ///@name Private Operators
        ///@{


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
        KratosFiniteCellStructuralApplication& operator=(KratosFiniteCellStructuralApplication const& rOther);

        /// Copy constructor.
        KratosFiniteCellStructuralApplication(KratosFiniteCellStructuralApplication const& rOther);


        ///@}

    }; // Class KratosFiniteCellStructuralApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos

#endif // KRATOS_FINITE_CELL_STRUCTURAL_APPLICATION_H_INCLUDED defined

