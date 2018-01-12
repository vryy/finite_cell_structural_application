// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/condition.h"
#include "custom_python/add_custom_conditions_to_python.h"
#include "custom_conditions/ghost_penalty_displacement_gradient_condition.h"
#include "custom_conditions/ghost_penalty_stress_condition.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellStructuralApplication_AddCustomConditionsToPython()
{

    class_< GhostPenaltyDisplacementGradientCondition, GhostPenaltyDisplacementGradientCondition::Pointer, bases< GhostPenaltyCondition > >
    ("GhostPenaltyDisplacementGradientCondition", init<>() )
    .def(self_ns::str(self))
    ;

    class_< GhostPenaltyStressCondition, GhostPenaltyStressCondition::Pointer, bases< GhostPenaltyCondition > >
    ("GhostPenaltyStressCondition", init<>() )
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

}  // namespace Kratos.

