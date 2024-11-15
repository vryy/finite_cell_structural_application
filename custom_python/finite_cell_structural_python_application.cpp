//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Feb 10, 2017 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)

// Project includes
#include "includes/define_python.h"
#include "finite_cell_structural_application.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_conditions_to_python.h"

namespace Kratos
{

namespace Python
{

    BOOST_PYTHON_MODULE(KratosFiniteCellStructuralApplication)
    {

        using namespace boost::python;

        class_<KratosFiniteCellStructuralApplication, KratosFiniteCellStructuralApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosFiniteCellStructuralApplication");

        FiniteCellStructuralApplication_AddCustomUtilitiesToPython();
        FiniteCellStructuralApplication_AddCustomConditionsToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_MAGNITUDE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PHYSICAL_STRESS_OFFSET_PARAMETER )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRESS_STABILIZATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( GHOST_PENALTY_STABILIZATION_ORDER )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( GHOST_PENALTY_STABILIZATION_FACTOR )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
