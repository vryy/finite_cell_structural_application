// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "finite_cell_application/custom_algebra/function.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/immersed_boundary_load_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellStructuralApplication_AddCustomUtilitiesToPython()
{

    void(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForce)(
        ModelPart&, const FunctionR1R3&, const double&, const double&, const int&,
        const FunctionR1R3&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce;

    class_<ImmersedBoundaryLoadUtility, ImmersedBoundaryLoadUtility::Pointer, boost::noncopyable>
    ("ImmersedBoundaryLoadUtility", init<>())
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForce)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

