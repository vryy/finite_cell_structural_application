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
#include "finite_cell_application/custom_algebra/function/function.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/immersed_boundary_load_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellStructuralApplication_AddCustomUtilitiesToPython()
{

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForce)(
        ModelPart&, const FunctionR1R3::Pointer&, const double&, const double&, const int&,
        const FunctionR1R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<0>;

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceWithBin)(
        ModelPart&, const FunctionR1R3::Pointer&, const double&, const double&, const int&,
        const FunctionR1R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<1>;

    void(ImmersedBoundaryLoadUtility::*pointer_to_Clean)(
        ModelPart&, ModelPart::ConditionsContainerType&, const int&) const = &ImmersedBoundaryLoadUtility::Clean;

    class_<ImmersedBoundaryLoadUtility, ImmersedBoundaryLoadUtility::Pointer, boost::noncopyable>
    ("ImmersedBoundaryLoadUtility", init<>())
    .def("InitializeBinning", &ImmersedBoundaryLoadUtility::InitializeBinning)
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForce)
    .def("SetupImmersedPointForceWithBin", pointer_to_SetupImmersedPointForceWithBin)
    .def("Clean", pointer_to_Clean)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

