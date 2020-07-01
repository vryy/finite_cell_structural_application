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
#include "brep_application/custom_algebra/function/function.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "finite_cell_application/custom_utilities/immersed_boundary_utility.h"
#include "custom_utilities/immersed_boundary_load_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<class TEntityType, int TSearchType>
boost::python::list SetupImmersedPointForceElement2(
        ImmersedBoundaryLoadUtility& rDummy,
        ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pMasterElements,
        typename TEntityType::Pointer& pEntity,
        std::size_t lastNodeId,
        std::size_t lastCondId,
        Properties::Pointer pProperties,
        const int& integration_order,               // integration order
        const FunctionR3R3::Pointer& pLoadFunction, // the distributed load function
        const bool& export_point,
        const int& echo_level)
{
    ModelPart::ConditionsContainerType NewConds = rDummy.SetupImmersedPointForce<TEntityType, TSearchType>(
            r_model_part,
            pMasterElements,
            pEntity,
            lastNodeId,
            lastCondId,
            pProperties,
            integration_order,
            pLoadFunction,
            export_point,
            echo_level);

    boost::python::list Output;
    Output.append(lastNodeId);
    Output.append(lastCondId);
    Output.append(NewConds);

    return Output;
}

void FiniteCellStructuralApplication_AddCustomUtilitiesToPython()
{

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceCurve)(
        ModelPart&, const FunctionR1R3::Pointer&, const double&, const double&, Properties::Pointer, const int&,
        const FunctionR1R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<0>;

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceCurveWithBin)(
        ModelPart&, const FunctionR1R3::Pointer&, const double&, const double&, Properties::Pointer, const int&,
        const FunctionR1R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<1>;

    ////////////////////////////////////////////////////

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceSurface)(
        ModelPart&, const FunctionR2R3::Pointer&, const double&, const double&, const double&, const double&, Properties::Pointer, const int&,
        const FunctionR2R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<0>;

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceSurfaceWithBin)(
        ModelPart&, const FunctionR2R3::Pointer&, const double&, const double&, const double&, const double&, Properties::Pointer, const int&,
        const FunctionR2R3::Pointer&, const bool&, const int&) const = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<1>;

    ////////////////////////////////////////////////////

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceElement)(
        ModelPart&, ModelPart::ElementsContainerType&, Element::Pointer&, Properties::Pointer, const int&, const FunctionR3R3::Pointer&, const bool&, const int&) const
        = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<Element, 0>;

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceElementWithBin)(
        ModelPart&, ModelPart::ElementsContainerType&, Element::Pointer&, Properties::Pointer, const int&, const FunctionR3R3::Pointer&, const bool&, const int&) const
        = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<Element, 1>;

    ////////////////////////////////////////////////////

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceCondition)(
        ModelPart&, ModelPart::ElementsContainerType&, Condition::Pointer&, Properties::Pointer, const int&, const FunctionR3R3::Pointer&, const bool&, const int&) const
        = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<Condition, 0>;

    ModelPart::ConditionsContainerType(ImmersedBoundaryLoadUtility::*pointer_to_SetupImmersedPointForceConditionWithBin)(
        ModelPart&, ModelPart::ElementsContainerType&, Condition::Pointer&, Properties::Pointer, const int&, const FunctionR3R3::Pointer&, const bool&, const int&) const
        = &ImmersedBoundaryLoadUtility::SetupImmersedPointForce<Condition, 1>;

    ////////////////////////////////////////////////////

    class_<ImmersedBoundaryLoadUtility, ImmersedBoundaryLoadUtility::Pointer, bases<ImmersedBoundaryUtility>, boost::noncopyable>
    ("ImmersedBoundaryLoadUtility", init<>())
//    .def("InitializeBinning", &ImmersedBoundaryLoadUtility::InitializeBinning)
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForceCurve)
    .def("SetupImmersedPointForceWithBin", pointer_to_SetupImmersedPointForceCurveWithBin)
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForceSurface)
    .def("SetupImmersedPointForceWithBin", pointer_to_SetupImmersedPointForceSurfaceWithBin)
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForceElement)
    .def("SetupImmersedPointForceWithBin", pointer_to_SetupImmersedPointForceElementWithBin)
    .def("SetupImmersedPointForce", SetupImmersedPointForceElement2<Element, 0>)
    .def("SetupImmersedPointForceWithBin", SetupImmersedPointForceElement2<Element, 1>)
    .def("SetupImmersedPointForce", pointer_to_SetupImmersedPointForceCondition)
    .def("SetupImmersedPointForceWithBin", pointer_to_SetupImmersedPointForceConditionWithBin)
//    .def("GetElements", pointer_to_PyGetElements) // this is moved to FiniteCellAuxilliaryUtility
//    .def("Clean", pointer_to_Clean) // this is moved to FiniteCellAuxilliaryUtility
    ;

}

}  // namespace Python.

}  // namespace Kratos.

