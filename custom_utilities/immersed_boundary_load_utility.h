//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            25 Feb 2017
//


#if !defined(KRATOS_IMMERSED_BOUNDARY_LOAD_UTILITY_H_INCLUDED )
#define  KRATOS_IMMERSED_BOUNDARY_LOAD_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/point_3d.h"
#include "brep_application/custom_algebra/function/function.h"
#include "finite_cell_application/custom_utilities/immersed_boundary_utility.h"
#include "custom_conditions/immersed_point_force.h"
#include "custom_utilities/finite_cell_auxiliary_utility.h"
#include "brep_application/brep_application.h"


namespace Kratos
{
///@addtogroup FiniteCellStructuralApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Abstract class for load application on immersed boundary
*/
class ImmersedBoundaryLoadUtility : public ImmersedBoundaryUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImmersedBoundaryLoadUtility
    KRATOS_CLASS_POINTER_DEFINITION(ImmersedBoundaryLoadUtility);

    typedef ImmersedBoundaryUtility BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ImmersedBoundaryLoadUtility() {}

    /// Destructor.
    virtual ~ImmersedBoundaryLoadUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create and add the immersed point force condition to the model_part, providing the curve and the load function
    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            const FunctionR1R3::Pointer& pCurve,        // the parametric curve, there the distributed load is applied
            const double& tmin, const double& tmax,     // the segment where the load is acting
            Properties::Pointer pProperties,
            const int& integration_order,               // integration order
            const FunctionR1R3::Pointer& pLoadFunction, // the distributed load function
            const bool& export_point,
            const int& echo_level
    ) const
    {
        // sampling the points on the interval
        std::vector<double> Coordinates;
        std::vector<PointType> Points;
        std::vector<double> Dlength;
        std::vector<double> Weights;

        ImmersedBoundaryUtility::SamplingCurve(*pCurve, tmin, tmax, integration_order,
                Coordinates, Points, Dlength, Weights);

        // compute the force array
        std::vector<array_1d<double, 3> > Forces(Points.size());
        for(std::size_t point = 0; point < Points.size(); ++point)
        {
//            KRATOS_WATCH(pLoadFunction->GetValue(Coordinates[point]))
            Weights[point] *= Dlength[point];
            Forces[point] = pLoadFunction->GetValue(Coordinates[point]);
        }

        // setup the immersed point force
        return SetupImmersedPointForce<TSearchType>(r_model_part, r_model_part.Elements(), pProperties, Points, Weights, Forces, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing the surface and the load function
    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            const FunctionR2R3::Pointer& pSurface,      // the parametric surface, there the distributed load is applied
            const double& t1min, const double& t1max,   // the segment where the load is acting in 1st-direction
            const double& t2min, const double& t2max,   // the segment where the load is acting in 2nd-direction
            Properties::Pointer pProperties,
            const int& integration_order,               // integration order
            const FunctionR2R3::Pointer& pLoadFunction, // the distributed load function
            const bool& export_point,
            const int& echo_level
    ) const
    {
        // sampling the points on the interval
        std::vector<FunctionR2R3::InputType> Coordinates;
        std::vector<PointType> Points;
        std::vector<double> Darea;
        std::vector<double> Weights;

        ImmersedBoundaryUtility::SamplingSurface(*pSurface, t1min, t1max, t2min, t2max, integration_order,
                Coordinates, Points, Darea, Weights);

        // compute the force array
        std::vector<array_1d<double, 3> > Forces(Points.size());
        for(std::size_t point = 0; point < Points.size(); ++point)
        {
//            KRATOS_WATCH(pLoadFunction->GetValue(Coordinates[point]))
            Weights[point] *= Darea[point];
            Forces[point] = pLoadFunction->GetValue(Coordinates[point]);
        }

        // setup the immersed point force
        return SetupImmersedPointForce<TSearchType>(r_model_part, r_model_part.Elements(), pProperties, Points, Weights, Forces, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing an entity (i.e element/condition)
    template<class TEntityType, int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            typename TEntityType::Pointer& pEntity,
            Properties::Pointer pProperties,
            const int& integration_order,               // integration order
            const FunctionR3R3::Pointer& pLoadFunction, // the distributed load function
            const bool& export_point,
            const int& echo_level
    ) const
    {
        return SetupImmersedPointForce<TEntityType, TSearchType>(r_model_part, r_model_part.Elements(), pEntity, pProperties, integration_order, pLoadFunction, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing an entity (i.e element/condition)
    template<class TEntityType, int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            ModelPart::ElementsContainerType& pMasterElements,
            typename TEntityType::Pointer& pEntity,
            Properties::Pointer pProperties,
            const int& integration_order,               // integration order
            const FunctionR3R3::Pointer& pLoadFunction, // the distributed load function
            const bool& export_point,
            const int& echo_level
    ) const
    {
        // find the maximum condition id
        std::size_t lastCondId = FiniteCellAuxiliaryUtility::GetLastConditionId(r_model_part);

        // find the maximum node Id
        std::size_t lastNodeId = FiniteCellAuxiliaryUtility::GetLastNodeId(r_model_part);

        // find the maximum properties Id
        std::size_t lastPropId = FiniteCellAuxiliaryUtility::GetLastPropertiesId(r_model_part);

        // setup the immersed point force
        return SetupImmersedPointForce<TEntityType, TSearchType>(r_model_part, pMasterElements, pEntity, lastNodeId, lastCondId, pProperties, integration_order, pLoadFunction, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing an entity (i.e element/condition)
    template<class TEntityType, int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            ModelPart::ElementsContainerType& pMasterElements,
            typename TEntityType::Pointer& pEntity,
            std::size_t& lastNodeId,
            std::size_t& lastCondId,
            Properties::Pointer pProperties,
            const int& integration_order,               // integration order
            const FunctionR3R3::Pointer& pLoadFunction, // the distributed load function
            const bool& export_point,
            const int& echo_level
    ) const
    {
        // get the integration point on the geometry
        GeometryData::IntegrationMethod ThisIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        const GeometryType::IntegrationPointsArrayType& integration_points
            = pEntity->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

        std::vector<double> DetJ;
        Function<double, double>::ComputeDetJ(DetJ, pEntity->GetGeometry(), integration_points);

        std::vector<PointType> Points(integration_points.size());
        std::vector<double> Weights(integration_points.size());
        std::vector<array_1d<double, 3> > Forces(integration_points.size());

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            pEntity->GetGeometry().GlobalCoordinates(Points[point], integration_points[point]);
            Weights[point] = DetJ[point] * integration_points[point].Weight();
            Forces[point] = pLoadFunction->GetValue(Points[point]);
            // if(echo_level > 3)
            // {
            //     KRATOS_WATCH(Points[point])
            //     KRATOS_WATCH(pLoadFunction->GetValue(Points[point]))
            //     KRATOS_WATCH(DetJ[point])
            //     KRATOS_WATCH(integration_points[point])
            // }
        }

        // setup the immersed point force
        return SetupImmersedPointForce<TSearchType>(r_model_part, pMasterElements, lastNodeId, lastCondId, pProperties, Points, Weights, Forces, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing the list of points and corresponding weights
    /// r_model_part    the model_part to add the conditions
    /// pMasterElements the pool of master elements to search
    /// rPoints         list of points to create the immersed load
    /// rForces         immersed load value at point
    /// export_point    export the point for post-processing/debugging
    /// echo_level      printing level
    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            ModelPart::ElementsContainerType& pMasterElements,
            Properties::Pointer pProperties,
            const std::vector<PointType>& rPoints,
            const std::vector<double>& rWeights,
            const std::vector<array_1d<double, 3> >& rForces,
            const bool& export_point,
            const int& echo_level) const
    {
        // find the maximum condition id
        std::size_t lastCondId = FiniteCellAuxiliaryUtility::GetLastConditionId(r_model_part);

        // find the maximum node Id
        std::size_t lastNodeId = FiniteCellAuxiliaryUtility::GetLastNodeId(r_model_part);

        // find the maximum properties Id
        std::size_t lastPropId = FiniteCellAuxiliaryUtility::GetLastPropertiesId(r_model_part);

        return SetupImmersedPointForce<TSearchType>(r_model_part, pMasterElements, lastNodeId, lastCondId, pProperties, rPoints, rWeights, rForces, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing the list of points and corresponding weights
    /// r_model_part    the model_part to add the conditions
    /// pMasterElements the pool of master elements to search
    /// lastNodeId      this information is used to assigned the node id incrementally
    /// lastCondId      this information is used to assigned the condition id incrementally
    /// pProperties     the properties used to setup the new condition
    /// rPoints         list of points to create the immersed load
    /// rForces         immersed load value at point
    /// export_point    export the point for post-processing/debugging
    /// echo_level      printing level
    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            ModelPart::ElementsContainerType& pMasterElements,
            std::size_t& lastNodeId,
            std::size_t& lastCondId,
            Properties::Pointer pProperties,
            const std::vector<PointType>& rPoints,
            const std::vector<double>& rWeights,
            const std::vector<array_1d<double, 3> >& rForces,
            const bool& export_point,
            const int& echo_level) const
    {
        // array to contain the new point force conditions
        ModelPart::ConditionsContainerType ImmersedConditions;

        // loop through all the points, find the corresponding containing element
        PointType LocalPoint;
        GeometryType::Pointer pTempGeometry;
        for(std::size_t i = 0; i < rPoints.size(); ++i)
        {
            Element::Pointer pElem;

            bool found;
            bool force_active = true;
            if(TSearchType == 0)
                found = ImmersedBoundaryUtility::SearchPartner(rPoints[i], pMasterElements, pElem, LocalPoint, force_active, echo_level);
            else if(TSearchType == 1)
                found = this->SearchPartnerWithBin(rPoints[i], pMasterElements, pElem, LocalPoint, force_active, echo_level);

            if(found)
            {
                if(pElem->Is(ACTIVE) || (pElem->GetValue(IS_INACTIVE) == false) || (pElem->GetValue(CUT_STATUS) <= 0))
                {
                    if(!export_point)
                    {
                        pTempGeometry = GeometryType::Pointer( new GeometryType() );
                    }
                    else
                    {
                        NodeType::Pointer pNewNode(new NodeType(++lastNodeId, rPoints[i]));
                        pNewNode->SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
                        pNewNode->SetBufferSize(r_model_part.GetBufferSize());
                        r_model_part.AddNode(pNewNode);
                        pTempGeometry = GeometryType::Pointer( new Point3D<NodeType>(pNewNode) );
                    }
//                    Condition::Pointer pNewCond = Condition::Pointer(
//                            new ImmersedPointForce(++lastCondId, pTempGeometry, rWeights[i], rForces[i], pElem, LocalPoint) );
                    Condition::Pointer pNewCond = Condition::Pointer(
                            new ImmersedPointForce(++lastCondId, pTempGeometry, pProperties, rWeights[i], rForces[i], pElem, LocalPoint) );

                    pNewCond->SetValue(IS_INACTIVE, false);
                    pNewCond->Set(ACTIVE, true);

                    ImmersedConditions.push_back(pNewCond);
                    if((echo_level & _REPORT_CONDITION_CREATED_) == _REPORT_CONDITION_CREATED_)
                        std::cout << "Point force is added in element " << pElem->Id() << std::endl;
                }
                else
                {
                    if((echo_level & _WARNING_ELEMENT_IS_INACTIVE_) == _WARNING_ELEMENT_IS_INACTIVE_)
                    {
                        std::cout << "WARNING: Point force at point " << rPoints[i]
                                  << ": the parent element is inactive"
                                  << std::endl;
                    }
                }
            }
            else
            {
                if((echo_level & _WARNING_FOUND_NO_ELEMENT_) == _WARNING_FOUND_NO_ELEMENT_)
                {
                    std::cout << "WARNING: Point force at point " << rPoints[i]
                              << ": no parent element can be found"
                              << std::endl;
                }
            }
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = ImmersedConditions.ptr_begin();
                it != ImmersedConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        if((echo_level & _REPORT_NUMBER_OF_CREATED_CONDITIONS_) == _REPORT_NUMBER_OF_CREATED_CONDITIONS_)
            std::cout << "Setup point forces completed, "
                  << ImmersedConditions.size() << " point force conditions was added to model_part" << std::endl;

        return ImmersedConditions;
    }


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
        return "Immersed Boundary Load Utility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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
    ImmersedBoundaryLoadUtility& operator=(ImmersedBoundaryLoadUtility const& rOther);

    /// Copy constructor.
    ImmersedBoundaryLoadUtility(ImmersedBoundaryLoadUtility const& rOther);


    ///@}

}; // Class ImmersedBoundaryLoadUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ImmersedBoundaryLoadUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ImmersedBoundaryLoadUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef ESTIMATE_RCOND
#undef CHECK_SOLUTION

#endif // KRATOS_IMMERSED_BOUNDARY_LOAD_UTILITY_H_INCLUDED  defined
