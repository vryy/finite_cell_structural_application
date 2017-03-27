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
#include "finite_cell_application/custom_algebra/function/function.h"
#include "finite_cell_application/custom_utilities/immersed_boundary_utility.h"
#include "custom_conditions/immersed_point_force.h"
#include "finite_cell_application.h"


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

    typedef Function<double, PointType> FunctionR1R3;

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


    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            const FunctionR1R3::Pointer& pCurve,    // the parametric curve, there the distributed load is applied
            const double& tmin, const double& tmax, // the segment where the load is acting
            const int& integration_order,           // integration order
            const FunctionR1R3::Pointer& pLoadFunction,  // the distributed load function
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
            Forces[point] = Dlength[point] * Weights[point] * pLoadFunction->GetValue(Coordinates[point]);
        }

        // setup the immersed point force
        return SetupImmersedPointForce<TSearchType>(r_model_part, Points, Forces, export_point, echo_level);
    }


    /// Create and add the immersed point force condition to the model_part, providing the list of points and corresponding weights
    template<int TSearchType>
    ModelPart::ConditionsContainerType SetupImmersedPointForce(ModelPart& r_model_part,
            const std::vector<PointType>& rPoints,
            const std::vector<array_1d<double, 3> >& rForces,
            const bool& export_point,
            const int& echo_level) const
    {
        // find the maximum condition id
        std::size_t lastCondId = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            if((*it)->Id() > lastCondId)
                lastCondId = (*it)->Id();
        }

        // find the maximum node Id
        std::size_t lastNodeId = 0;
        for(typename ModelPart::NodesContainerType::ptr_iterator it = r_model_part.Nodes().ptr_begin();
                it != r_model_part.Nodes().ptr_end(); ++it)
        {
            if((*it)->Id() > lastNodeId)
                lastNodeId = (*it)->Id();
        }

        // find the maximum properties Id
        std::size_t lastPropId = 0;
        for(typename ModelPart::PropertiesContainerType::ptr_iterator it = r_model_part.rProperties().ptr_begin();
                it != r_model_part.rProperties().ptr_end(); ++it)
        {
            if((*it)->Id() > lastPropId)
                lastPropId = (*it)->Id();
        }

        // create new properties
        Properties::Pointer pProperties = Properties::Pointer(new Properties(++lastPropId));
        r_model_part.AddProperties(pProperties);

        // array to contain the new point force conditions
        ModelPart::ConditionsContainerType ImmersedConditions;

        // loop through all the points, find the corresponding containing element
        PointType LocalPoint;
        GeometryType::Pointer pTempGeometry;
        for(std::size_t i = 0; i < rPoints.size(); ++i)
        {
            Element::Pointer pElem;

            bool found;
            if(TSearchType == 0)
                found = ImmersedBoundaryUtility::SearchPartner(rPoints[i], r_model_part.Elements(), pElem, LocalPoint);
            else if(TSearchType == 1)
                found = this->SearchPartnerWithBin(rPoints[i], r_model_part.Elements(), pElem, LocalPoint);

            if(!found)
                std::cout << "WARNING: can't find the containing element for point " << rPoints[i] << std::endl;
            else
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
                    Condition::Pointer pNewCond = Condition::Pointer(
                            new ImmersedPointForce(++lastCondId, pTempGeometry, rForces[i], pElem, LocalPoint) );
//                    Condition::Pointer pNewCond = Condition::Pointer(
//                            new ImmersedPointForce(++lastCondId, pTempGeometry, pProperties, rForces[i], pElem, LocalPoint) );

                    pNewCond->Set(ACTIVE, true);

                    ImmersedConditions.push_back(pNewCond);
                    if(echo_level > 1)
                        std::cout << "Point force is added in element " << pElem->Id() << std::endl;
                }
            }
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = ImmersedConditions.ptr_begin();
                it != ImmersedConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        if(echo_level > 0)
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
