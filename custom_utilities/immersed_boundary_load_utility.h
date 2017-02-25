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
#include "finite_cell_application/custom_algebra/function.h"
#include "finite_cell_application/custom_utilities/immersed_boundary_utility.h"
#include "custom_conditions/immersed_point_force.h"


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
class ImmersedBoundaryLoadUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImmersedBoundaryLoadUtility
    KRATOS_CLASS_POINTER_DEFINITION(ImmersedBoundaryLoadUtility);

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


    void SetupImmersedPointForce(ModelPart& r_model_part, const FunctionR1R3& rCurve,
            const double& tmin, const double& tmax, const int& integration_order,
            const FunctionR1R3& rLoadFunction) const
    {
        // sampling the points on the interval
        std::vector<double> Coordinates;
        std::vector<PointType> Points;
        std::vector<double> Dlength;
        std::vector<double> Weights;
        ImmersedBoundaryUtility::SamplingCurve(rCurve, tmin, tmax, integration_order,
                Coordinates, Points, Dlength, Weights);

        // compute the force array
        std::vector<array_1d<double, 3> > Forces(Points.size());
        for(std::size_t point = 0; point < Points.size(); ++point)
        {
            Forces[point] = Dlength[point] * Weights[point] * rLoadFunction.GetValue(Coordinates[point]);
        }

        // setup the immersed point force
        SetupImmersedPointForce(r_model_part, Points, Forces);
    }


    /// Create and add the immersed point force condition to the model_part, providing the list of points and corresponding weights
    static void SetupImmersedPointForce(ModelPart& r_model_part,
            const std::vector<PointType>& rPoints, const std::vector<array_1d<double, 3> >& rForces)
    {
        // find the maximum condition id
        std::size_t lastCondId = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            // find the maximum condition id
            if((*it)->Id() > lastCondId)
                lastCondId = (*it)->Id();
        }

        // array to contain the new point force conditions
        ModelPart::ConditionsContainerType NewConditions;

        // loop through all the points, find the corresponding containing element
        PointType LocalPoint;
        GeometryType::Pointer pTempGeometry;
        for(std::size_t i = 0; i < rPoints.size(); ++i)
        {
            Element::Pointer pElem;
            bool found = ImmersedBoundaryUtility::SearchPartner(rPoints[i], r_model_part.Elements(), pElem, LocalPoint);
            if(!found)
                std::cout << "WARNING: can't find the containing element for point " << rPoints[i] << std::endl;
            else
            {
                pTempGeometry = GeometryType::Pointer( new GeometryType() );
                Condition::Pointer pNewCond = Condition::Pointer(
                        new ImmersedPointForce(++lastCondId, pTempGeometry, rForces[i], pElem, LocalPoint) );

                pNewCond->Set(ACTIVE, true);

                NewConditions.push_back(pNewCond);
            }
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin();
                it != NewConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        std::cout << "Setup point forces completed, "
              << NewConditions.size() << " point force conditions was added to model_part" << std::endl;
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
{}

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
