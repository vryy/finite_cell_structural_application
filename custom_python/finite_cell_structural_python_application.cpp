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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "finite_cell_structural_application.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosFiniteCellStructuralApplication)
    {

        class_<KratosFiniteCellStructuralApplication, KratosFiniteCellStructuralApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosFiniteCellStructuralApplication");

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
