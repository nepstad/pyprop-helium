
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <superlu.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_superlu_wrapper()
{
    class_< SuperLUSolver<2>, boost::noncopyable >("SuperLUSolver_2", init<  >())
        .def("Setup", &SuperLUSolver<2>::Setup)
        .def("Solve", &SuperLUSolver<2>::Solve)
        .def("PrintStatistics", &SuperLUSolver<2>::PrintStatistics)
    ;

}

