
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <coupledvelocity.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct CustomPotential_LaserVelocityDerivativeR1_3_Wrapper: CustomPotential_LaserVelocityDerivativeR1<3>
{
    CustomPotential_LaserVelocityDerivativeR1_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserVelocityDerivativeR1<3>& p0):
        CustomPotential_LaserVelocityDerivativeR1<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserVelocityDerivativeR1_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserVelocityDerivativeR1<3>(), py_self(py_self_) {}

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotential_LaserVelocityDerivativeR1<3>::SetBasisPairs(p0, p1);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserVelocityDerivativeR1<3>::UpdatePotentialData(p0, p1, p2, p3);
    }

    PyObject* py_self;
};

struct CustomPotential_LaserVelocityDerivativeR2_3_Wrapper: CustomPotential_LaserVelocityDerivativeR2<3>
{
    CustomPotential_LaserVelocityDerivativeR2_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserVelocityDerivativeR2<3>& p0):
        CustomPotential_LaserVelocityDerivativeR2<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserVelocityDerivativeR2_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserVelocityDerivativeR2<3>(), py_self(py_self_) {}

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotential_LaserVelocityDerivativeR2<3>::SetBasisPairs(p0, p1);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserVelocityDerivativeR2<3>::UpdatePotentialData(p0, p1, p2, p3);
    }

    PyObject* py_self;
};

struct CustomPotential_LaserVelocity_3_Wrapper: CustomPotential_LaserVelocity<3>
{
    CustomPotential_LaserVelocity_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserVelocity<3>& p0):
        CustomPotential_LaserVelocity<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserVelocity_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserVelocity<3>(), py_self(py_self_) {}

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotential_LaserVelocity<3>::SetBasisPairs(p0, p1);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserVelocity<3>::UpdatePotentialData(p0, p1, p2, p3);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_coupledvelocity_wrapper()
{
    class_< CustomPotential_LaserVelocityDerivativeR1<3>, CustomPotential_LaserVelocityDerivativeR1_3_Wrapper >("CustomPotential_LaserVelocityDerivativeR1_3", init<  >())
        .def(init< const CustomPotential_LaserVelocityDerivativeR1<3>& >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocityDerivativeR1<3>::SetBasisPairs, &CustomPotential_LaserVelocityDerivativeR1_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocityDerivativeR1<3>::UpdatePotentialData, &CustomPotential_LaserVelocityDerivativeR1_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocityDerivativeR1<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocityDerivativeR1<3>::GetBasisPairList)
    ;

    class_< CustomPotential_LaserVelocityDerivativeR2<3>, CustomPotential_LaserVelocityDerivativeR2_3_Wrapper >("CustomPotential_LaserVelocityDerivativeR2_3", init<  >())
        .def(init< const CustomPotential_LaserVelocityDerivativeR2<3>& >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocityDerivativeR2<3>::SetBasisPairs, &CustomPotential_LaserVelocityDerivativeR2_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocityDerivativeR2<3>::UpdatePotentialData, &CustomPotential_LaserVelocityDerivativeR2_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocityDerivativeR2<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocityDerivativeR2<3>::GetBasisPairList)
    ;

    class_< CustomPotential_LaserVelocity<3>, CustomPotential_LaserVelocity_3_Wrapper >("CustomPotential_LaserVelocity_3", init<  >())
        .def(init< const CustomPotential_LaserVelocity<3>& >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocity<3>::SetBasisPairs, &CustomPotential_LaserVelocity_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocity<3>::UpdatePotentialData, &CustomPotential_LaserVelocity_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocity<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocity<3>::GetBasisPairList)
    ;

}

