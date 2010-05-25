
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <coupledlength.cpp>
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

struct CustomPotential_LaserLength_3_Wrapper: CustomPotential_LaserLength<3>
{
    CustomPotential_LaserLength_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserLength<3>& p0):
        CustomPotential_LaserLength<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserLength_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserLength<3>(), py_self(py_self_) {}

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength<3>::UpdatePotentialData(p0, p1, p2, p3);
    }

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection(p0);
    }

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotentialCoupledSphericalBase<3>::SetBasisPairs(p0, p1);
    }

    blitz::Array<int,2> GetBasisPairList(int p0) {
        return call_method< blitz::Array<int,2> >(py_self, "GetBasisPairList", p0);
    }

    blitz::Array<int,2> default_GetBasisPairList(int p0) {
        return CustomPotentialCoupledSphericalBase<3>::GetBasisPairList(p0);
    }

    PyObject* py_self;
};

struct CustomPotential_LaserLength_X_3_Wrapper: CustomPotential_LaserLength_X<3>
{
    CustomPotential_LaserLength_X_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserLength_X<3>& p0):
        CustomPotential_LaserLength_X<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserLength_X_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserLength_X<3>(), py_self(py_self_) {}

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength_X<3>::UpdatePotentialData(p0, p1, p2, p3);
    }

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection(p0);
    }

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotentialCoupledSphericalBase<3>::SetBasisPairs(p0, p1);
    }

    blitz::Array<int,2> GetBasisPairList(int p0) {
        return call_method< blitz::Array<int,2> >(py_self, "GetBasisPairList", p0);
    }

    blitz::Array<int,2> default_GetBasisPairList(int p0) {
        return CustomPotentialCoupledSphericalBase<3>::GetBasisPairList(p0);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_coupledlength_wrapper()
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

    class_< CustomPotential_LaserLength<3>, CustomPotential_LaserLength_3_Wrapper >("CustomPotential_LaserLength_3", init<  >())
        .def(init< const CustomPotential_LaserLength<3>& >())
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotential_LaserLength<3>::UpdatePotentialData, (void (CustomPotential_LaserLength_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotential_LaserLength_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", (void (CustomPotentialCoupledSphericalBase<3>::*)(const ConfigSection&) )&CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection, (void (CustomPotential_LaserLength_3_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_3_Wrapper::default_ApplyConfigSection)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotential_LaserLength_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_3_Wrapper::*)(int))&CustomPotential_LaserLength_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength<3>::Coefficient)
        .staticmethod("Coefficient")
    ;

    class_< CustomPotential_LaserLength_X<3>, CustomPotential_LaserLength_X_3_Wrapper >("CustomPotential_LaserLength_X_3", init<  >())
        .def(init< const CustomPotential_LaserLength_X<3>& >())
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength_X<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotential_LaserLength_X<3>::UpdatePotentialData, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotential_LaserLength_X_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", (void (CustomPotentialCoupledSphericalBase<3>::*)(const ConfigSection&) )&CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_X_3_Wrapper::default_ApplyConfigSection)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_X_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_X_3_Wrapper::*)(int))&CustomPotential_LaserLength_X_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength_X<3>::Coefficient)
        .staticmethod("Coefficient")
    ;

}

