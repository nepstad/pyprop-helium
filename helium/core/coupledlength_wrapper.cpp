
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <coupledlength.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct CustomPotential_LaserLength_3_Wrapper: CustomPotential_LaserLength<3>
{
    CustomPotential_LaserLength_3_Wrapper(PyObject* py_self_, const CustomPotential_LaserLength<3>& p0):
        CustomPotential_LaserLength<3>(p0), py_self(py_self_) {}

    CustomPotential_LaserLength_3_Wrapper(PyObject* py_self_):
        CustomPotential_LaserLength<3>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotential_LaserLength<3>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength<3>::UpdatePotentialData(p0, p1, p2, p3);
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

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotential_LaserLength_X<3>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength_X<3>::UpdatePotentialData(p0, p1, p2, p3);
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
    class_< CustomPotential_LaserLength<3>, CustomPotential_LaserLength_3_Wrapper >("CustomPotential_LaserLength_3", init<  >())
        .def(init< const CustomPotential_LaserLength<3>& >())
        .def_readwrite("Scaling", &CustomPotential_LaserLength<3>::Scaling)
        .def("ApplyConfigSection", (void (CustomPotential_LaserLength<3>::*)(const ConfigSection&) )&CustomPotential_LaserLength<3>::ApplyConfigSection, (void (CustomPotential_LaserLength_3_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_3_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotential_LaserLength<3>::UpdatePotentialData, (void (CustomPotential_LaserLength_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotential_LaserLength_3_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotential_LaserLength_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_3_Wrapper::*)(int))&CustomPotential_LaserLength_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength<3>::Coefficient)
        .staticmethod("Coefficient")
    ;

    class_< CustomPotential_LaserLength_X<3>, CustomPotential_LaserLength_X_3_Wrapper >("CustomPotential_LaserLength_X_3", init<  >())
        .def(init< const CustomPotential_LaserLength_X<3>& >())
        .def_readwrite("Scaling", &CustomPotential_LaserLength_X<3>::Scaling)
        .def("ApplyConfigSection", (void (CustomPotential_LaserLength_X<3>::*)(const ConfigSection&) )&CustomPotential_LaserLength_X<3>::ApplyConfigSection, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_X_3_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength_X<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotential_LaserLength_X<3>::UpdatePotentialData, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotential_LaserLength_X_3_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotential_LaserLength_X_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_X_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_X_3_Wrapper::*)(int))&CustomPotential_LaserLength_X_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength_X<3>::Coefficient)
        .staticmethod("Coefficient")
    ;

}

