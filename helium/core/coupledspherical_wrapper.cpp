
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <coupledspherical.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct CustomPotentialEvaluationR12_3_Wrapper: CustomPotentialEvaluationR12<3>
{
    CustomPotentialEvaluationR12_3_Wrapper(PyObject* py_self_, const CustomPotentialEvaluationR12<3>& p0):
        CustomPotentialEvaluationR12<3>(p0), py_self(py_self_) {}

    CustomPotentialEvaluationR12_3_Wrapper(PyObject* py_self_):
        CustomPotentialEvaluationR12<3>(), py_self(py_self_) {}

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CustomPotentialEvaluationR12<3>::UpdatePotentialData(p0, p1, p2, p3);
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

struct CoupledSphericalKineticEnergyEvaluator_3_Wrapper: CoupledSphericalKineticEnergyEvaluator<3>
{
    CoupledSphericalKineticEnergyEvaluator_3_Wrapper(PyObject* py_self_, const CoupledSphericalKineticEnergyEvaluator<3>& p0):
        CoupledSphericalKineticEnergyEvaluator<3>(p0), py_self(py_self_) {}

    CoupledSphericalKineticEnergyEvaluator_3_Wrapper(PyObject* py_self_):
        CoupledSphericalKineticEnergyEvaluator<3>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CoupledSphericalKineticEnergyEvaluator<3>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        CoupledSphericalKineticEnergyEvaluator<3>::UpdatePotentialData(p0, p1, p2, p3);
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
void Export_coupledspherical_wrapper()
{
    class_< CustomPotentialEvaluationR12<3>, CustomPotentialEvaluationR12_3_Wrapper >("CustomPotentialEvaluationR12_3", init<  >())
        .def(init< const CustomPotentialEvaluationR12<3>& >())
        .def("UpdatePotentialData", (void (CustomPotentialEvaluationR12<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotentialEvaluationR12<3>::UpdatePotentialData, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotentialEvaluationR12_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", (void (CustomPotentialCoupledSphericalBase<3>::*)(const ConfigSection&) )&CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(const ConfigSection&))&CustomPotentialEvaluationR12_3_Wrapper::default_ApplyConfigSection)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotentialEvaluationR12_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotentialEvaluationR12_3_Wrapper::*)(int))&CustomPotentialEvaluationR12_3_Wrapper::default_GetBasisPairList)
        .def("CoefficientR12", &CustomPotentialEvaluationR12<3>::CoefficientR12)
        .def("CondonShortleyPhase", &CustomPotentialEvaluationR12<3>::CondonShortleyPhase)
        .staticmethod("CondonShortleyPhase")
        .staticmethod("CoefficientR12")
    ;

    class_< CoupledSphericalKineticEnergyEvaluator<3>, CoupledSphericalKineticEnergyEvaluator_3_Wrapper >("CoupledSphericalKineticEnergyEvaluator_3", init<  >())
        .def(init< const CoupledSphericalKineticEnergyEvaluator<3>& >())
        .def_readwrite("Mass", &CoupledSphericalKineticEnergyEvaluator<3>::Mass)
        .def("ApplyConfigSection", (void (CoupledSphericalKineticEnergyEvaluator<3>::*)(const ConfigSection&) )&CoupledSphericalKineticEnergyEvaluator<3>::ApplyConfigSection, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(const ConfigSection&))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CoupledSphericalKineticEnergyEvaluator<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CoupledSphericalKineticEnergyEvaluator<3>::UpdatePotentialData, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(int))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_GetBasisPairList)
    ;

}

