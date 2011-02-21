
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <diatomicpotential.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct DiatomicCoulomb_3_Wrapper: DiatomicCoulomb<3>
{
    DiatomicCoulomb_3_Wrapper(PyObject* py_self_, const DiatomicCoulomb<3>& p0):
        DiatomicCoulomb<3>(p0), py_self(py_self_) {}

    DiatomicCoulomb_3_Wrapper(PyObject* py_self_):
        DiatomicCoulomb<3>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        DiatomicCoulomb<3>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,3> p0, boost::shared_ptr<Wavefunction<3> > p1, std::complex<double> p2, double p3) {
        DiatomicCoulomb<3>::UpdatePotentialData(p0, p1, p2, p3);
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
void Export_diatomicpotential_wrapper()
{
    class_< DiatomicCoulomb<3>, DiatomicCoulomb_3_Wrapper >("DiatomicCoulomb_3", init<  >())
        .def(init< const DiatomicCoulomb<3>& >())
        .def_readwrite("R", &DiatomicCoulomb<3>::R)
        .def_readwrite("ThetaR", &DiatomicCoulomb<3>::ThetaR)
        .def_readwrite("pluss", &DiatomicCoulomb<3>::pluss)
        .def("ApplyConfigSection", (void (DiatomicCoulomb<3>::*)(const ConfigSection&) )&DiatomicCoulomb<3>::ApplyConfigSection, (void (DiatomicCoulomb_3_Wrapper::*)(const ConfigSection&))&DiatomicCoulomb_3_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (DiatomicCoulomb<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&DiatomicCoulomb<3>::UpdatePotentialData, (void (DiatomicCoulomb_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&DiatomicCoulomb_3_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (DiatomicCoulomb_3_Wrapper::*)(int, const blitz::Array<int,2>&))&DiatomicCoulomb_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (DiatomicCoulomb_3_Wrapper::*)(int))&DiatomicCoulomb_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &DiatomicCoulomb<3>::Coefficient)
        .def("MultipoleCoeff", &DiatomicCoulomb<3>::MultipoleCoeff)
        .def("CondonShortleyPhase", &DiatomicCoulomb<3>::CondonShortleyPhase)
        .staticmethod("Coefficient")
        .staticmethod("MultipoleCoeff")
        .staticmethod("CondonShortleyPhase")
    ;

}

