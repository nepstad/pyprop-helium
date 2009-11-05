
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <coupledlength.cpp>
#include <coupledspherical.cpp>
#include <coupledvelocity.cpp>
#include <potential.cpp>
#include <superlu.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct CustomPotential_LaserVelocityDerivativeR1_3_Wrapper: CustomPotential_LaserVelocityDerivativeR1<3>
{
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

struct CustomPotentialEvaluationR12_3_Wrapper: CustomPotentialEvaluationR12<3>
{
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
BOOST_PYTHON_MODULE(libheliumcore)
{
    class_< SuperLUSolver<2>, boost::noncopyable >("SuperLUSolver_2", init<  >())
        .def("Setup", &SuperLUSolver<2>::Setup)
        .def("Solve", &SuperLUSolver<2>::Solve)
        .def("PrintStatistics", &SuperLUSolver<2>::PrintStatistics)
    ;

    class_< CustomPotential_LaserVelocityDerivativeR1<3>, boost::noncopyable, CustomPotential_LaserVelocityDerivativeR1_3_Wrapper >("CustomPotential_LaserVelocityDerivativeR1_3", init<  >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocityDerivativeR1<3>::SetBasisPairs, &CustomPotential_LaserVelocityDerivativeR1_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocityDerivativeR1<3>::UpdatePotentialData, &CustomPotential_LaserVelocityDerivativeR1_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocityDerivativeR1<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocityDerivativeR1<3>::GetBasisPairList)
    ;

    class_< CustomPotential_LaserVelocityDerivativeR2<3>, boost::noncopyable, CustomPotential_LaserVelocityDerivativeR2_3_Wrapper >("CustomPotential_LaserVelocityDerivativeR2_3", init<  >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocityDerivativeR2<3>::SetBasisPairs, &CustomPotential_LaserVelocityDerivativeR2_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocityDerivativeR2<3>::UpdatePotentialData, &CustomPotential_LaserVelocityDerivativeR2_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocityDerivativeR2<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocityDerivativeR2<3>::GetBasisPairList)
    ;

    class_< CustomPotential_LaserVelocity<3>, boost::noncopyable, CustomPotential_LaserVelocity_3_Wrapper >("CustomPotential_LaserVelocity_3", init<  >())
        .def("SetBasisPairs", &CustomPotential_LaserVelocity<3>::SetBasisPairs, &CustomPotential_LaserVelocity_3_Wrapper::default_SetBasisPairs)
        .def("UpdatePotentialData", &CustomPotential_LaserVelocity<3>::UpdatePotentialData, &CustomPotential_LaserVelocity_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", &CustomPotential_LaserVelocity<3>::ApplyConfigSection)
        .def("GetBasisPairList", &CustomPotential_LaserVelocity<3>::GetBasisPairList)
    ;

    class_< CustomPotential_LaserLength<3>, boost::noncopyable, CustomPotential_LaserLength_3_Wrapper >("CustomPotential_LaserLength_3", init<  >())
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotential_LaserLength<3>::UpdatePotentialData, (void (CustomPotential_LaserLength_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotential_LaserLength_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", (void (CustomPotentialCoupledSphericalBase<3>::*)(const ConfigSection&) )&CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection, (void (CustomPotential_LaserLength_3_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_3_Wrapper::default_ApplyConfigSection)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotential_LaserLength_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_3_Wrapper::*)(int))&CustomPotential_LaserLength_3_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength<3>::Coefficient)
        .def("CondonShortleyPhase", &CustomPotential_LaserLength<3>::CondonShortleyPhase)
        .staticmethod("Coefficient")
        .staticmethod("CondonShortleyPhase")
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>, boost::noncopyable >("KineticEnergyPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>, boost::noncopyable >("CoupledSphericalCoulombPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CoupledSphericalCoulombPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< CustomPotentialEvaluationR12<3>, boost::noncopyable, CustomPotentialEvaluationR12_3_Wrapper >("CustomPotentialEvaluationR12_3", init<  >())
        .def("UpdatePotentialData", (void (CustomPotentialEvaluationR12<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CustomPotentialEvaluationR12<3>::UpdatePotentialData, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CustomPotentialEvaluationR12_3_Wrapper::default_UpdatePotentialData)
        .def("ApplyConfigSection", (void (CustomPotentialCoupledSphericalBase<3>::*)(const ConfigSection&) )&CustomPotentialCoupledSphericalBase<3>::ApplyConfigSection, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(const ConfigSection&))&CustomPotentialEvaluationR12_3_Wrapper::default_ApplyConfigSection)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CustomPotentialEvaluationR12_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotentialEvaluationR12_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CustomPotentialEvaluationR12_3_Wrapper::*)(int))&CustomPotentialEvaluationR12_3_Wrapper::default_GetBasisPairList)
        .def("CoefficientR12", &CustomPotentialEvaluationR12<3>::CoefficientR12)
        .def("CondonShortleyPhase", &CustomPotentialEvaluationR12<3>::CondonShortleyPhase)
        .staticmethod("CondonShortleyPhase")
        .staticmethod("CoefficientR12")
    ;

    class_< CoupledSphericalKineticEnergyEvaluator<3>, boost::noncopyable, CoupledSphericalKineticEnergyEvaluator_3_Wrapper >("CoupledSphericalKineticEnergyEvaluator_3", init<  >())
        .def_readwrite("Mass", &CoupledSphericalKineticEnergyEvaluator<3>::Mass)
        .def("ApplyConfigSection", (void (CoupledSphericalKineticEnergyEvaluator<3>::*)(const ConfigSection&) )&CoupledSphericalKineticEnergyEvaluator<3>::ApplyConfigSection, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(const ConfigSection&))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CoupledSphericalKineticEnergyEvaluator<3>::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double) )&CoupledSphericalKineticEnergyEvaluator<3>::UpdatePotentialData, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(blitz::Array<std::complex<double>,3>, boost::shared_ptr<Wavefunction<3> >, std::complex<double>, double))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialCoupledSphericalBase<3>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialCoupledSphericalBase<3>::SetBasisPairs, (void (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(int, const blitz::Array<int,2>&))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialCoupledSphericalBase<3>::*)(int) )&CustomPotentialCoupledSphericalBase<3>::GetBasisPairList, (blitz::Array<int,2> (CoupledSphericalKineticEnergyEvaluator_3_Wrapper::*)(int))&CoupledSphericalKineticEnergyEvaluator_3_Wrapper::default_GetBasisPairList)
    ;

    class_< DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>, boost::noncopyable >("ComplexAbsorbingPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<OverlapPotential<3>,3>, boost::noncopyable >("OverlapPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<OverlapPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<BoxNormPotential<3>,3>, boost::noncopyable >("BoxNormPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<BoxNormPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<SingleIonizationBox<3>,3>, boost::noncopyable >("SingleIonizationBox_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<SingleIonizationBox<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>, boost::noncopyable >("DoubleIonizationBox_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<DoubleIonizationBox<3>,3>::CalculateExpectationValue)
    ;

}

