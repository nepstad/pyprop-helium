
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_wrapper()
{
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

