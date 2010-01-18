
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <analysis.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libheliumanalysis)
{
    def("CalculateProjectionRadialProductStates", CalculateProjectionRadialProductStates);
    def("CalculatePopulationRadialProductStates", CalculatePopulationRadialProductStates);
    def("GetWavefunctionParticleExchange",  GetWavefunctionParticleExchange);
    def("AddSingleAngularProjectionAvgPhi",  AddSingleAngularProjectionAvgPhi);
    def("AddDoubleAngularProjectionAvgPhi",  AddDoubleAngularProjectionAvgPhi);
    def("AddDoubleAngularProjectionCoplanar",  AddDoubleAngularProjectionCoplanar);
    def("GetCoulombPhase",  GetCoulombPhase);
    def("RemoveProjectionRadialProductStates",  RemoveProjectionRadialProductStates);
    def("SetRadialCoulombWave",  SetRadialCoulombWave);
}

