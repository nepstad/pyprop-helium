// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_coupledlength_wrapper();
void Export_coupledspherical_wrapper();
void Export_coupledvelocity_wrapper();
void Export_hydrogenmolecule_wrapper();
void Export_superlu_wrapper();
void Export_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libheliumcore)
{
    Export_coupledlength_wrapper();
    Export_coupledspherical_wrapper();
    Export_coupledvelocity_wrapper();
    Export_hydrogenmolecule_wrapper();
    Export_superlu_wrapper();
    Export_wrapper();
}
