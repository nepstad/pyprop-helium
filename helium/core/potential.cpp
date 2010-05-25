#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

/*
 * Radial Kinetic Energy
 */
template<int Rank>
class KineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return - 1. / (2. * mass);
	}

};

/* 
 * Coulomb Potential
 */
template<int Rank>
class CoupledSphericalCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Z;
	int RadialRank1;
	int RadialRank2;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(RadialRank1);
		double r2 = pos(RadialRank2);
		return - Z/r1 - Z/r2;
	}
};


template<int Rank>
class ComplexAbsorbingPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int radialRank1;
	int radialRank2;
	double scalingReal;
	double scalingImag;
	double factorReal;
	double factorImag;
	double absorberStart;
	double absorberLength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank_1", radialRank1);
		// Stian og Sigurd har vaert her
		config.Get("radial_rank_2", radialRank2);
		config.Get("absorber_start", absorberStart);
		config.Get("absorber_length", absorberLength);
		config.Get("scaling_real", scalingReal);
		config.Get("scaling_imag", scalingImag);
		config.Get("factor_real", factorReal);
		config.Get("factor_imag", factorImag);
	}

	/*
	 * Called for every grid point 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(radialRank1);
		double r2 = pos(radialRank2);
		cplx V = 0;
		if (r1 > absorberStart)
		{
			double curLength1 = (r1 - absorberStart) / absorberLength;
			double Vr1 = factorReal * std::pow(curLength1, scalingReal);
			double Vi1 = factorImag * std::pow(curLength1, scalingImag);
			V += cplx(Vr1 , Vi1);
		}
		if (r2 > absorberStart)
		{
			double curLength2 = (r2 - absorberStart) / absorberLength;
			double Vr2 = factorReal * std::pow(curLength2, scalingReal);
			double Vi2 = factorImag * std::pow(curLength2, scalingImag);
			V += cplx(Vr2 , Vi2);
		}
		return V;
	}
};


template<int Rank>
class OverlapPotential : public PotentialBase<Rank>
{
public:
        //Required by DynamicPotentialEvaluator
        cplx TimeStep;
        double CurTime;


        void ApplyConfigSection(const ConfigSection &config)
        {
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
                return 1.0;
        }
};


template<int Rank>
class BoxNormPotential : public PotentialBase<Rank>
{
public:
		int RadialRank1;
		int RadialRank2;
		int BoxSize;

        void ApplyConfigSection(const ConfigSection &config)
        {
			config.Get("radial_rank_1", RadialRank1);
			config.Get("radial_rank_2", RadialRank2);
			config.Get("box_size", BoxSize);
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
			double r1 = pos(RadialRank1);
			double r2 = pos(RadialRank2);
			double r = std::sqrt(r1 * r1 + r2 * r2);

			if (r < BoxSize)
				return 1.0;
			else 
				return 0.0;
        }
};

template<int Rank>
class SingleIonizationBox : public PotentialBase<Rank>
{
public:
		int RadialRank1;
		int RadialRank2;
		double InnerBoxSize;
		double Width;

        void ApplyConfigSection(const ConfigSection &config)
        {
			config.Get("radial_rank_1", RadialRank1);
			config.Get("radial_rank_2", RadialRank2);
			config.Get("inner_box_size", InnerBoxSize);
			config.Get("width", Width);
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
			double r1 = pos(RadialRank1);
			double r2 = pos(RadialRank2);

			if (r1 > InnerBoxSize && r2 < Width)
				return 1.0;
			else if (r2 > InnerBoxSize and r1 < Width)
				return 1.0;
			else 
				return 0.0;
        }
};

template<int Rank>
class DoubleIonizationBox : public PotentialBase<Rank>
{
public:
		int RadialRank1;
		int RadialRank2;
		double InnerBoxSize;
		double Width;

        void ApplyConfigSection(const ConfigSection &config)
        {
			config.Get("radial_rank_1", RadialRank1);
			config.Get("radial_rank_2", RadialRank2);
			config.Get("inner_box_size", InnerBoxSize);
			config.Get("width", Width);
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
			double r1 = pos(RadialRank1);
			double r2 = pos(RadialRank2);

			if (r1 > InnerBoxSize && r2 > InnerBoxSize)
				return 1.0;
			else 
				return 0.0;
        }
};

