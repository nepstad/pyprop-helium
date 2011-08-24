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


/*
 * Manolopoulos transmission free absorbing potential
 *
 *   See D. E. Manolopoulos, J. Chem. Phys 117, 9552, 2002.
 *
 *   energy_cutoff: wavepacket components with energy greater than this is
 *                  absorbed
 *   grid_max:      last point on the grid absorberwise
 *   delta:         accuracy parameter, determines width of absorber
 *
 */
template<int Rank>
class ManolopoulosAbsorber : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Rank info
	int RadialRank1, RadialRank2;

	//Potential parameters
	double EnergyCutoff;
	double GridMax;
	double Delta;
	double Start;

	//CAP constants
	double A, B, C;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("energy_cutoff", EnergyCutoff);
		config.Get("grid_max", GridMax);
		config.Get("absorber_start", Start);
		config.Get("radial_rank_1", RadialRank1);
		config.Get("radial_rank_2", RadialRank2);

		A = 0.112449;
		B = 0.0082735;
		C = 2.62206;

		//Calculate absorber start
		double kmin = std::sqrt(2*EnergyCutoff);
		Start = GridMax - C / (2 * Delta * kmin);

		//Calculate delta
		double kmin = std::sqrt(2*EnergyCutoff);
		Delta = C / (2 * kmin * (GridMax - Start));
		
		cout << "Manolopoulos absorber" << endl;
		cout << "  Absorber starts at " << Start << std::endl;
		cout << "  Delta parameter = " << Delta << std::endl;		
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(RadialRank1);
		double r2 = pos(RadialRank2);
		cplx cap = 0;

		// Calculate absorber for r1
		if (r1 >= Start)
		{
			double u = 2 * Delta * std::sqrt(2*EnergyCutoff) * (r1 - Start);
			double y = A*u - B*u*u*u + 4.0/((C-u)*(C-u)) - 4.0/((C+u)*(C+u));
			cap += -I * EnergyCutoff * y;
		}

		// Calculate absorber for r2
		if (r2 >= Start)
		{
			double u = 2 * Delta * std::sqrt(2*EnergyCutoff) * (r2 - Start);
			double y = A*u - B*u*u*u + 4.0/((C-u)*(C-u)) - 4.0/((C+u)*(C+u));
			cap += -I * EnergyCutoff * y;
		}

		return cap;
	}
};

