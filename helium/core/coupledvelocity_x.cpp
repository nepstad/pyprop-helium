#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/representation/combinedrepresentation.h>
#include <core/representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

#include "laserhelper.h"
#include <gsl/gsl_sf_gamma.h>
using namespace CoupledSpherical;



/* First part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics
 *
 * <Ylm | - \frac{1}{r} \frac{\sin \phi}{\sin \theta} \partialdiff{}{\phi}
 *  	  + \frac{1}{r} \cos \phi \cos \theta \partialdiff{}{\theta}
 *        - \frac{1}{r} \cos \phi \sin \theta | Yl'm'>
 *
 */
template<int Rank>
class CustomPotential_LaserVelocity_X
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocity_X() {}
	virtual ~CustomPotential_LaserVelocity_X() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);
	

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers 
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			/*
			 * A selection rule based on the fact that \Delta m = \pm 1 and 
			 * the kronecter deltas in the matrix element <Y_L'^M'|H(r_1) + H(r_2) |Y_L^M>.
			 */
			if (std::abs(Mp - M) != 1) continue;

			double coupling1 = 0;
			double coupling2 = 0;
			double I1_1 = 0;
			double I1_2 = 0;
			double I2_1 = 0; 
			double I2_2 = 0;
			double I3_1 = 0;
			double I3_2 = 0;
			double dlta1 = 0;
			double dlta2 = 0;
			double norms = 0;
			double temp;
			double E_lm;
			double F_lm;
			double G_lm;
			double H_lm;
			double I_lm;
			double J_lm;
			double K_lm;
			double dlta_m;
			double K2_term1;
			double K2_term2;
			double K1_term1;
			double K1_term2;
			double K1_term3;

			double dlta_l1 = LaserHelper::kronecker(l1p,l1);
			double dlta_l2 = LaserHelper::kronecker(l2p,l2);


			for (int m1=-l1; m1<=l1; m1++) 
			{
				for (int m1p = (m1 - 1); m1p <= (m1 + 1); m1p++)
				{
					if (m1p != m1)
					{
						int m2 = M - m1;
						int m2p = Mp - m1p;

						if (std::abs(m1) > l1) continue;
						if (std::abs(m1p) > l1p) continue;
						if (std::abs(m2) > l2) continue;
						if (std::abs(m2p) > l2p) continue;

						dlta1 = LaserHelper::kronecker(m1p - 1, m1);
						dlta2 = LaserHelper::kronecker(m1p + 1, m1);	

					
						double curCg = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, Mp);
						/*
						 * 
						 */
						if ((l1 - l1p) % 2 == 1)
						{
							/*
							 * Integral I2 in r1.
							 */
							norms = LegendreNorm(l1,m1) * LegendreNorm(l1p,m1p);

							temp = M_PI * m1 * norms * (dlta1 - dlta2); 
							temp *= VelGaugeXIntegrals::K1(l1p,m1p,l1,m1);


							I2_1 += temp * dlta_l2;

							/*
							 * Integral J1 in r1.
							 */
							F_lm = LaserHelperVelGaugeX::F(l1,m1);
							G_lm = LaserHelperVelGaugeX::G(l1,m1);
							H_lm = LaserHelperVelGaugeX::H(l1,m1);
							I_lm = LaserHelperVelGaugeX::I(l1,m1);

							K1_term1 = VelGaugeXIntegrals::K1(l1p,std::abs(m1p),l1-2,std::abs(m1));
							K1_term2 = VelGaugeXIntegrals::K1(l1p,std::abs(m1p),l1,std::abs(m1));
							K1_term3 = VelGaugeXIntegrals::K1(l1p,std::abs(m1p),l1+2,std::abs(m1));
				
							temp = LegendreNorm(l1-2,m1) * F_lm * K1_term1;
							temp += LegendreNorm(l1,m1) * (G_lm + H_lm) * K1_term2;
							temp += LegendreNorm(l1+2,m1) * I_lm * K1_term3;
							temp *= M_PI * (dlta1 + dlta2) * LegendreNorm(l1p,m1p);
						
							I3_1 += std::abs(m1) * temp * dlta_l2;

							/*
							 * Integral J2 in r1.
							 */
							J_lm = LaserHelperVelGaugeX::J(l1,m1);
							K_lm = LaserHelperVelGaugeX::K(l1,m1);
							E_lm = LaserHelperVelGaugeX::E(l1,m1);
						
							dlta_m = LaserHelperVelGaugeX::delta_m(m1);

							K2_term1 = VelGaugeXIntegrals::K2(l1-1,std::abs(m1+dlta_m),l1p,std::abs(m1p));
							K2_term2 = VelGaugeXIntegrals::K2(l1+1,std::abs(m1+dlta_m),l1p,std::abs(m1p));

							temp = LegendreNorm(l1-1,m1+dlta_m) * J_lm * K2_term1;
							temp += LegendreNorm(l1+1,m1+dlta_m) * K_lm * K2_term2;
							temp *= M_PI * LegendreNorm(l1p,m1p) * (dlta1 + dlta2);

							I3_1 -= dlta_m * E_lm * temp * dlta_l2;
						}


						if ((l2 - l2p) % 2 == 1)
						{
							/*
							 * Integral I2 in r2.
							 */
							norms = LegendreNorm(l2,m2) * LegendreNorm(l2p,m2p);
							
							temp = M_PI * m2 * norms * (dlta1 - dlta2);
							temp *= VelGaugeXIntegrals::K1(l2p,m2p,l2,m2); 

							I2_2 += temp * dlta_l1;

							/*
							 * Integral J1 in r2.
							 */
							F_lm = LaserHelperVelGaugeX::F(l2,m2);
							G_lm = LaserHelperVelGaugeX::G(l2,m2);
							H_lm = LaserHelperVelGaugeX::H(l2,m2);
							I_lm = LaserHelperVelGaugeX::I(l2,m2);
			
							K1_term1 = VelGaugeXIntegrals::K1(l2p,std::abs(m2p),l2-2,std::abs(m2));
							K1_term2 = VelGaugeXIntegrals::K1(l2p,std::abs(m2p),l2,std::abs(m2));
							K1_term3 = VelGaugeXIntegrals::K1(l2p,std::abs(m2p),l2+2,std::abs(m2));

							temp = LegendreNorm(l2-2,m2) * F_lm * K1_term1;
							temp += LegendreNorm(l2,m2) * (G_lm + H_lm) * K1_term2;
							temp += LegendreNorm(l2+2,m2) * I_lm * K1_term3;
							temp *= M_PI * (dlta1 + dlta2) * LegendreNorm(l2p,m2p);
	
							I3_2 += std::abs(m2) * temp * dlta_l1;


							/*
							 * Integral J2 in r2.
							 */
							J_lm = LaserHelperVelGaugeX::J(l2,m2);
							K_lm = LaserHelperVelGaugeX::K(l2,m2);
							E_lm = LaserHelperVelGaugeX::E(l2,m2);
						
							dlta_m = LaserHelperVelGaugeX::delta_m(m2);

							K2_term1 = VelGaugeXIntegrals::K2(l2-1,std::abs(m2+dlta_m),l2p,std::abs(m2p));
							K2_term2 = VelGaugeXIntegrals::K2(l2+1,std::abs(m2+dlta_m),l2p,std::abs(m2p));

							temp = LegendreNorm(l2-1,m2+dlta_m) * J_lm * K2_term1;
							temp += LegendreNorm(l2+1,m2+dlta_m) * K_lm * K2_term2;
							temp *= M_PI * LegendreNorm(l2p,m2p) * (dlta1 + dlta2);

							I3_2 -= dlta_m * E_lm * temp * dlta_l1;
						}


						if (l1p >= std::abs(l1-1) && l1p <= l1 + 1)
						{
							/*
							 * Integral I1 in r1.
							 */
							double Coefficient = VelGaugeXIntegrals::Coefficient(l1p,l1);
							temp = Coefficient * cg(l1,1,0,0,l1p,0);
							temp *= (cg(l1,1,m1,-1,l1p,m1p) - cg(l1,1,m1,1,l1p,m1p));
							I1_1 += temp * dlta_l2;
						}


						if (l2p >= std::abs(l2-1) && l2p <= l2 + 1)
						{
							/*
							 * Integral I2 in r2.
							 */
							double Coefficient = VelGaugeXIntegrals::Coefficient(l2p,l2); 
							temp = Coefficient * cg(l2,1,0,0,l2p,0);
							temp *= (cg(l2,1,m2,-1,l2p,m2p) - cg(l2,1,m2,1,l2p,m2p));
							I1_2 += temp * dlta_l1;
						}
						
						coupling1 += (-I1_1 - I2_1 + I3_1) * curCg;
						coupling2 += (-I1_2 - I2_2 + I3_2) * curCg;


					}
				}
			}

			//cout << M << " " << Mp << endl;
			double tol = 0.0000000000001; 
			if (std::abs(coupling1-0.0)< tol && std::abs(coupling2-0.0)< tol)
			//if (std::abs(I2_1-0.0)< tol && std::abs(I2_2-0.0)< tol) 
			{
				cout << "coupling1 " << l1 << "," << l2 << "," << L << "," << M << " ||| " << l1p << "," << l2p << "," << Lp << "," << Mp<< ": " << coupling1 << " " << coupling2 << endl;
			}
			else
			{
				cout << "xxxxxxxxx " << l1 << "," << l2 << "," << L << "," << M << " ||| " << l1p << "," << l2p << "," << Lp << "," << Mp << ": " << coupling1 << " " << coupling2 << endl;
			}

			//if (std::abs(coupling2-0.0)< 0.0000000000001 )
			//{
			//	cout << "coupling2 " << l1 << "," << l2 << "," << L << "," << M << " ||| " << l1p << "," << l2p << "," << Lp << "," << Mp << endl;
			//}



			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;
				double r1 = localr1(ri1);

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
					double r2 = localr2(ri2);
				
					data(index) = - IM * (coupling1/r1 + coupling2/r2);
				}
			}
		}
	}

	

	static double LegendreNorm(double l, double m)
	{
		/*
		 * Y_{l,m} = LegendreNorm(l,m) * P_l^m
		 */

		if (std::abs(m)<=l && l>=0)
		{
			double norm;
			norm = std::pow(-1.,.5 * (m + std::abs(m))) * sqrt((2. * l + 1.) / (4. * M_PI));
			norm *= std::sqrt(exp(gsl_sf_lnfact(l - std::abs(m)) - gsl_sf_lnfact(l + std::abs(m))));
			return norm;
		}
		else
		{
			return 0.;
		}
	}

	//static double Coefficient(double lp, double l)
	//{
	//	return std::sqrt(.5 * (2. * l + 1.) / (2. * lp + 1.) );
	//}


};





/* Second part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * Should be used with first order differentiation in r1
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR1_X
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocityDerivativeR1_X() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR1_X() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers 
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			if (std::abs(Mp - M) != 1) continue;

			double coupling1 = 0;
			double I1_1 = 0;
			double temp;

		 	double dlta_l2 = LaserHelper::kronecker(l2p,l2);

			for (int m1=-l1; m1<=l1; m1++)
			{
				for (int m1p = (m1 - 1); m1p <= (m1 + 1); m1p++)
				{
					if (m1p != m1)
					{
						int m2 = M - m1;
						int m2p = Mp - m1p;

						if (std::abs(m1) > l1) continue;
						if (std::abs(m1p) > l1p) continue;
						if (std::abs(m2) > l2) continue;
						if (std::abs(m2p) > l2p) continue;



						double curCg = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, Mp);
						
						if (l1p >= std::abs(l1-1) && l1p <= l1 + 1)
						{
							/*
							 * Integral I1 in r1.
							 */
							double Coefficient = VelGaugeXIntegrals::Coefficient(l1p,l1);
							temp = Coefficient * cg(l1,1,0,0,l1p,0);
							temp *= (cg(l1,1,m1,-1,l1p,m1p) - cg(l1,1,m1,1,l1p,m1p));
							I1_1 += temp * dlta_l2;
						}
						
						coupling1 += I1_1 * curCg;
			 		}
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
				
					data(index) = - IM * (coupling1);
				}
			}
		}
	}
};

/* Third part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * Should be used with first order differentiation in r2
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR2_X
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocityDerivativeR2_X() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR2_X() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers 
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			if (std::abs(Mp - M) != 1) continue;

			double coupling2 = 0;
			double I1_2 = 0;
			double temp;

			double dlta_l1 = LaserHelper::kronecker(l1p,l1);

			for (int m1=-l1; m1<=l1; m1++)
			{
				for (int m1p = (m1 - 1); m1p <= (m1 + 1); m1p++)
				{
					if (m1p != m1)
					{
						int m2 = M - m1;
						int m2p = Mp - m1p;

						if (std::abs(m1) > l1) continue;
						if (std::abs(m1p) > l1p) continue;
						if (std::abs(m2) > l2) continue;
						if (std::abs(m2p) > l2p) continue;

						double curCg = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, Mp);

						if (l2p >= std::abs(l2-1) && l2p <= l2 + 1)
						{
							/*
							 * Integral I2 in r2.
							 */
							double Coefficient = VelGaugeXIntegrals::Coefficient(l2p,l2);
							temp = Coefficient * cg(l2,1,0,0,l2p,0);
							temp *= (cg(l2,1,m2,-1,l2p,m2p) - cg(l2,1,m2,1,l2p,m2p));
							I1_2 += temp * dlta_l1;
						}
						
						coupling2 += I1_2 * curCg;
					}
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
				
					data(index) = - IM * (coupling2);
				}
			}
		}
	}
};
