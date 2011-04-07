
//including interfacer
#include "coupledbase.h"

//Including legendre polynomials
#include<gsl/gsl_sf_legendre.h>

//Including laserhelper for kroenicker
#include "laserhelper.h"

/*
 * -1 / |r + R/2| - 1/|r - R/2|
 */

using namespace CoupledSpherical;

template<int Rank>
class DiatomicCoulomb : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

	DiatomicCoulomb() {}
	virtual ~DiatomicCoulomb() {}
		
	double R;
	double ThetaR;
	int pluss;
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialCoupledSphericalBase<Rank>::ApplyConfigSection
			(config);
		config.Get("inter_nuclear_r", R);
		config.Get("theta_inter_nucl", ThetaR);
		config.Get("pluss", pluss);
	}
	
	virtual void UpdatePotentialData(typename blitz::Array<cplx,Rank> data,
	   typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		
		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSpherical::CoupledSphericalHarmonicRepresentation
			CplHarmRepr;
		typename CmbRepr::Ptr Repr = boost::static_pointer_cast<CmbRepr>
			(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr=boost::static_pointer_cast<CplHarmRepr>
			(Repr->GetRepresentation(this->AngularRank));
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount= data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->
			GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->
			GetLocalGrid(this->RadialRank2);
	
		ClebschGordan cg;
		BasisPairList angBasisPairs = GetBasisPairList
			(this->AngularRank);
		if(data.extent(this->RadialRank1) != r1Count)
			throw std::runtime_error("Invalid r1 size");
		if(data.extent(this->RadialRank2) != r2Count)
			throw std::runtime_error("Invalid r2 size");
		if(data.extent(this->AngularRank) != angBasisPairs.extent(0))
			throw std::runtime_error("Invalid ang size");
		
		blitz::TinyVector<int, Rank> index;
		data = 0;
		double 	R_half =  R/2.0;
		//Looping over angular indices
		for(int angIndex = 0; angIndex < angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;
			
			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex= angBasisPairs(angIndex, 1);

			CoupledIndex left = angRepr->Range.GetCoupledIndex
				(leftIndex);
			CoupledIndex right =angRepr->Range.GetCoupledIndex
				(rightIndex);

			double cosTheta = std::cos(ThetaR);

			int Lp = left.L;
			int Mp = left.M;
			int l1p =left.l1;
			int l2p =left.l2;
			
			int L = right.L;
			int M = right.M;
			int l1 =right.l1;
			int l2 =right.l2;

			int minL3 = std::abs(l1 - l1p);
			int maxL3 = l1 + l1p;
			for(int l3 = minL3; l3<=maxL3; l3++)
			{
				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l1,l1p);
				l3Coeff *= cg(l1, l3, 0, 0, l1p, 0);
				l3Coeff *= LaserHelper::kronecker(l2,l2p);
			
				double l3Sum = 0;
				
				for(int m1p = -l1p; m1p <= l1p; m1p++)
				{
					int m2p = Mp - m1p;
					for(int m1 = -l1; m1 <= l1; m1++)
					{
						int m2 = M - m1;
						int m3 = m1p - m1;
						
						if(l3 % 2 ==1) continue;
						if(std::abs(m1) > l1) continue;
						if(std::abs(m1p)>l1p) continue;
						if(std::abs(m2) > l2) continue;
						if(std::abs(m2p)>l2p) continue;
						if(std::abs(m3) > l3) continue;
			
						double cur = 1; 
	
						cur *= gsl_sf_legendre_sphPlm
						   (l3,std::abs(m3),cosTheta);
						cur *= 2.0;
						cur *=CondonShortleyPhase(-m3);
						cur *= MultipoleCoeff(l3);
						
						cur *= cg(l1p,l2p,
							m1p,m2p,Lp,Mp);
						cur *= cg(l1,l2,m1,m2,L,M);
						cur *= cg(l1,l3,m1,m3,l1p,m1p);
						
						cur *= LaserHelper::kronecker
							(m2, m2p);
						
						l3Sum += cur;
					}
				}
				
				l3Sum *= l3Coeff;
			
				for(int ri1=0; ri1 < r1Count; ri1++)
				{
					
					double r1 = localr1(ri1);
					index(this->RadialRank1) = ri1;
					for(int ri2=0; ri2 < r2Count; ri2++)
					{
						//double r2 = localr2(ri2);
						index(this->RadialRank2) = ri2;
						
						double rmin=std::min(r1,R_half);
						double rmax=std::max(r1,R_half);
						double rfrac = rmin / rmax;
						
						data(index) += -1. * l3Sum * 
						  std::pow(rfrac, l3) / rmax;
					}
				}
			}
		
			// r_1 <=> r_2
	
			minL3 = std::abs(l2 - l2p);
			maxL3 = l2 + l2p;
			
			for(int l3 = minL3; l3 <= maxL3; l3++)
			{
				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l2,l2p);
				l3Coeff *= cg(l2,l3,0,0,l2p,0);
				l3Coeff *= LaserHelper::kronecker(l1,l1p);
			
				double l3Sum = 0;
				
				for(int m1p = -l1p; m1p <= l1p; m1p++)
				{
					int m2p = Mp - m1p;
					
					for(int m1= -l1; m1 <= l1; m1++)
					{
						int m2 = M - m1;
						int m3 = m2p - m2;

						//makes sure that only even ls are non-zero
						if((l3) % 2 == 1) continue;
						if(std::abs(m1)>l1) continue;
						if(std::abs(m1p)>l1p) continue;
						if(std::abs(m2) > l2) continue;
						if(std::abs(m2p)>l2p) continue;
						if(std::abs(m3) > l3) continue;
					
						
						double cur = 1; 
						cur *= gsl_sf_legendre_sphPlm
						  (l3,std::abs(m3),cosTheta);
						cur *= CondonShortleyPhase(-m3);
						cur *= MultipoleCoeff(l3);
						cur *= 2.0;
						cur *=cg(l1p,l2p,m1p,m2p,Lp,Mp);
						cur *= cg(l1,l2,m1,m2,L,M);
						cur *= cg(l2,l3,m2,m3,l2p,m2p);
						
						cur *= LaserHelper::kronecker
							(m1,m1p);
					
						l3Sum += cur;
					}
				}
				l3Sum *= l3Coeff;
			
				for(int ri1 = 0; ri1 < r1Count; ri1++)
				{
					//double r1 = localr1(ri1);
					index(this->RadialRank1) = ri1;
					
					for(int ri2 = 0; ri2 < r2Count; ri2++)
					{
						double r2 = localr2(ri2);
						index(this->RadialRank2) = ri2;
						
						double rmin=std::min(r2,R_half);
						double rmax=std::max(r2,R_half);
						double rfrac = rmin / rmax;
						data(index) += -1. *l3Sum *
						 std::pow(rfrac,l3) /rmax;
					}
				}
			}
		}
	}

	static double Coefficient(int a, int b)
	{
		return std::sqrt((2. *a + 1.) / (2. *b +1.));
	}

	static double MultipoleCoeff(int c)
	{
		return std::sqrt((4. * M_PI) / (2. * c + 1.));
	}
	
	static double CondonShortleyPhase(int m)
	{
		if(m < 0) return 1.0;
		return std::pow(-1.0, m);
	}

};
