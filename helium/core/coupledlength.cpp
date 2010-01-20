#include "coupledbase.h"

/*
 * Potential evaluator for linearly polarized length gauge electric field (z-direction)
 */
using namespace CoupledSpherical;

template<int Rank>
class CustomPotential_LaserLength : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotential_LaserLength() {}
	virtual ~CustomPotential_LaserLength() {}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{

		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSphericalHarmonicRepresentation CplHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr = boost::static_pointer_cast< CplHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(this->RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (left.M != right.M) continue;
	
			// "Left" quantum numbers
			int L = left.L;
			int M = left.M;
			int l1 = left.l1;
			int l2 = left.l2;
			
			// "Right" quantum numbers (Mp = M)
			int Lp = right.L;
			int l1p = right.l1;
			int l2p = right.l2;

			double I1 = 0;
			double I2 = 0;
			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			for (int m1=-lStop; m1<=lStop; m1++)
			{
				int m2 = M - m1;
				int m1p = m1;
				int m2p = m2;

				if (std::abs(m1) > l1) continue;
				if (std::abs(m1p) > l1p) continue;
				if (std::abs(m2) > l2) continue;
				if (std::abs(m2p) > l2p) continue;
				
				// Based on Bransden/Joachain 2nd ed. p. 205.
				double cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, M); 
				cur *= cg(l1, 1, 0, 0, l1p, 0) * cg(l1, 1, m1, 0, l1p, m1p);
				cur *= Coefficient(l1, l1p);
				I1 += cur;
				
				cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= cg(l2, 1, 0, 0, l2p, 0) * cg(l2, 1, m2, 0, l2p, m2p);
				cur *= Coefficient(l2, l2p);
				I2 += cur;

			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				double r1 = localr1(ri1);
				index(this->RadialRank1) = ri1;

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					double r2 = localr2(ri2);
					index(this->RadialRank2) = ri2;

					data(index) += (I1 * r1 + I2 * r2);
				}
			}
		}
	}

	static double Coefficient(int l, int lp)
	{
		//sqrt(4pi/3) included
		return std::sqrt((2 * l + 1.0) / (2 * lp + 1.0));
	}
};


/*
 * Potential evaluator for linearly polarized length gauge electric field (x-direction)
 */
template<int Rank>
class CustomPotential_LaserLength_X : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotential_LaserLength_X() {}
	virtual ~CustomPotential_LaserLength_X() {}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{

		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSphericalHarmonicRepresentation CplHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr = boost::static_pointer_cast< CplHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(this->RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (std::abs(left.M - right.M) != 1) continue;

			// "Left" quantum numbers
			int L = left.L;
			int M = left.M;
			int l1 = left.l1;
			int l2 = left.l2;
			
			// "Right" quantum numbers 
			int Lp = right.L;
			int l1p = right.l1;
			int l2p = right.l2;
			int Mp = right.M;

			double I1 = 0;
			double I2 = 0;

			for (int m1 = -l1; m1 <= l1; m1++)
			{
				int m2 = M - m1;
				
				for (int m1p = (m1 - 1); m1p <= (m1 + 1); m1p++)
				{	
					int m2p = Mp - m1p;

					if (std::abs(m2p - m2) > 1) continue;

					//For given Y_lm, |m|<=l
					if (std::abs(m1) > l1) continue;
					if (std::abs(m1p) > l1p) continue;
					if (std::abs(m2) > l2) continue;
					if (std::abs(m2p) > l2p) continue;

					double cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, Mp);
					cur *= Coefficient(l1, l1p);
					cur *= cg(l1, 1, 0, 0, l1p, 0); 
					cur *= (cg(l1, 1, m1, -1, l1p, m1p) - cg(l1, 1, m1, 1, l1p, m1p));	
					I1 += cur;

					cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, Mp);
					cur *= Coefficient(l2, l2p);
					cur *= cg(l2, 1, 0, 0, l2p, 0);
				        cur *= (cg(l2, 1, m2, -1, l2p, m2p) - cg(l2, 1, m2, 1, l2p, m2p));
					I2 += cur;
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				double r1 = localr1(ri1);
				index(this->RadialRank1) = ri1;

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					double r2 = localr2(ri2);
					index(this->RadialRank2) = ri2;

					data(index) += (I1 * r1 + I2 * r2);
				}
			}
		}
	}

	static double Coefficient(int l, int lp)
	{
		return std::sqrt((2 * l + 1.0 ) / (2*(2 * lp + 1.0)));
	}
};

