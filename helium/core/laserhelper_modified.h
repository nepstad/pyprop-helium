#ifndef LASERHELPER_H
#define LASERHELPER_H

#include <gsl/gsl_sf_gamma.h>

class LaserHelper
{
public:
	static double C(double l, double m)
	{
		return l * std::sqrt( ((l+1+m)*(l+1-m)) / ((2*l+1)*(2*l+3)) );
	}

	static double D(double l, double m)
	{
		return -(l+1) * std::sqrt( ((l+m)*(l-m)) / ((2*l+1)*(2*l-1)) );
	}

	static double E(double l, double m)
	{
		return std::sqrt( ((l+1+m)*(l+1-m)) / ((2*l+1)*(2*l+3)) );
	}

	static double F(double l, double m)
	{
		return std::sqrt( ((l+m)*(l-m)) / ((2*l+1)*(2*l-1)) );
	}

	static int kronecker(int a, int b)
	{
		return a == b ? 1 : 0;
	}
};


class LaserHelperVelGaugeX
/*
* The names are according to the note 
* "Velocity gauge for a field along the x axis"
*/
{
public:

	static double E(double l, double m)
	{
		double j = (l - std::abs(m)) * (l + std::abs(m) + 1.); 
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( (l - std::abs(m)) * (l + std::abs(m) + 1.) );
		}
	}

	static double F(double l, double m)
	{
		double j = ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.)) );
		}
	}

	static double G(double l, double m)
	{
		double j = ((l+m)*(l-m)) / ((2.*l+1.)*(2.*l-1.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m)*(l-m)) / ((2.*l+1.)*(2.*l-1.)) );
		}
	}

	static double H(double l, double m)
	{
		double j = ((l+m+1.)*(l-m+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+1.)*(l-m+1.)) / ((2.*l+1.)*(2.*l+3.)) );
		}
	}

	static double I(double l, double m)
	{
		double j = ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.)) );
		}
	}

	static double delta_m(double m)
	{
		return m >= 0. ? 1. : -1.;
	}


	static double J(double l, double m)
	{
		double j = ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.)) );
		}
	}



	static double K(double l, double m)
	{
		double j = ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < 0)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.)) );
		}
	}

};

class VelGaugeXIntegrals
{
public:

	static double K1(int l, int m, int p, int q)
	{
		int c1 = l-p;
		int c2 = m-q;
		int nu = std::min(l,p);
		int mu = std::min(m,q);

		if ((l+p) % 2==1)
		{
			if (((c1 < 0) && (c2 < 0)) || ((c1 > 0) && (c2 > 0)))
			{
				return -2.0 * exp(gsl_sf_lnfact(nu+mu) - gsl_sf_lnfact(nu-mu));
			}
			else
			{
				return 0.0;
			}
		}
		else
		{
			return 0.;
		}
	}


	static double K2(int l, int m, int p, int q)
	{
		if ((l+p-m-q) % 2 == 0)
		{
			double outerSum = 0.0;
			int outerMax = std::floor((l-m)/2.);
			int innerMax = std::floor((p-q)/2.);

			for (int i=0; i<=outerMax; i++)
			{
				double innerSum = .0;
				for (int j=0; j<=innerMax; j++)
				{
					innerSum += Cconstant(p,q,j);
					innerSum *= exp(gsl_sf_lngamma(.5 * (l+p-m-q -2.*(i+j)+1.)) + gsl_sf_lngamma(.5 * (m+q+2.*(i+j+1.))) - gsl_sf_lngamma(.5*(l+p+3.)));
				}
				outerSum += innerSum * Cconstant(l,m,i); 
			}
			return outerSum;
		}
		else
		{
		 	return 0.0;	
		}
	}

	static double Cconstant(double alpha, double beta, double gamma)
	{
		double C;
		C = std::pow(-1.,gamma) * std::pow(2.,-(beta + 2 * gamma)); 
		C *= exp(gsl_sf_lnfact(alpha + beta) - gsl_sf_lnfact(beta + gamma)-gsl_sf_lnfact(gamma) - gsl_sf_lnfact(alpha - beta - 2. * gamma));
		return C;
	}

	static double Coefficient(double lp, double l)
	{
		return std::sqrt(.5 * (2. * l + 1.) / (2. * lp + 1.) );
	}



};

#endif

