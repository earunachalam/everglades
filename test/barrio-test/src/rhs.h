#ifndef RHS_H
#define RHS_H

#include <vector>

#include "Field.h"

struct constants
{
	// constants for reaction and diffusion terms
	double alpha,
		   beta,
		   gamma,
		   delta,
		   r1,
		   r2,
		   d;

	// absolute time used for t-dep source
	double time;
} c;

double rhs(std::vector<std::vector<std::shared_ptr<Field> > >& tseries, uint tprev_idx, uint tlag_idx, uint field_idx, uint x1idx, uint x2idx, constants p)
{

	// values for target cell and neighbors for field c = u or v
	// **at current time** - used for diffusion calculations
	double c 	= tseries.at(tprev_idx).at(field_idx)->at(x1idx,x2idx);				// concentration in cell of interest
	double clx1 = tseries.at(tprev_idx).at(field_idx)->x1less(x1idx, x2idx);		// in cell one index lower in x1 direction
	double cmx1 = tseries.at(tprev_idx).at(field_idx)->x1more(x1idx, x2idx);		// in cell one index higher in x1 direction
	double clx2 = tseries.at(tprev_idx).at(field_idx)->x2less(x1idx, x2idx);		// in cell one index lower in x2 direction
	double cmx2 = tseries.at(tprev_idx).at(field_idx)->x2more(x1idx, x2idx);		// in cell one index higher in x2 direction

	double caxp = tseries.at(tprev_idx).at(field_idx)->ax1;						// lattice constant for x direction
	double cayp = tseries.at(tprev_idx).at(field_idx)->ax2;						// for y direction

	// terms in gradient
	double d2c_dx2_lx1 = (clx1-c)/std::pow(caxp,2.0);					// rate in from cell one index lower in x1 direction
	double d2c_dx2_mx1 = (cmx1-c)/std::pow(caxp,2.0);					// one index higher
	double d2c_dx2_lx2 = (clx2-c)/std::pow(cayp,2.0);					// rate in from cell one index lower in x2 direction
	double d2c_dx2_mx2 = (cmx2-c)/std::pow(cayp,2.0);					// one index higher

	double grad2c = d2c_dx2_lx1 + d2c_dx2_mx1 + d2c_dx2_lx2 + d2c_dx2_mx2;		// total rate in from diffusion

	double c_kin		 = tseries.at(tlag_idx).at(field_idx)->at(x1idx,x2idx);	// concentration of species w/index field_idx in cell of interest
	double other_c_kin	 = tseries.at(tlag_idx).at(1-field_idx)->at(x1idx,x2idx);	// concentration of other species in cell of interest

	// eventual return value
	double rhs_operator_value;

	// account for differences in reaction and diffusion between the two species
	if (field_idx == 0 /* is field u */)
	{
		rhs_operator_value = p.delta*p.d*grad2c + p.alpha*c_kin*(1.00 - p.r1*std::pow(other_c_kin,2.00)) + other_c_kin*(1.00 - p.r2*c_kin);
	}
	else if (field_idx == 1 /* is field v */)
	{
		rhs_operator_value = p.delta*grad2c + p.beta*c_kin*(1.00 + p.alpha*p.r1*c_kin*other_c_kin/p.beta) + other_c_kin*(p.gamma + p.r2*c_kin);
	}
	else
	{
		std::cerr << __FILE__ << " L " << __LINE__ << " Error: nonexistent field index called." << std::endl;
		abort();
	}

	return rhs_operator_value;
}

#endif
