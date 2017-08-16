#include <armadillo>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>
#include <omp.h>
#include <vector>

#include "Field.h"

double rhs(std::vector<std::vector<std::shared_ptr<Field> > >& tseries, uint t, uint field_idx, uint x1idx, uint x2idx)
{
	// constants for reaction and diffusion terms
	double alpha = 0.899,
		   beta = -0.91,
		   gamma = -0.899,
		   delta = 2.0,
		   r1 = 0.02,
		   r2 = 0.2,
		   d = 0.516;

	// values for target cell and neighbors for field c = u or v
	// **at current time** - used for diffusion calculations
	double c 	= tseries.at(t).at(field_idx)->at(x1idx,x2idx);			// concentration in cell of interest
	double clx1 = tseries.at(t).at(field_idx)->x1less(x1idx, x2idx);		// in cell one index lower in x1 direction
	double cmx1 = tseries.at(t).at(field_idx)->x1more(x1idx, x2idx);		// in cell one index higher in x1 direction
	double clx2 = tseries.at(t).at(field_idx)->x2less(x1idx, x2idx);		// in cell one index lower in x2 direction
	double cmx2 = tseries.at(t).at(field_idx)->x2more(x1idx, x2idx);		// in cell one index higher in x2 direction
	
	double caxp = tseries.at(t).at(field_idx)->ax1;						// lattice constant for x direction
	double cayp = tseries.at(t).at(field_idx)->ax2;						// for y direction

	// terms in gradient
	double d2c_dx2_lx1 = (clx1-c)/std::pow(caxp,2.0);					// rate in from cell one index lower in x1 direction
	double d2c_dx2_mx1 = (cmx1-c)/std::pow(caxp,2.0);					// one index higher
	double d2c_dx2_lx2 = (clx2-c)/std::pow(cayp,2.0);					// rate in from cell one index lower in x2 direction
	double d2c_dx2_mx2 = (cmx2-c)/std::pow(cayp,2.0);					// one index higher

	double grad2c = d2c_dx2_lx1 + d2c_dx2_mx1 + d2c_dx2_lx2 + d2c_dx2_mx2;		// total rate in from diffusion
	
	double other_c = tseries.at(t).at(1-field_idx)->at(x1idx,x2idx);			// concentration of other species in cell of interest
	
	// eventual return value
	double rhs_operator_value;

	// account for differences in reaction and diffusion between the two species
	if (field_idx == 0 /* is field u */)
	{
		//rhs_operator_value = delta*grad2c;
		rhs_operator_value = delta*d*grad2c + alpha*c*(1.00 - r1*std::pow(other_c,2.00)) + other_c*(1.00 - r2*c);
	}
	else if (field_idx == 1 /* is field v */)
	{
		//rhs_operator_value = delta*grad2c;
		rhs_operator_value = delta*grad2c + beta*c*(1.00 + alpha*r1*c*other_c/beta) + other_c*(gamma + r2*c);
	}
	else
	{
		std::cerr << __FILE__ << " L " << __LINE__ << " Error: nonexistent field index called." << std::endl;
		abort();
	}

	return rhs_operator_value;
}


int main()
{
	double dt = 0.001;
	double ttot = 2001;
	uint Ntstep = static_cast<uint>(ttot/dt);
	uint Ntstep_progress = static_cast<uint>(Ntstep/100);
	uint tstep_start_write = 0;
	uint tstep_end_write = Ntstep;
	uint tstep_intvl_write = 200;

	// dimensions of system (number of cells)
	uint Nx1 = 100, Nx2 = 100;
	
	// concentration fields
	std::vector<std::shared_ptr<Field> > field_ptrs;
	field_ptrs.push_back(std::make_shared<Field>(Nx1, Nx2, 1.0, 1.0));

	field_ptrs.push_back(std::make_shared<Field>(Nx1, Nx2, 1.0, 1.0));
	
	// store entire trajectory
	std::vector<std::vector<std::shared_ptr<Field> > > tseries;
	tseries.push_back(field_ptrs);

	// write initial configuration
	for (uint ifield = 0; ifield < 2; ++ ifield)
	{
		std::string filename = "dat/pnl_c/t0_f" + std::to_string(ifield) +".dat";
		std::ofstream datafile(filename);
		datafile << tseries.at(0).at(ifield)->value << std::endl;
		datafile.close();
	}

	// start timestepping
	for (uint tstep = 1; tstep < Ntstep; ++tstep)
	{
		tseries.push_back(tseries.back());
		
		for (uint ifield = 0; ifield < 2; ++ifield)
		{
			// zero-flux boundary condition
			tseries.back().at(ifield)->bound_x1less = tseries.back().at(ifield)->value.row(0);
			tseries.back().at(ifield)->bound_x2less = tseries.back().at(ifield)->value.col(0);
			tseries.back().at(ifield)->bound_x1more = tseries.back().at(ifield)->value.row(Nx1-1);
			tseries.back().at(ifield)->bound_x2more = tseries.back().at(ifield)->value.col(Nx2-1);

			// advance in time, computing rows in parallel
			#pragma omp parallel for
			for (uint x1idx = 0; x1idx < Nx1; ++x1idx)
			{
				for (uint x2idx = 0; x2idx < Nx2; ++x2idx)
				{
					tseries.at(tstep).at(ifield)->value(x1idx,x2idx) += dt*rhs(tseries, tstep-1, ifield, x1idx, x2idx);
				}
			}
			
			// write data
			if (((tstep-tstep_start_write)%tstep_intvl_write == 0) && (tstep < tstep_end_write))
			{
				std::string filename = "dat/pnl_c/t" + std::to_string(static_cast<int>(tstep*dt)) + "_f" + std::to_string(ifield) +".dat";
				std::ofstream datafile(filename);
				datafile << tseries.at(tstep).at(ifield)->value << std::endl;
				datafile.close();
			}
		}

		// display progress
		if (tstep%Ntstep_progress == 0)
		{
			printf("\r %2.0f%% complete", 100.*tstep/Ntstep); fflush(stdout);
		}
	}
	printf("\r 100%% complete"); fflush(stdout);

	return 0;
}
