#include <armadillo>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "Field.h"
#include "rk4_step.h"

double dt = 0.005;
double ttot = 30.0;


double diffusion(std::vector<std::vector<Field> >& tseries, uint t, uint field_idx, uint x1idx, uint x2idx)
{
	// constants for reaction and diffusion terms
	double alpha = 0.899,
		   beta = -0.91,
		   gamma = -0.899,
		   delta = 2.0,
		   r1 = 0.02,
		   r2 = 0.2,
		   d = 0.536;

	// values for target cell and neighbors for field c = u or v
	// **at current time** - used for diffusion calculations
	double c 	= tseries.at(t).at(field_idx).at(x1idx,x2idx);			// concentration in cell of interest
	double clx1 = tseries.at(t).at(field_idx).x1less(x1idx, x2idx);		// in cell one index lower in x1 direction
	double cmx1 = tseries.at(t).at(field_idx).x1more(x1idx, x2idx);		// in cell one index higher in x1 direction
	double clx2 = tseries.at(t).at(field_idx).x2less(x1idx, x2idx);		// in cell one index lower in x2 direction
	double cmx2 = tseries.at(t).at(field_idx).x2more(x1idx, x2idx);		// in cell one index higher in x2 direction
	
	double caxp = tseries.at(t).at(field_idx).m_ax1;					// lattice constant for x direction
	double cayp = tseries.at(t).at(field_idx).m_ax2;					// for y direction

	// terms in gradient
	double d2c_dx2_lx1 = (c-clx1)/std::pow(caxp,2.0);
	double d2c_dx2_mx1 = (cmx1-c)/std::pow(caxp,2.0);
	double d2c_dx2_lx2 = (c-clx2)/std::pow(cayp,2.0);
	double d2c_dx2_mx2 = (c-cmx2)/std::pow(cayp,2.0);

	double grad2c = d2c_dx2_lx1 + d2c_dx2_mx1 + d2c_dx2_lx2 + d2c_dx2_mx2;
	
	double runtotal;

	if (field_idx == 0 /* is field u */)
	{

		runtotal = delta*d*grad2c + alpha*c*(;
	}
}


int main(/*int argc, char** argv*/)
{
	//double (*cro_inst)(std::vector<std::vector<double> >& fieldvalues, uint field_idx, uint cell_xidx, uint cell_yidx) = &cro_test;
	
	// dimensions of system (number of cells)
	uint Nx1 = 10, Nx2 = 10;
	
	Field u(Nx1, Nx2, 1.0, 1.0);
	Field v(Nx1, Nx2, 1.0, 1.0);
	
	std::vector<Field> fields = {u, v};
	
	std::vector<std::vector<Field> > tseries;
	tseries.push_back(fields);

	for (uint timestep = 1; timestep < static_cast<uint>(ttot/dt); ++timestep)
	{
		rk4_step(tseries, dt, diffusion);
	}

	//FILE* outfile = fopen("dat/y_of_t.dat","w");
	//for (uint ct_i_timestep = 0; ct_i_timestep < y1.size(); ++ct_i_timestep)
	//{
		//fprintf(outfile, "%lf %lf %lf\n", static_cast<double>(ct_i_timestep)*dt, y1.at(ct_i_timestep), y2.at(ct_i_timestep));
	//}
	
	return 0;
}
