#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "rk4_step.h"

double dt = 0.005;
double ttot = 30.0;

double cro_test(std::vector<std::vector<double> >& t0, unsigned int field_idx, unsigned int cell_idx)
{
	double Du = 1,
		   Dv = 40,
		   c1 = 0.1,
		   c2 = 0.9,
		   cm1 = 1,
		   c3 = 1;

	double u = t0.at(/* field idx = */ 0).at(cell_idx);
	double v = t0.at(/* field idx = */ 1).at(cell_idx);

	// use if's to use different expressions for different fields!
	if (field_idx == 0)
	{
		return c1 - cm1*u + c3*u*u*v;
	}
	else if (field_idx == 1)
	{
		return c2 - c3*u*u*v;
	}
	else
	{
		printf("Error\n");
		return 0.00;
	}
}

int main(int argc, char** argv)
{
	double (*cro_inst)(std::vector<std::vector<double> >& t0, unsigned int field_idx, unsigned int cell_idx) = &cro_test;
	
	std::vector<std::vector<double> > t0;
	t0.push_back({2});
	t0.push_back({0});
	
	std::vector<double> y1, y2;

	for (unsigned int timestep = 1; timestep < static_cast<unsigned int>(ttot/dt); ++timestep)
	{
		y1.push_back(t0.at(0).front());
		y2.push_back(t0.at(1).front());
		rk4_step(t0, cro_inst, dt);
	}

	FILE* outfile = fopen("dat/y_of_t.dat","w");
	for (unsigned int ct_i_timestep = 0; ct_i_timestep < y1.size(); ++ct_i_timestep)
	{
		fprintf(outfile, "%lf %lf %lf\n", static_cast<double>(ct_i_timestep)*dt, y1.at(ct_i_timestep), y2.at(ct_i_timestep));
	}
	
	return 0;
}
