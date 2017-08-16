#include <cmath>
#include <csignal>
#include <iostream>
#include <vector>

#define GC(A) {int status; char * demangled = abi::__cxa_demangle(typeid(A).name(),0,0,&status); std::cout << __LINE__ << ": " #A << "\t" << demangled <<"\    n";}                                                                                                                                                       
#define P2S(a) std::cout << __LINE__ << ": " << #a << ": " << (a) << std::endl

void rk4_step(std::vector<std::vector<double> >& t0, double (*change_rate_operator)(std::vector<std::vector<double> >& t0, unsigned int field_idx, unsigned int cell_idx), double dt)
{
	// initialize arrays for the four approximations
	std::vector<std::vector<double> >
		k1 = t0,
		k1mod = t0,
		k2 = t0,
		k2mod = t0,
		k3 = t0,
		k3mod = t0,
		k4 = t0;

	// calculate k1
	for (unsigned int ct_i_field = 0; ct_i_field < t0.size(); ++ct_i_field)
	{
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k1.at(ct_i_field).at(ct_j_cell) = change_rate_operator(t0, ct_i_field, ct_j_cell);
		}
	}

	// calculate k2
	for (unsigned int ct_i_field = 0; ct_i_field < t0.size(); ++ct_i_field)
	{
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k1mod.at(ct_i_field).at(ct_j_cell) = t0.at(ct_i_field).at(ct_j_cell) + (k1.at(ct_i_field).at(ct_j_cell)*(dt/2.00));
		}

		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k2.at(ct_i_field).at(ct_j_cell) = change_rate_operator(k1mod, ct_i_field, ct_j_cell);
		}
	}
	
	// calculate k3
	for (unsigned int ct_i_field = 0; ct_i_field < t0.size(); ++ct_i_field)
	{
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k2mod.at(ct_i_field).at(ct_j_cell) = t0.at(ct_i_field).at(ct_j_cell) + (k2.at(ct_i_field).at(ct_j_cell)*(dt/2.00));
		}
		
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k3.at(ct_i_field).at(ct_j_cell) = change_rate_operator(k2mod, ct_i_field, ct_j_cell);
		}
	}
	
	// calculate k4
	for (unsigned int ct_i_field = 0; ct_i_field < t0.size(); ++ct_i_field)
	{
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k3mod.at(ct_i_field).at(ct_j_cell) = t0.at(ct_i_field).at(ct_j_cell) + (k3.at(ct_i_field).at(ct_j_cell)*dt);
		}
		
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			k4.at(ct_i_field).at(ct_j_cell) = change_rate_operator(k3mod, ct_i_field, ct_j_cell);
		}
	}

	// calculate weighted sum
	for (unsigned int ct_i_field = 0; ct_i_field < t0.size(); ++ct_i_field)
	{
		for (unsigned int ct_j_cell = 0; ct_j_cell < t0.at(ct_i_field).size(); ++ct_j_cell)
		{
			// calculate weighted change in each value
			double m = (1.00/6.00)*(k1.at(ct_i_field).at(ct_j_cell) + 2.00*k2.at(ct_i_field).at(ct_j_cell) + 2.00*k3.at(ct_i_field).at(ct_j_cell) + k4.at(ct_i_field).at(ct_j_cell));

			// update value by adding change
			t0.at(ct_i_field).at(ct_j_cell) += m*dt;
		}
	}

}
