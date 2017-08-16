#include <armadillo>
#include <cmath>
#include <csignal>
#include <iostream>
#include <vector>

#include "Field.h"

#define GC(A) {int status; char * demangled = abi::__cxa_demangle(typeid(A).name(),0,0,&status); std::cout << __LINE__ << ": " #A << "\t" << demangled <<"\n";}                                                                                           
#define P2S(a) std::cout << __LINE__ << ": " << #a << ": " << (a) << std::endl

void rk4_step(std::vector<std::vector<Field> >& tseries, uint t, double (*change_rate_operator)(std::vector<std::vector<Field> >& tseries, uint field_idx, uint cell_x1idx, uint cell_x2idx), double dt)
{
	// initialize arrays for the four approximations
	std::vector<Field>
		k1 = tseries.at(t),
		k1mod = tseries.at(t),
		k2 = tseries.at(t),
		k2mod = tseries.at(t),
		k3 = tseries.at(t),
		k3mod = tseries.at(t),
		k4 = tseries.at(t);

	// calculate k1
	for (uint fieldidx = 0; fieldidx < tseries.at(t).size(); ++fieldidx)
	{
		for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		{
			for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			{
				k1.at(fieldidx).at(x1idx,x2idx) = change_rate_operator(tseries, fieldidx, x1idx, x2idx);
			}
		}
	}

	// calculate k2
	for (uint fieldidx = 0; fieldidx < tseries.at(t).size(); ++fieldidx)
	{
		k1mod.at(fieldidx).m_value = tseries.at(t).at(fieldidx).m_value + (k1.at(fieldidx).m_value*(dt/2.00));
		//for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		//{
			//for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			//{
				//k1mod.at(fieldidx).at(ct_j_cell) = tseries.at(t).at(fieldidx).at(ct_j_cell) + (k1.at(fieldidx).at(ct_j_cell)*(dt/2.00));
			//}
		//}

		for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		{
			for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			{
				k2.at(fieldidx).at(x1idx,x2idx) = change_rate_operator(k1mod, fieldidx, x1idx, x2idx);
			}
		}
	}
	
	// calculate k3
	for (uint fieldidx = 0; fieldidx < tseries.at(t).size(); ++fieldidx)
	{
		k2mod.at(fieldidx).m_value = tseries.at(t).at(fieldidx).m_value + (k2.at(fieldidx).m_value*(dt/2.00));
		//for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		//{
			//for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			//{
				//k2mod.at(fieldidx).at(ct_j_cell) = tseries.at(t).at(fieldidx).at(ct_j_cell) + (k2.at(fieldidx).at(ct_j_cell)*(dt/2.00));
			//}
		//}
		
		for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		{
			for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			{
				k3.at(fieldidx).at(x1idx,x2idx) = change_rate_operator(k2mod, fieldidx, x1idx, x2idx);
			}
		}
	}
	
	// calculate k4
	for (uint fieldidx = 0; fieldidx < tseries.at(t).size(); ++fieldidx)
	{
		//for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		//{
			//for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			//{
				//k3mod.at(fieldidx).at(ct_j_cell) = tseries.at(t).at(fieldidx).at(ct_j_cell) + (k3.at(fieldidx).at(ct_j_cell)*dt);
			//}
		//}
		
		for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		{
			for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			{
				k4.at(fieldidx).at(x1idx,x2idx) = change_rate_operator(k3mod, fieldidx, x1idx, x2idx);
			}
		}
	}

	// calculate weighted sum
	for (uint fieldidx = 0; fieldidx < tseries.at(t).size(); ++fieldidx)
	{
		for (uint x1idx = 0; x1idx < tseries.at(t).at(fieldidx).m_Nx1; ++x1idx)
		{
			tseries.at(t).at(fieldidx).add_value(dt*(1.00/6.00)*(k1.at(fieldidx).m_value + 2.00*k2.at(fieldidx).m_value + 2.00*k3.at(fieldidx).m_value + k4.at(fieldidx).m_value));

			//for (uint x2idx = 0; x2idx < tseries.at(t).at(fieldidx).m_Nx2; ++x2idx)
			//{
				//// calculate weighted change in each value
				//double m = (1.00/6.00)*(k1.at(fieldidx).at(ct_j_cell) + 2.00*k2.at(fieldidx).at(ct_j_cell) + 2.00*k3.at(fieldidx).at(ct_j_cell) + k4.at(fieldidx).at(ct_j_cell));

				//// update value by adding change
				//tseries.at(t).at(fieldidx).at(ct_j_cell) += m*dt;
			//}
		}
	}

}
