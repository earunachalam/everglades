#include <cassert>
#include <iostream>
#include <vector>

int main()
{
	double dt = 0.01;										// timestep size in system (natural) units
	double ttot = 0.33;										// total amount of time to integrate
	uint Ntstep = static_cast<uint>(ttot/dt);				// number of timesteps necessary to reach this time

	double tau = 0.03;										// size of time delay in natural units
	uint tau_tstep = static_cast<uint>(tau/dt);				// corresponding number of timesteps

	std::vector<int> v(tau_tstep+1);
	// create initial configuration and history
	for (int tstep = -tau_tstep; tstep < 1; ++tstep)
	{
		v.at(tstep+tau_tstep) = tstep;
	}
	for (auto& elem: v) std::cout << elem << "\t"; std::cout << std::endl;

	// start timestepping
	for (uint tstep = 1; tstep < Ntstep; ++tstep)
	{
		for (auto& elem: v) std::cout << elem << "\t";
		std::cout << "ts = " << tstep << " ";
		std::cout << "prv = " << ((tstep-1+tau_tstep)%(tau_tstep+1)) << " ";
		std::cout << "lag = " << ((tstep-1)%(tau_tstep+1)) << std::endl;
		v.at((tstep-1)%(tau_tstep+1)) = tstep;
	}

	return 0;
}
