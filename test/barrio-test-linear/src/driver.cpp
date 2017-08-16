#include <armadillo>
#include <cassert>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <omp.h>
#include <vector>

#include "Field.h"
#include "rhs.h"

#include <cxxabi.h>
#define GC(A) {int status; char * demangled = abi::__cxa_demangle(typeid(A).name(),0,0,&status); std::cout << __LINE__ << ": " #A << "\t" << demangled <<"\n";}                                                                                                                                                       

// loads the configuration file given by 'filename' and returns a map of key-value pairs
std::map<std::string,std::string> load_config(FILE* fp)
{
	std::map<std::string,std::string> opts;

	char key[100];
	char value[100];
	int line = 1;
    while(!feof(fp))
	{
		if (fscanf(fp, "%s %s\n", key, value) == 2)
		{
			std::string strkey(key), strval(value);
			opts[strkey] = std::string(strval);
			++line;
		}
		else
		{
			printf("Error reading configuration file at line %d\n.", line);
			abort();
		}
	}
	return opts;
}

int main(int argc, char* argv[])
{
	if (argc<2) { printf("Error: need to provide name of configuration file.\n"); abort(); }

	FILE* fp = fopen(argv[1], "r");
	if (fp == nullptr) { printf("Error: File could not be opened."); abort(); }

	std::map<std::string,std::string> params = load_config(fp);

	double dt = std::stod(params.find("dt")->second);		// timestep size in system (natural) units
	double ttot = std::stod(params.find("ttot")->second);	// total amount of time to integrate
	uint Ntstep = static_cast<uint>(ttot/dt);				// number of timesteps necessary to reach this time
	uint Ntstep_progress = static_cast<uint>(Ntstep/100);	// number of timesteps between successive progress updates (/100 = every percent, etc.)

	double tau = std::stod(params.find("tau")->second);		// size of time delay in natural units
	uint tau_tstep = static_cast<uint>(tau/dt);				// corresponding number of timesteps

	uint tstep_start_write = std::stoi(params.find("tstep_start_write")->second);
	uint tstep_end_write = Ntstep - std::stoi(params.find("tstep_endb4_write")->second);
	uint tstep_intvl_write = std::stoi(params.find("tstep_intvl_write")->second);
	std::string datadir = params.find("datadir")->second;

	constants c;
	c.alpha = std::stod(params.find("alpha")->second);
	c.beta = std::stod(params.find("beta")->second);
	c.gamma = std::stod(params.find("gamma")->second);
	c.delta = std::stod(params.find("delta")->second);
	c.r1 = std::stod(params.find("r1")->second);
	c.r2 = std::stod(params.find("r2")->second);
	c.d = std::stod(params.find("d")->second);

	// dimensions of system (number of cells)
	uint Nx1 = 100, Nx2 = 100;

	// store entire trajectory
	uint Nstored = tau_tstep + 1;
	std::vector<std::vector<std::shared_ptr<Field> > > tseries(Nstored);

	// create initial configuration and history
	for (int tstep = -tau_tstep; tstep < 1; ++tstep)
	{
		std::vector<std::shared_ptr<Field> > field_ptrs;
		for (uint ifield = 0; ifield < 2; ++ ifield)
		{
			// create concentration fields
			field_ptrs.push_back(std::make_shared<Field>(Nx1, Nx2, 1.0, 1.0, tstep));
			field_ptrs.push_back(std::make_shared<Field>(Nx1, Nx2, 1.0, 1.0, tstep));

			// write configuration to file
			if (tstep == 0)
			{
				std::string filename = datadir + "t" + std::to_string(tstep) + "_f" + std::to_string(ifield) +".dat";
				std::ofstream datafile(filename);
				datafile << field_ptrs.at(ifield)->value << std::endl;
				datafile.close();
			}
		}
		tseries.at(tstep+tau_tstep) = field_ptrs;
	}

	// start timestepping
	for (uint tstep = 1; tstep < Ntstep; ++tstep)
	{
		c.time = static_cast<double>(tstep)*dt;
		unsigned int tprev_idx = (tstep-1+tau_tstep)%Nstored;	// index of previous timestep (within timeseries vector)
		unsigned int tlag_idx = (tstep-1)%Nstored;				// index of timestep used for lag (within timeseries vector)
		std::vector<std::shared_ptr<Field> > field_ptrs;		// temporary container for fields to be inserted into timeseries

		for (uint ifield = 0; ifield < 2; ++ifield)
		{
			// create concentration field (copy prev)
			field_ptrs.push_back(std::make_shared<Field>(*tseries.at(tprev_idx).at(ifield)));

			// zero-flux boundary condition
			field_ptrs.at(ifield)->bound_x1less = field_ptrs.at(ifield)->value.row(0);
			field_ptrs.at(ifield)->bound_x2less = field_ptrs.at(ifield)->value.col(0);
			field_ptrs.at(ifield)->bound_x1more = field_ptrs.at(ifield)->value.row(Nx1-1);
			field_ptrs.at(ifield)->bound_x2more = field_ptrs.at(ifield)->value.col(Nx2-1);

			// advance time, computing rows in parallel
			#pragma omp parallel for
			for (uint x1idx = 0; x1idx < Nx1; ++x1idx)
			{
				for (uint x2idx = 0; x2idx < Nx2; ++x2idx)
				{
					field_ptrs.at(ifield)->value(x1idx,x2idx) += dt*rhs(tseries, tprev_idx, tlag_idx, ifield, x1idx, x2idx, c);
				}
			}

			// write data
			if (((tstep-tstep_start_write)%tstep_intvl_write == 0) && (tstep < tstep_end_write))
			{
				std::string time = std::to_string(static_cast<int>(tstep*dt));
				std::string filename = datadir + "t" + time + "_f" + std::to_string(ifield) +".dat";
				std::ofstream datafile(filename);
				datafile << field_ptrs.at(ifield)->value << std::endl;
				datafile.close();
			}
		}

		tseries.at(tlag_idx) = field_ptrs;

		// display progress
		if (tstep%Ntstep_progress == 0)
		{
			printf("\r ||");
			uint nfullbars = static_cast<uint>(0.5*100.*tstep/Ntstep);
			for (uint pct = 0; pct < nfullbars; ++pct) printf("=");
			for (uint pct = nfullbars; pct < 50; ++pct) printf(" ");
			printf("|| %2.0f%% complete", 100.*tstep/Ntstep); fflush(stdout);
		}
	}
	printf("\r ||"); fflush(stdout);
	for (uint pct = 0; pct < 50; ++pct) printf("=");
	printf("|| 100%% complete\n"); fflush(stdout);

	return 0;
}
