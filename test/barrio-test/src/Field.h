#ifndef FIELD_H
#define FIELD_H

#include <armadillo>
#include <iostream>
#include <random>
#include <string>
#include <utility>


// note that the directions are those that you would use if the entire system was rotated 90 degrees clockwise

class Field
{
	public:
		uint 			Nx1;					// number of cells in x direction
		uint 			Nx2;					// number of cells in y direction
		arma::mat 		value;					// value of concentration field at each (x,y) point
		arma::rowvec 	bound_x1less;			// boundary values on x1less (left) edge
		arma::rowvec 	bound_x1more;			// on x1more (right) edge
		arma::colvec 	bound_x2less;			// on x2less (lower) edge
		arma::colvec 	bound_x2more;			// on x2more (upper) edge
		double 			ax1;					// lattice const x direction
		double 			ax2;					// lattice const y direction
		
		std::mt19937 	rng;					// prng engine

		Field(uint Nx1, uint Nx2, double ax1, double ax2)
		{
			this->Nx1 = Nx1;
			this->Nx2 = Nx2;
			this->ax1 = ax1;
			this->ax2 = ax2;
			this->value.resize(this->Nx1, this->Nx2);
			this->bound_x1less.resize(this->Nx2);
			this->bound_x1more.resize(this->Nx2);
			this->bound_x2less.resize(this->Nx1);
			this->bound_x2more.resize(this->Nx1);

			// initial condition
			//for (uint icellx = 0; icellx < this->Nx1; ++icellx)
			//{
				//for (uint jcelly = 0; jcelly < this->Nx2; ++jcelly)
				//{
					//this->value(icellx,jcelly) = 0.00;
				//}
			//}

			//for (uint icellx = static_cast<uint>(this->Nx1/4); icellx < 3*static_cast<uint>(this->Nx1/4); ++icellx)
			//{
				//for (uint jcelly = static_cast<uint>(this->Nx2/4); jcelly < 3*static_cast<uint>(this->Nx2/4); ++jcelly)
				//{
					//this->value(icellx,jcelly) = 1.00;
				//}
			//}
			std::uniform_real_distribution<double> dist(-0.50, 0.50);
			this->value.imbue([&](){ return dist(rng); });

			//this->bound_x2more.imbue([&](){ return dist(rng); });
			//this->bound_x1more.imbue([&](){ return dist(rng); });
			//this->bound_x2less.imbue([&](){ return dist(rng); });
			//this->bound_x1less.imbue([&](){ return dist(rng); });
			//std::cout << this->bound_x1less << std::endl;
			//std::cout << this->value << std::endl;
			//std::cout << this->bound_x1more << std::endl;
			//std::cout << this->bound_x2more << std::endl;
			//std::cout << this->bound_x2less << std::endl;
		}

		double& at(uint x1idx, uint x2idx)
		{
			return this->value(x1idx, x2idx);
		}

		double& x1less(uint x1idx, uint x2idx)
		{
			return (x1idx == 0) ? this->bound_x1less(x2idx) : this->value(x1idx-1,x2idx);
		}

		double& x1more(uint x1idx, uint x2idx)
		{
			return (x1idx == (this->Nx1-1)) ? this->bound_x1more(x2idx) : this->value(x1idx+1,x2idx);
		}

		double& x2less(uint x1idx, uint x2idx)
		{
			return (x2idx == 0) ? this->bound_x2less(x1idx) : this->value(x1idx,x2idx-1);
		}

		double& x2more(uint x1idx, uint x2idx)
		{
			return (x2idx == (this->Nx2-1)) ? this->bound_x2more(x1idx) : this->value(x1idx,x2idx+1);
		}

		void add_value(arma::mat addl)
		{
			this->value += addl;
		}

};

#endif
