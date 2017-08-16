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
		
		int 			t;

		Field(uint Nx1, uint Nx2, double ax1, double ax2, int t = 0)
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
			this->t = t;

			//for (unsigned int i = 0; i < Nx1; ++i)
				//for (unsigned int j = 0; j < Nx2; ++j)
					//this->value(i,j) = cos(0.45*i)

			std::uniform_real_distribution<double> dist(-0.50, 0.50);
			std::random_device rd;					// seed init
			std::mt19937 rng(rd());					// prng engine
			this->value.imbue([&](){ return 0.01*dist(rng); });
		}

		Field(const Field& f)
		{
			this->Nx1 = f.Nx1;
			this->Nx2 = f.Nx2;
			this->ax1 = f.ax1;
			this->ax2 = f.ax2;
			this->value	= f.value;
			this->bound_x1less.resize(this->Nx2);
			this->bound_x1more.resize(this->Nx2);
			this->bound_x2less.resize(this->Nx1);
			this->bound_x2more.resize(this->Nx1);
			this->t = f.t + 1;
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
