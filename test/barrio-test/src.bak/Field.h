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
		uint 			m_Nx1;					// number of cells in x direction
		uint 			m_Nx2;					// number of cells in y direction
		arma::mat 		m_value;				// value of concentration field at each (x,y) point
		arma::rowvec 	m_bound_x1less;			// boundary values on x1less (left) edge
		arma::rowvec 	m_bound_x1more;			// on x1more (right) edge
		arma::colvec 	m_bound_x2less;			// on x2less (lower) edge
		arma::colvec 	m_bound_x2more;			// on x2more (upper) edge
		double 			m_ax1;					// lattice const x direction
		double 			m_ax2;					// lattice const y direction
		
		std::mt19937 	rng;					// prng engine

		Field(uint Nx1, uint Nx2, double ax1, double ax2)
		{
			m_Nx1 = Nx1;
			m_Nx2 = Nx2;
			m_ax1 = ax1;
			m_ax2 = ax2;
			m_value.resize(m_Nx1, m_Nx2);
			m_bound_x1less.resize(m_Nx2);
			m_bound_x1more.resize(m_Nx2);
			m_bound_x2less.resize(m_Nx1);
			m_bound_x2more.resize(m_Nx1);

			// initial condition
			for (uint icellx = 0; icellx < m_Nx1; ++icellx)
			{
				for (uint jcelly = 0; jcelly < m_Nx2; ++jcelly)
				{
					m_value(icellx,jcelly) = 1.00;
				}
			}

			//std::uniform_real_distribution<double> dist(0.00, 1.00);
			//m_value.imbue([&](){ return dist(rng); });
			//m_bound_x2more.imbue([&](){ return dist(rng); });
			//m_bound_x1more.imbue([&](){ return dist(rng); });
			//m_bound_x2less.imbue([&](){ return dist(rng); });
			//m_bound_x1less.imbue([&](){ return dist(rng); });
			//std::cout << m_bound_x1less << std::endl;
			//std::cout << m_value << std::endl;
			//std::cout << m_bound_x1more << std::endl;
			//std::cout << m_bound_x2more << std::endl;
			//std::cout << m_bound_x2less << std::endl;
		}

		double& at(uint x1idx, uint x2idx)
		{
			return m_value(x1idx, x2idx);
		}

		double& x1less(uint x1idx, uint x2idx)
		{
			return (x1idx == 0) ? m_bound_x1less(x2idx) : m_value(x1idx-1,x2idx);
		}

		double& x1more(uint x1idx, uint x2idx)
		{
			return (x1idx == (m_Nx1-1)) ? m_bound_x1more(x2idx) : m_value(x1idx+1,x2idx);
		}

		double& x2less(uint x1idx, uint x2idx)
		{
			return (x2idx == 0) ? m_bound_x2less(x1idx) : m_value(x1idx,x2idx-1);
		}

		double& x2more(uint x1idx, uint x2idx)
		{
			return (x2idx == (m_Nx2-1)) ? m_bound_x2more(x1idx) : m_value(x1idx,x2idx+1);
		}

		void add_value(arma::mat addl)
		{
			m_value += addl;
		}

};

#endif
