/*
 * constants.hpp
 *
 *  Created on: 23.07.2020
 *      Author: Megacrafter127
 */
#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#include "constmath.hpp"
#include <limits>

template<typename T> struct constants {
	typedef T scalar;
	constexpr static const scalar PI = 0x3.243F6A8885A308D31319p0;
	constexpr static const unsigned infint = std::numeric_limits<unsigned>::max();
	scalar m_e=4.18493492472587e-3, //mass of electron
		m_p=7.68417979427673, //mass of proton
		m_n=7.69477146703062, //mass of neutron
		h_bar=1, //reduced Plank constant
		epsilon_0=0.25/PI, //electric field constant
		e=1, //elementary electic charge
		c=1; //speed of light in vacuum
	constexpr scalar h() const { //Plank constant
		return 2*PI*h_bar;
	}
	constexpr scalar a_B() const { //Bohr radius
		return 4*PI*epsilon_0*pow(h_bar,2) / (m_e*pow(e,2));
	}
	constexpr scalar mu_0() const { //magnetic field constant
		return 1/(pow(c,2)*epsilon_0);
	}
	constexpr scalar mu_K() const { //nuclear magneton
		return h_bar*e/(2*m_p);
	}
	constexpr scalar alpha() const { //fine structure constant
		return pow(e,2) / (2*c*epsilon_0*h());
	}
	constexpr scalar R_y() const { //Rydberg energy
		return pow(alpha()*c,2)*m_e / 2;
	}
};

#endif /* CONSTANTS_HPP_ */
