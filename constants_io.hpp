/*
 * constants_io.hpp
 *
 *  Created on: 27.07.2020
 *      Author: Megacrafter127
 */
#ifndef CONSTANTS_IO_HPP_
#define CONSTANTS_IO_HPP_

#include "constants.hpp"

#include <cstdio>

template<typename T> constexpr static const char constantFmt[]="\nm_e: %f\nm_p: %f\nm_n: %f\nh_bar: %f\nepsilon_0: %f\ne: %f\nc: %f";
template<> constexpr const char constantFmt<double>[]="\nm_e: %lf\nm_p: %lf\nm_n: %lf\nh_bar: %lf\nepsilon_0: %lf\ne: %lf\nc: %lf";
template<> constexpr const char constantFmt<long double>[]="\nm_e: %Lf\nm_p: %Lf\nm_n: %Lf\nh_bar: %Lf\nepsilon_0: %Lf\ne: %Lf\nc: %Lf";

template<typename T> int load_constants(constants<T> &c,FILE *f) {
	return fscanf(f,constantFmt<T>,&c.m_e,&c.m_p,&c.m_n,&c.h_bar,&c.epsilon_0,&c.e,&c.c);
}
template<typename T> int save_constants(const constants<T> &c,FILE *f) {
	return fprintf(f,constantFmt<T>,c.m_e,c.m_p,c.m_n,c.h_bar,c.epsilon_0,c.e,c.c);
}

#endif /* CONSTANTS_IO_HPP_ */
