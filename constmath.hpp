/*
 * constmath.hpp
 *
 *  Created on: 23.07.2020
 *      Author: Megacrafter127
 */
#ifndef CONSTMATH_HPP_
#define CONSTMATH_HPP_

#include <cmath>

template<typename T> constexpr T pow(T a, unsigned p) { //a^p
	return p? (p%2? (a*pow(a,p-1)) : pow(a,p/2)*pow(a,p/2)) :1;
}

template<typename T> constexpr T max(T a, T b) {
	return a<b?b:a;
}

template<typename T> constexpr T sqrts(T a) {
	return sqrtl(a);
}
template<> constexpr float sqrts<float>(float a) {
	return sqrtf(a);
}
template<> constexpr double sqrts<double>(double a) {
	return sqrt(a);
}

template<typename T> constexpr T cbrts(T a) {
	return cbrtl(a);
}
template<> constexpr float cbrts<float>(float a) {
	return cbrtf(a);
}
template<> constexpr double cbrts<double>(double a) {
	return cbrt(a);
}

constexpr unsigned absu(int a) {
	return a<0?-a:a;
}

#endif /* CONSTMATH_HPP_ */
