/*
 * orbitals.hpp
 *
 *  Created on: 23.07.2020
 *      Author: Megacrafter127
 */
#ifndef ORBITALS_HPP_
#define ORBITALS_HPP_

#include "constants.hpp"

template<typename T> struct orbital {
	typedef T scalar;
	unsigned n=0; //0<=n //principal quantum number -1
	unsigned l=0; //l<=n //
	int m_l=-l; //-l<=m_l<=l
	bool s=false; //false=-1/2;true=1/2
	typedef enum {
		SHELL,
		SUBSHELL,
		PAIR,
		INDIVIDUAL
	} orbitalGroup;
	constexpr orbital(unsigned n, unsigned l, int m_l, bool s) :
					n(n), l(l), s(s), m_l(m_l) {}
	constexpr orbital()=default;
	constexpr orbital(const orbital&)=default;
	/**
	 * @return combined angular momentum number+1/2
	 */
	constexpr unsigned jph() const {
		return l?l+(s?1:0):1;
	}
	/**
	 * how many proton charges of the nucleus are obscured from this orbital by the other orbital if it is occupied by an electron
	 * @param other - the other orbital
	 * @return the number of proton charges obscured.
	 */
	constexpr scalar S_i(const orbital<scalar> &other) const {
		if(*this==other) return 0;
		if(n<other.n) return 0;
		if(other.n+1<n) return 1;
		if(other.n<n) return l>1?1:0.85;
		if(other.l>l && (other.l!=1 || l!=0)) return 0;
		if(l>1 && l>other.l) return 1;
		return 0.35;
	}
	/**
	 * @param c - the natural constants
	 * @return the absolute angular momentum of this orbital.
	 */
	constexpr scalar L(const constants<scalar> &c) const { //absolute angular momentum
		return c.h_bar*sqrts<scalar>(l*(l+1));
	}
	/**
	 *  Compute 1/(r^3), where r is the average distance between the electron and the center of the nucleus.
	 * @param c - the natural constants
	 * @param Z - the charge of the nucleus
	 * @return 1/(r^3)
	 */
	constexpr scalar invR3(const constants<scalar> &c, scalar Z) const {
		return pow( Z / (c.a_B() * (n+1)) ,3) / (l*(l+0.5)*(l+1));
	}
	/**
	 * @param c - the natural constants
	 * @param Z - the charge of the nucleus
	 * @return strength of the magnetic field in the Z direction induced by this orbital.
	 */
	constexpr scalar B(const constants<scalar> &c, scalar Z) const {
		if(l) return c.mu_0()*c.e*L(c)*m_l*invR3(c,Z) / (4*c.PI*c.m_e);
		else return 0;
	}
	/**
	 * @param c - the natural constants
	 * @param Z - the charge of the nucleus
	 * @return the principal energy level of this orbital.
	 */
	constexpr scalar E_n(const constants<scalar> &c, scalar Z) const {
		return -(c.m_e * pow(c.e,4) * pow(Z,2)) / (2 * pow( (n+1) * c.h() * 2 * c.epsilon_0 ,2));
	}
	/**
	 * @param c - the natural constants
	 * @param Z - the charge of the nucleus
	 * @return the change in energy caused by fine structure splitting.
	 */
	constexpr scalar dE_FS(const constants<scalar> &c, scalar Z) const {
		return ((1/jph()) - (3 / (4*n+4))) * pow(Z*c.alpha(),2) / (n+1);
	}
	/**
	 * @param c - the natural constants
	 * @param B - the external magnetic field in the Z direction
	 * @return the change in energy caused by the Zeeman effect.
	 */
	constexpr scalar dE_mag(const constants<scalar> &c, scalar B) const {
		return -c.h_bar*c.e*m_l*B / (2*c.m_e);
	}
	/**
	 * @param c - the natural constants
	 * @param Z - the charge of the nucleus
	 * @param B - the external magnetic field in the Z direction
	 * @return the total energy level of this orbital.
	 */
	constexpr scalar E(const constants<scalar> &c, scalar Z, scalar B) const {
		return E_n(c,Z) * (1+dE_FS(c,Z)) + dE_mag(c,B);
	}
	constexpr bool operator==(const orbital<scalar> &other) const noexcept {
		return n==other.n && l==other.l && s==other.s && m_l==other.m_l;
	}
	constexpr bool operator!=(const orbital<scalar> &other) const noexcept {
		return !operator==(other);
	}
	constexpr bool operator<(const orbital<scalar> &other) const noexcept {
		if(n!=other.n) return n<other.n;
		else if(l!=other.l) return l<other.l;
		else if(m_l!=other.m_l)return m_l<other.m_l;
		else return s<other.s;
	}
	template<orbitalGroup group> constexpr bool sameGroup(const orbital &other) const noexcept {
		switch(group) {
		default:
		case INDIVIDUAL:
			if(s!=other.s) return false;
		case PAIR:
			if(m_l!=other.m_l) return false;
		case SUBSHELL:
			if(l!=other.l) return false;
		case SHELL:
			return n==other.n;
		};
	}
	template<orbitalGroup group> constexpr orbital lower_limit() const noexcept {
		switch(group) {
		default:
		case INDIVIDUAL:
			return *this;
		case PAIR:
			return orbital(n,l,m_l,false);
		case SUBSHELL:
			return orbital(n,l,-l,false);
		case SHELL:
			return orbital(n,0,0,false);
		}
	}
	template<orbitalGroup group> constexpr orbital upper_limit() const noexcept {
		switch(group) {
		default:
		case INDIVIDUAL:
			return *this;
		case PAIR:
			return orbital(n,l,m_l,true);
		case SUBSHELL:
			return orbital(n,l,l,true);
		case SHELL:
			return orbital(n,n,n,true);
		}
	}
	/**
	 * @return whether or not this orbital actually exists.
	 */
	constexpr bool valid() const noexcept {
		return l<=n && absu(m_l)<=l;
	}
	constexpr orbital &operator++() noexcept {
		s=!s;
		if(!s) {
			m_l++;
			if(absu(m_l)>l) {
				l++;
				if(l>n)  {
					l=0;
					n++;
				}
				m_l=-l;
			}
		}
		return *this;
	}
	constexpr orbital &operator--() noexcept {
		s=!s;
		if(s) {
			m_l--;
			if(absu(m_l)>l) {
				if(l==0) {
					n--;
					l=n;
				} else {
					l--;
				}
				m_l=l;
			}
		}
		return *this;
	}
	constexpr inline orbital operator++(int) noexcept {
		orbital cpy=*this;
		this->operator++();
		return cpy;
	}
	constexpr orbital operator--(int) noexcept {
		orbital cpy=*this;
		this->operator--();
		return cpy;
	}
	
	template<typename U> constexpr operator orbital<U>() const noexcept {
		return orbital<U>(n,l,s,m_l);
	}
};

#endif /* ORBITALS_HPP_ */
