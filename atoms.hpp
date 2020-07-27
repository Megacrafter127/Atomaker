/*
 * atoms.hpp
 *
 *  Created on: 26.07.2020
 *      Author: Megacrafter127
 */
#ifndef ATOMS_HPP_
#define ATOMS_HPP_

#include "orbitals.hpp"

#include <set>

template<typename T> struct atom {
	const unsigned Z; //number of protons in nucleus
	std::set<orbital<T>> orbitals; //the set of electrons in orbitals
	constexpr atom(unsigned Z) : Z(Z),orbitals() {}
	constexpr atom()=delete;
	constexpr atom(const atom &a)=default;
	constexpr atom(const atom &a, const std::pair<orbital<T>,orbital<T>> &swap) : Z(a.Z),orbitals(a.orbitals) { //create copy of an atom where 1 electron is moved to a different orbital
		orbitals.erase(swap.first);
		orbitals.insert(swap.second);
	}
	
	constexpr T S(const orbital<T> &orb, const orbital<T> &exclude=orbital<T>(-1,-1,false,0)) const { //compute total shielding effect for a given orbital. an orbital can be excluded from the computation
		T ret=0;
		for(auto i=orbitals.cbegin();i!=orbitals.cend();i++) {
			 if(*i!=exclude) ret+=orb.S_i(*i);
		}
		return ret;
	}
	
	constexpr T B(const constants<T> &c, const orbital<T> &exclude=orbital<T>(-1,-1,false,0)) const { //compute total magnetic field in the Z direction
		T ret=0;
		for(auto i=orbitals.cbegin();i!=orbitals.cend();i++) {
			if(*i!=exclude) ret+=i->B(c,Z-S(*i,exclude));
		}
		return ret;
	}
	
	constexpr T E_i(const constants<T> &c, const orbital<T> &orb, const orbital<T> &exclude=orbital<T>(-1,-1,false,0)) const { //compute energy level of a given orbital. an orbital can be excluded from affecting the calculations
		const T B=0;//this->B(c,orb);
		return orb.E(c,Z-S(orb,exclude),B);
	}
	
	constexpr T E(const constants<T> &c) const { //compute total energy level of electrons
		T ret=0;
		for(auto i=orbitals.cbegin();i!=orbitals.cend();i++) {
			ret+=E_i(c,*i);
		}
		return ret;
	}
	
	constexpr orbital<T> reseatSpot(const constants<T> &c, const orbital<T> &old) const { //recompute the lowest energy orbital for the given electron
		orbital<T> minO=old;
		T minE=E(c);
		bool hit=true;
		for(orbital<T> orb;;orb++) {
			if(orb.l==0 && !orb.s) {
				if(hit) hit=false;
				else break;
			}
			if(orbitals.count(orb)) { //pauli principle forbids 2 electrons from occupying the same orbital
				hit=true;
				continue;
			}
			const T E=atom(*this,std::make_pair(old,orb)).E(c);
			if(minE>E) {
				hit=true;
				minO=orb;
				minE=E;
			}			
		}
		return minO;
	}
	
	std::pair<orbital<T>,orbital<T>> reseat(const constants<T> &c) { //test if any electron's energy levels need recomputing, and recompute it
		for(auto i=orbitals.crbegin();i!=orbitals.crend();i++) {
			const orbital<T> res=reseatSpot(c,*i);
			if(res!=*i) {
				const orbital<T> old=*i;
				orbitals.erase(--i.base());
				orbitals.insert(res);
				return std::make_pair(old,res);
			}
		}
		return std::make_pair(orbital<T>(-1,-1,false,0),orbital<T>(-1,-1,false,0));
	}
	
	orbital<T> populate(const constants<T> &c) { //add an electron to the atom
		constexpr orbital<T> free=orbital<T>(-1,-1,false,0);
		const orbital<T> ret=reseatSpot(c,free);
		if(ret==free) return free;
		orbitals.insert(ret);
		return ret;
	}
};



#endif /* ATOMS_HPP_ */
