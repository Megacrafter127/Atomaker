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
#include <map>

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
	
	/**
	 * Compute the total shielding effect for a given orbital.
	 * @param orb - the orbital to calculate the shielding effect for
	 * @param exclude - the orbital to exclude from the calculations in the shielding effect
	 * @return total shielding effect for the given orbital, in proton charges
	 */
	constexpr T S(const orbital<T> &orb, const orbital<T> &exclude=orbital<T>(-1,-1,0,false)) const {
		T ret=0;
		for(auto &&orbital:orbitals) {
			 if(orbital!=exclude) ret+=orb.S_i(orbital);
		}
		return ret;
	}
	
	/**
	 * Compute total magnetic field in the Z direction.
	 * @param c - the natural constants
	 * @param exclude - an optional orbital to exclude from the calculations
	 * @return the total magnetic field in the Z direction.
	 */
	constexpr T B(const constants<T> &c, const orbital<T> &exclude=orbital<T>(-1,-1,0,false)) const {
		T ret=0;
		for(auto &&orb:orbitals) {
			if(orb!=exclude) ret+=orb.B(c,Z-S(orb,exclude));
		}
		return ret;
	}
	
	/**
	 * Compute ionization energy of a given orbital.
	 * @param c - the natural constants
	 * @param orb - the orbital to calculate the ionization energy for
	 * @param exclude - an optional orbital to exclude from the calculations
	 * @return the ionization energy for the given orbital.
	 */
	constexpr T E_i(const constants<T> &c, const orbital<T> &orb, const orbital<T> &exclude=orbital<T>(-1,-1,0,false)) const {
		const T B=0;//this->B(c,orb);
		return orb.E(c,Z-S(orb,exclude),B);
	}
	
	/**
	 * Compute the total binding energy of all electrons.
	 * @param c - the natural constants
	 * @return the total binding energy of all electrons.
	 */
	constexpr T E(const constants<T> &c) const {
		T ret=0;
		for(auto &&orb:orbitals) {
			ret+=E_i(c,orb);
		}
		return ret;
	}
	
	/**
	 * Recompute the lowest energy orbital for the given electron.
	 * @param c - the natural constants
	 * @param old - the electron in question
	 * @return the lowest energy orbital for the given electron.
	 */
	constexpr orbital<T> reseatSpot(const constants<T> &c, const orbital<T> &old) const {
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
	
	/**
	 * Test if any electron's energy levels need recomputing, and recompute it.
	 * @param c - the natural constants
	 * @return a pair consisting of the previous orbital and the new orbital, if a change occurred. otherwise the exact return value isn't well defined, except that comparing for equality will be true.
	 */
	std::pair<orbital<T>,orbital<T>> reseat(const constants<T> &c) {
		for(auto i=orbitals.crbegin();i!=orbitals.crend();i++) {
			const orbital<T> res=reseatSpot(c,*i);
			if(res!=*i) {
				const orbital<T> old=*i;
				orbitals.erase(--i.base());
				orbitals.insert(res);
				return std::make_pair(old,res);
			}
		}
		return std::make_pair(orbital<T>(-1,-1,0,false),orbital<T>(-1,-1,0,false));
	}
	
	/**
	 * Add an electron to the atom.
	 * A.k.a: find the unoccupied orbital with the highest binding energy, and insert an electron into it.
	 * @param c - the natural constants
	 * @return the orbital the electron was inserted into.
	 */
	orbital<T> populate(const constants<T> &c) { //
		constexpr orbital<T> free=orbital<T>(-1,-1,0,false);
		const orbital<T> ret=reseatSpot(c,free);
		if(ret==free) return free;
		orbitals.insert(ret);
		return ret;
	}
	
	/**
	 * @return the set of valence orbitals for this atom.
	 */
	constexpr std::set<orbital<T>> valenceOrbitals() const {
		std::set<orbital<T>> ret;
		std::map<unsigned,unsigned> max;
		for(auto i=orbitals.crbegin();i!=orbitals.crend();i++) {
			if(max[i->l]==0) {
				max[i->l]=i->n+1;
			}
			if(max[i->l]==i->n+1) {
				ret.insert(*i);
			}
		}
		return ret;
	}
	
	/**
	 * Calculates the minimum required energy to ionize this atom to the desired degree.
	 * @param c - the natural constants
	 * @param level - how many electrons the atom needs to be stripped of
	 * @return the minimum required energy to remove the desired number of electrons.
	 */
	constexpr T ionizationEnergy(const constants<T> &c, unsigned level) const {
		std::set<orbital<T>> removed;
		T ret=0;
		for(unsigned c=0;c<level && c<orbitals.size();c++) {
			const orbital<T> min;
			T minE=-INFINITY;
			for(auto i=orbitals.crbegin();i!=orbitals.crend();i++) {
				if(removed.count(*i)) continue;
				T E=E_i(c,*i);
				if(E>minE) {
					minE=E;
					min=*i;
				}
			}
			removed.insert(min);
			ret+=minE;
		}
		return -ret;
	}
};



#endif /* ATOMS_HPP_ */
