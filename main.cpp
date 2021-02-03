/*
 * main.cpp
 *
 *  Created on: 25.07.2020
 *      Author: Megacrafter127
 */

#include "atoms.hpp"

#include <cstdio>

#include <cstdlib>

#include <set>

#include <getopt.h>

template<typename T> struct Ecmp {
	constants<T> c;
	unsigned Z=1;
	constexpr bool operator()(const orbital<T> &a, const orbital<T> &b) const {
		return a.E(c,Z,0)<b.E(c,Z,0);
	}
};

template<typename T> constexpr const char printRep[]="n: %u\tl: %u\ts: %s\tm_l: %d\t\n\tE: %f\n";
template<> constexpr const char printRep<double>[]="n: %u\tl: %u\ts: %s\tm_l: %d\t\n\tE: %lf\n";
template<> constexpr const char printRep<long double>[]="n: %u\tl: %u\ts: %s\tm_l: %d\t\n\tE: %Lf\n";

template<typename T> constexpr const char printDelt[]="dn: %d\tdl: %d\tds: %d\tdm_l: %d\t\n\tdE: %f\n";
template<> constexpr const char printDelt<double>[]="dn: %d\tdl: %d\tds: %d\tdm_l: %d\t\n\tdE: %lf\n";
template<> constexpr const char printDelt<long double>[]="dn: %d\tdl: %d\tds: %d\tdm_l: %d\t\n\tdE: %Lf\n";


template<typename T> void orbitalConfig(Ecmp<T> cmp, unsigned n) {
	atom<T> at(cmp.Z);
	for(unsigned i=0;i<n;i++) {
		const orbital<T> orb=at.populate(cmp.c);
		if(orb==orbital<T>(-1,-1,0,false)) {
			printf("Electron %u was rejected\n",i);
			break;
		}
		const T E=at.E_i(cmp.c,orb);
		printf(printRep<T>,orb.n,orb.l,orb.s?"1/2":"-1/2",orb.m_l,E);
		bool reseating;
		orbital<T> o2=orb;
		for(reseating=false;;reseating=true) {
			std::pair<orbital<T>,orbital<T>> res=at.reseat(cmp.c);
			if(res.first==res.second) break;
			if(res.first==o2) o2=res.second;
			printf("Reseated:\n");
			const T E1=at.E_i(cmp.c,res.first,res.second),E2=at.E_i(cmp.c,res.second);
			printf(printRep<T>,res.first.n,res.first.l,res.first.s?"1/2":"-1/2",res.first.m_l,E1);
			printf(printRep<T>,res.second.n,res.second.l,res.second.s?"1/2":"-1/2",res.second.m_l,E2);
			printf("delta:\t");
			printf(printDelt<T>,res.first.n-res.second.n,res.first.l-res.second.l,(res.first.s?1:0)-(res.second.s?1:0),res.first.m_l-res.second.m_l,E1-E2);
		}
		if(reseating) {
			printf("Actual Energy:\n");
			const T E2=at.E_i(cmp.c,o2);
			printf(printRep<T>,o2.n,o2.l,o2.s?"1/2":"-1/2",o2.m_l,E2);
			printf("delta:\t");
			printf(printDelt<T>,orb.n-o2.n,orb.l-o2.l,(orb.s?1:0)-(o2.s?1:0),orb.m_l-o2.m_l,E-E2);
		}
	}
	const std::set<orbital<T>> valence=at.valenceOrbitals();
	printf("Valence electrons:\n");
	for(auto &&orb:valence) {
		T e=at.E_i(cmp.c,orb);
		printf(printRep<T>,orb.n,orb.l,(orb.s?"1/2":"-1/2"),orb.m_l,e);
	}
}

template<typename T> void calcEnergy(Ecmp<T> cmp, unsigned n) {
	orbital<T> orb(n,0,0,false);
	for(unsigned i=0;orb.n==n;orb++,i++) {
		printf("%u\t",i);
		printf(printRep<T>,orb.n,orb.l,orb.s?"1/2":"-1/2",orb.m_l,orb.E(cmp.c,cmp.Z,0));
	}
}

enum : int {
	DEFAULT=0,
		ORBITAL_CONFIG=1,
		CALC_ENERGY=2,
		HELP=4,
		PRINT_CONSTANTS=8,
};

#include "constants_io.hpp"

int main(int vc, char **va) {
	int mode=DEFAULT;
	Ecmp<long double> cmp;
	int c;
	static struct option long_options[]={
	                                     {"orbitalConfig", required_argument, &mode, ORBITAL_CONFIG},
	                                     {"calcEnergies", required_argument, &mode, CALC_ENERGY},
	                                     {"help", no_argument, &mode, HELP},
	                                     {"printConstants",no_argument, &mode, PRINT_CONSTANTS}
	};
	unsigned n=0;
	int opt_idx;
	while((c=getopt_long(vc,va,"Z:z:C:c:h",long_options,&opt_idx)) != -1) {
		FILE *f;
		switch(c) {
		case 0:
			switch(long_options[opt_idx].val) {
			case ORBITAL_CONFIG:
			case CALC_ENERGY:
				n=strtoul(optarg,NULL,0);
				break;
			case PRINT_CONSTANTS:
				save_constants(cmp.c,stdout);
				break;
			}
			break;
			break;
			case 'Z':
			case 'z':
				cmp.Z=strtoul(optarg,NULL,0);
				break;
			case 'C':
			case 'c':
				f=fopen(optarg,"r");
				if(!f) {
					fprintf(stderr,"Unable to open file: %s",optarg);
					return 1;
				}
				load_constants(cmp.c,f);
				fclose(f);
				break;
			case 'h':
				mode=HELP;
				break;
		}
	}
	switch(mode) {
	case DEFAULT:
		do{
			printf("\nNumber of electrons: ");
		} while(scanf(" %u",&n)==0);
	case ORBITAL_CONFIG:
		orbitalConfig(cmp,n);
		break;
	case CALC_ENERGY:
		calcEnergy(cmp,n);
		break;
	case HELP:
		printf("Atomaker [options...]\n\nOptions:\n");
		printf("\t-Z -z <number of protons>\n\t\tSpecify how many protons the nucleus should contain.\n\n");
		printf("\t-C -c <file>\n\t\tLoad the physical constants from the given file.\n\n");
		printf("Mode Options:\n");
		printf("\t--help\n\t\tPrint this message.\n\n");
		printf("\t--calcEnergies <n>\n\t\tCalculate the energies of all orbitals of the given principle quantum number.\n\n");
		printf("\t--orbitalConfig <electrons>\n\t\tAttempt to calculate the electron configuration of an atom, by adding 1 electron at a time, and minimizing the total energy level of the atom.\n\n");
		break;
	}
	return 0;
}
