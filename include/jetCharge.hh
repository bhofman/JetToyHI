#ifndef jetCharge_h
#define jetCharge_h

#include "TPDGCode.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class calculates the jet charge
// ** WORK IN PROGRESS **  Need some work on adding PDG code from pu14
//---------------------------------------------------------------

class jetCharge {

private:
	std::vector<double> jetcharge_ ;
	double num_ ;
	double den_ ;
	double kappa_ ; 

	std::vector<fastjet::PseudoJet> fjInputs_;

	std::map<int, int> PDGCharges = {
	    {211, 1},    // pi+
	    {-211, -1},  // pi-
	    {321, 1},    // K+
	    {-321, -1},  // K-
	    {2212, 1},   // p
	    {-2212, -1}, // p-
	    {11, -1},    // e-
	    {-11, 1},    // e+
	    {13, -1},    // mu-
	    {-13, 1},    // mu+
	    {22, 0},     // gamma
	    {111, 0},    // pi0
	    {130, 0},    // K0L
	    {2112, 0},   // n
	    {-2112, 0},  // nbar
	    {311, 0},    // K0
	    {12, 0},     // nue
	    {-12, 0},    // nuebar
	    {14, 0},     // numu
	    {-14, 0},    // numubar
	    {16, 0},     // nutau
	    {-16, 0}     // nutaubar
	};

public:
	jetCharge(std::vector<fastjet::PseudoJet> fjInputs, double kappa):
		fjInputs_(fjInputs),
		kappa_(kappa)
	{
	}

	std::vector<double> getCharge()
	{
		jetcharge_.reserve(fjInputs_.size());
		
		for(fastjet::PseudoJet& jet : fjInputs_) {

			if (!jet.has_constituents()) {
			    std::cout<<"Can not calculate jet charge if jet has no constituents!"<<std::endl;
			    jetcharge_.push_back(-99999.);
			}
			for (auto part : jet.constituents())
			{
			    den_ = den_ + part.perp(); // pow(part.perp(), kappa);
			    int _pdgid = part.user_info().pdg_id;

			    std::cout<<_pdgid<<std::endl;

			    if (PDGCharges.count(_pdgid) == 1)
			    {
			        num_ += PDGCharges[_pdgid] * pow(part.perp(), kappa_);
			    }
			    else
			    {
			        cout << " Charge for PDG id " << _pdgid << " is not defined, using 0." << endl;
			    }
			}

			den_ = pow(den_, kappa_);
			jetcharge_.push_back(num_/den_);
		}

		return jetcharge_;
	}
};	
#endif