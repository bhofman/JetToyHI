#ifndef thermalAlice_h
#define thermalAlice_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <random>

#include "extraInfo.hh"

//ROOT stuff
#include <TRandom3.h>
#include <TMath.h>
#include "TF1.h"

using namespace fastjet;
using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a thermal event following Boltzman distribution
// Author: M. Verweij
//---------------------------------------------------------------

class thermalAlice {

private :
  TF1              *funcThrm_;

public :
  thermalAlice()
  {
    
    funcThrm_ = new TF1("funcThrm_","TMath::Power(x, [0]-1)*TMath::Exp(-x/[1])", 0.15, 200.); // gamma function
    funcThrm_->SetParNames("alpha", "beta");
    funcThrm_->SetParameters(2.,0.4);
  }
  
  std::vector<fastjet::PseudoJet> createThermalEventAlice() {

    std::vector<fastjet::PseudoJet> particles;

    //double meanN = 2500.;
    //double sigmaN = 500.;
    //std::random_device rd {};
    //std::mt19937 gen {rd()};
    //std::normal_distribution<> d {meanN, sigmaN}; // mean , sigma
    
    unsigned int Nparticles = 2500;//std::round(d(gen));

    for(unsigned int i = 0; i<Nparticles; ++i) {
      //pt from gamam function
      double pt = funcThrm_->GetRandom();
      //random phi
      double phimin = 0.;
      double phimax = TMath::TwoPi();
      double phi = gRandom->Rndm() * (phimax - phimin) + phimin;
      //random rapidity
      double rap = gRandom->Rndm() * 0.9;
      //pion mass
      double mass = 0.1395;
      
      fastjet::PseudoJet p4;
      p4.reset_momentum_PtYPhiM(pt,rap,phi,mass);
      particles.push_back(p4);
    }
    return particles;
  }

};

#endif
