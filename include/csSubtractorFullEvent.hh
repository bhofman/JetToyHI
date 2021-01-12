#ifndef csSubtractorFullEvent_h
#define csSubtractorFullEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the full event constituent subtraction
// Author: M. Verweij
//---------------------------------------------------------------

class csSubtractorFullEvent {

private :
  double alpha_;
  double rParam_;
  double ghostArea_;
  double ghostRapMax_;
  double rho_;
  double rhom_;
  std::vector<fastjet::PseudoJet> fjInputs_;

  contrib::ConstituentSubtractor subtractor_;

  
public :
  csSubtractorFullEvent(double alpha = 1., double rParam = 0.25, double ghostArea = 0.005, double ghostRapMax = 3.0) :
    alpha_(alpha),
    rParam_(rParam),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    rho_(-1),
    rhom_(-1)
  {
    //init constituent subtractor
    subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
    subtractor_.set_max_distance(rParam_); //free parameter for the maximal allowed distance between particle i and ghost k
    subtractor_.set_alpha(alpha_); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
    subtractor_.set_scale_fourmomentum(); //Keep rapidity and pseudo-rapidity fixed (scale fourmomentum). Recommended - observed better performance than the mass correction. Use: subtractor.set_scale_fourmomentum();
    //subtractor_.set_do_mass_subtraction(); // No function argument needed, was (true);
  }

  void setAlpha(double a)     { alpha_ = a; }
  void setRParam(double r)    { rParam_ = r; }
  void setGhostArea(double a) { ghostArea_ = a; }

  void setRho(double r)       { rho_ = r; }
  void setRhom(double r)      { rhom_ = r; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoM() const { return rhom_; }
  
  std::vector<fastjet::PseudoJet> doSubtractionFullEvent() {

    fastjet::GridMedianBackgroundEstimator bkgd_estimator(3.,0.5);  //max_eta, grid spacing
    bkgd_estimator.set_particles(fjInputs_);  

    if(rho_<0.) {    
      cout << "rho < 0; with rho = " << rho_ <<"rhom = " << rhom_ << endl;
      rho_ = bkgd_estimator.rho();
      rhom_ = bkgd_estimator.rho_m();
      cout << "Now using rho = "<< rho_ << " using rhom = "<< rhom_ << endl;
      } else {
      cout << "External rho = " << rho_ << " External rhom = " << rhom_ << endl;
      }

    subtractor_ = contrib::ConstituentSubtractor(rho_,rhom_,alpha_,rParam_,contrib::ConstituentSubtractor::deltaR);
    subtractor_.set_background_estimator(&bkgd_estimator);
    subtractor_.set_max_eta(3.);
    subtractor_.initialize();
    
    //cout << subtractor_.description() << endl; // print info (optional)
    std::vector<fastjet::PseudoJet> corrected_event = subtractor_.subtract_event(fjInputs_);
    return corrected_event;
  }
};

#endif
