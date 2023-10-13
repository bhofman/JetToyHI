#ifndef __JetChargeDynamical_HH__
#define __JetChargeDynamical_HH__
#include <iostream>
#include <TRandom3.h>
//------------------------------------------------------------------------
/// Dynamical jet charge
///
/// This is defined as in https://arxiv.org/pdf/2101.04304.pdf

class JetChargeDynamical {
public:
  /// default ctor
  JetChargeDynamical(double Xi = 0.3, double kappaLower = 1.0,double kappaHigher = 0.3, double ptmin = -1.) :
    _Xi(Xi),
    _kappaLower(kappaLower),
    _kappaHigher(kappaHigher),
    _ptmin(ptmin)
   {}

  /// compute the function
  virtual double result(const fastjet::PseudoJet &jet) const {
    // check the jet is appropriate for computation
    if (!jet.has_constituents()) {
      Printf("Dynamical Jet charge calculation can only be applied on jets for which the constituents are known.");
      return -999.;
    }
    vector<fastjet::PseudoJet> constits = jet.constituents();
    double sumcharge = 0;
    double jetPt = jet.perp();
    if (jetPt == 0) return -999.;

    for(fastjet::PseudoJet p : constits) {
      if(p.perp()<_ptmin) continue;
      double zFrac = p.perp()/jetPt;
      double ch = p.user_info<PU14>().charge();
      if(zFrac <= _Xi) sumcharge += ch*std::pow(zFrac,_kappaLower);
      if(zFrac > _Xi) sumcharge += ch*std::pow(zFrac,_kappaHigher);
    }
    return sumcharge;
  }
  
protected:
  double _Xi;
  double _kappaLower;
  double _kappaHigher;
  double _ptmin;
};

#endif