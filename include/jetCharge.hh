#ifndef __JetCharge_HH__
#define __JetCharge_HH__

//------------------------------------------------------------------------
/// jet charge
///
/// This is defined as in 

class JetCharge {
public:
  /// default ctor
  JetCharge(double kappa = 0.5, double ptmin = -1.) :
    _kappa(kappa),
    _ptmin(ptmin)
   {}

  /// compute the function
  virtual double result(const fastjet::PseudoJet &jet) const {
    // check the jet is appropriate for computation
    if (!jet.has_constituents()) {
      Printf("Jet charge calculation can only be applied on jets for which the constituents are known.");
      return -999.;
    }
    vector<fastjet::PseudoJet> constits = jet.constituents();
    double sumcharge = 0.;
    double jetPt = jet.perp();
    if (jetPt == 0) return -999.;

    for(fastjet::PseudoJet p : constits) {
      if(p.perp()<_ptmin) continue;
      const double & ch = p.user_info<PU14>().charge(); //three_charge()
      double zFrac = p.perp()/jetPt;
      sumcharge += ch*std::pow(zFrac,_kappa);

    }
    return sumcharge;
  }
  
protected:
  double _kappa;
  double _ptmin;
};

#endif