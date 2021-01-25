#ifndef csSubFullEventIterative_h
#define csSubFullEventIterative_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/IterativeConstituentSubtractor.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the full event constituent subtraction iteratively
// ** WORK IN PROGRESS **
//---------------------------------------------------------------

class csSubFullEventIterative {

private:
	std::vector<fastjet::PseudoJet> fjInputs_;
	contrib::IterativeConstituentSubtractor subtractor_;
	// Parameters for csSubs
	vector<double> rParam_;
	vector<double> alpha_;
	double ghostArea_;
	double max_eta_;
	double ghostRapMax_;
	double rho_;
	double rhom_;

public:
	csSubFullEventIterative(vector<double> alpha = {1.,1.}, vector<double> rParam = {0.25,0.25},  double ghostArea = 0.005,	double ghostRapMax = 3.0) :
		alpha_(alpha),
		rParam_(rParam),
		ghostArea_(ghostArea),
		ghostRapMax_(ghostRapMax),
		rho_(-1),
		rhom_(-1)
	{
		subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
		subtractor_.set_ghost_removal(true);
		subtractor_.set_scale_fourmomentum(); 
	}

	void setAlpha(vector<double> a)     { alpha_ = a; }
	void setRParam(vector<double >r)    { rParam_ = r; }
	void setGhostArea(double a) { ghostArea_ = a; }

	void setRho(double r)       { rho_ = r; }
	void setRhom(double r)      { rhom_ = r; }

	void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

	double getRho()  const { return rho_; }
	double getRhoM() const { return rhom_; }

	void setMaxEta(double max_eta)     	{ max_eta_ = max_eta; subtractor_.set_max_eta(max_eta); }

	void setBackground() {
		GridMedianBackgroundEstimator bge_rho(max_eta_,0.5);
		bkgd_estimator.set_particles(fjInputs_);  
		subtractor_.set_background_estimator(&bge_rho); 
	}

	std::vector<fastjet::PseudoJet> doSubtractionFullEvent() {
		subtractor_.set_ghost_area(ghostArea_);
		subtractor_.set_parameters(max_distances_,alphas_);
		subtractor_.initialize();

		if(rho_<0.) {  cout << "Rho < 0 " << endl ; }

		std::vector<fastjet::PseudoJet> corrected_event = subtractor_.subtract_event(fjInputs_);
		return corrected_event;
  	}
};	
#endif