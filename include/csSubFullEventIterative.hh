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
	double max_eta_


public:
	csSubFullEventIterative(vector<double> alpha = {1.,1.}, vector<double> rParam = {0.25,0.25},  double ghostArea = 0.005,                ) : 
		alpha_(alpha),
		rParam_(rParam),
		ghostArea_(ghostArea),
	{
		subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
		subtractor_.set_parameters(max_distances_,alphas_);
		subtractor_.set_ghost_removal(true);
		subtractor_.set_ghost_area(ghostArea_);
	}

	void setMaxEta(double max_eta)     	{ max_eta_ = max_eta; subtractor_.set_max_eta(max_eta); }
	
	void setAlpha(vector<double> a)     { alpha_ = a; }
	void setRParam(vector<double >r)    { rParam_ = r; }

	void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }

	void setBackground() {
		GridMedianBackgroundEstimator bge_rho(max_eta_,0.5);
		bkgd_estimator.set_particles(fjInputs_);  
		subtractor_.set_background_estimator(&bge_rho); 
	}



	std::vector<fastjet::PseudoJet> doSubtractionFullEvent() {
		subtractor_.initialize();
		std::vector<fastjet::PseudoJet> corrected_event = subtractor_.subtract_event(fjInputs_);
		return corrected_event;
  	}






















};	
#endif