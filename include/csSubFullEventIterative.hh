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
	vector<double> alpha_;
	vector<double> rParam_;
	double ghostArea_;
	double ghostRapMax_;
	double rho_;
	double rhom_;
	double max_eta_;

	std::vector<fastjet::PseudoJet> fjInputs_;
	contrib::IterativeConstituentSubtractor subtractor_;

public:
	csSubFullEventIterative(vector<double> alpha = {0.,0.,0.}, vector<double> rParam = {0.1,0.1,.1},  double ghostArea = 0.005,	double ghostRapMax = 3.0) :
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
		subtractor_.set_parameters(rParam_,alpha_);
		subtractor_.set_ghost_area(ghostArea_);
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

	void setBackgroundGrid() { // Set the background rho using a grid instead of jets
		subtractor_.initialize();

		GridMedianBackgroundEstimator bge_rho(max_eta_,0.5);
		
		bge_rho.set_particles(fjInputs_);  
		subtractor_.set_background_estimator(&bge_rho); 
		rho_ = bge_rho.rho();
		rhom_ = bge_rho.rho_m();
	}

	void setBackground() {
		subtractor_.initialize();

		AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(5.)); // Ghost rho should go atleasy 2R beyond jet region
		JetDefinition jet_def_bkgd(kt_algorithm, 0.4);
   		Selector selector = SelectorAbsRapMax(3.) * (!SelectorNHardest(2));
		JetMedianBackgroundEstimator bge_rho(selector, jet_def_bkgd, area_def_bkgd);

		bge_rho.set_particles(fjInputs_);  
		subtractor_.set_background_estimator(&bge_rho); 
		rho_ = bge_rho.rho();
		rhom_ = bge_rho.rho_m();

		subtractor_.set_background_estimator(&bge_rho);
		subtractor_.set_common_bge_for_rho_and_rhom(true);
	}

	std::vector<fastjet::PseudoJet> doSubtractionFullEvent() {
		//cout << subtractor_.description() << endl;

		if(rho_<0.) {  cout << "Rho < 0 " << endl ; }
		//cout << "rho_: "<<rho_<<" rhom_: "<<rhom_<<endl;

		std::vector<fastjet::PseudoJet> corrected_event = subtractor_.subtract_event(fjInputs_);
		return corrected_event;
  	}
};	
#endif