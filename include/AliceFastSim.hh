#ifndef __AliceFastSim_HH__
#define __AliceFastSim_HH__

#include "TTree.h"
#include "TFile.h"
#include <TRandom.h>
#include <TH1D.h>

//---------------------------------------------------------------
// Description
// This class runs a fast simulation of the ALICE detector response
// Author: B. Hofman
//---------------------------------------------------------------

class AliceFastSim {
public:
  /// default ctor
  AliceFastSim(const char* filename = "include/tr_eff.root", const char* histname = "tr_eff"){
    file = new TFile(filename);
    histogram = (TH1D*)file->Get(histname);

    lowerLimit = histogram->GetXaxis()->GetXmin();
    upperLimit = histogram->GetXaxis()->GetXmax();
  }
  /// default dtor
  ~AliceFastSim() {
    file->Close();
  }

  void setInputEvent(vector<PseudoJet> event) {
    _fullEvent = event;
  }

  /// Returns event with detector acceptance applied
  /// Applies charged particles with pT > 0.15 GeV/c and |eta| < 0.9 cut
  /// This corresponds to the 'Truth' level MC
  vector<PseudoJet> AliceAcceptance() {
    _truthEvent.clear();
    for(fastjet::PseudoJet particles : _fullEvent) {
        if(particles.perp()>=0.15 && particles.eta()>=-0.9 && particles.eta()<=0.9 && particles.user_info<PU14>().charge()!=0) {
            _truthEvent.push_back(particles);
        }
    }
    return _truthEvent;
  }

  /// Return event with detector response applied
  /// Applies tracking efficiency and detector resolution to the truth level event
  /// This corresponds to the 'Detector' level MC
  vector<PseudoJet> AliceDetector() {
    _detectorEvent.clear();
    if (_truthEvent.size() == 0) {
        std::cout << "Warning: AliceAcceptance() has not been run. Running now..." << endl;
        AliceAcceptance();
    }
    
    for(fastjet::PseudoJet particles : _truthEvent) {
        double _pt = particles.perp();
        double randomNumber = randomGenerator.Rndm();
        double eff = 0;

        if (_pt > lowerLimit && _pt < upperLimit){
            int bin = histogram->FindBin(_pt);
            eff = histogram->GetBinContent(bin);
        }
        else if (_pt >= upperLimit){
            eff = histogram->GetBinContent(histogram->GetNbinsX());
        }

        if (randomNumber < eff) {
            double pt = smearedPt(_pt);
            double eta = particles.eta();
            double phi = particles.phi();
            double mass = 0.139; // pion mass
            float E = sqrt(pt*cos(phi)*pt*cos(phi)+pt*sin(phi)*pt*sin(phi)+pt*sinh(eta)*pt*sinh(eta)+mass*mass);

            PseudoJet particle = PseudoJet(pt*cos(phi),pt*sin(phi),pt*sinh(eta),E); //px,py,pz,E
            _detectorEvent.push_back(particle);
        }
    }

    return _detectorEvent;
  }
  
protected:
  vector<PseudoJet> _fullEvent;
  vector<PseudoJet> _truthEvent;
  vector<PseudoJet> _detectorEvent;

  TFile *file;
  TH1D *histogram;
  double lowerLimit, upperLimit;

  TRandom randomGenerator;

  double smearedPt (double _pt) { // Adapted from https://github.com/ezradlesser/pyjetty/blob/9e2eb21f0c1a74c1b576c1ec187b03fd447ddd87/pyjetty/alice_analysis/process/user/fastsim/eff_smear.py#L27
    double _pt_smeared;
    double sigma;

    if (_pt < 1){
        sigma =  _pt * (-0.035 * _pt + 0.04);
    }
    else if (_pt < 60) {
        sigma = _pt * (0.00085 * _pt + 0.00415);
    }
    else {
        sigma = _pt * (0.0015 * _pt - 0.035);
    }

    _pt_smeared = randomGenerator.Gaus(_pt, sigma);

    return _pt_smeared;
  }

};

#endif
