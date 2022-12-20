#include <iostream>
#include <chrono>
#include "TFile.h"
#include "TTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "include/ProgressBar.h"
#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"
#include "include/extraInfo.hh"
#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/Angularity.hh"
#include "include/dyGroomer.hh"
#include "include/csSubtractor.hh"
#include "include/csSubFullEventIterative.hh"
#include "include/jetCharge.hh"
#include "include/jetChargeDynamical.hh"
using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  int user_pt = cmdline.value<int>("-pt",10); 

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;     
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 2.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  Angularity Angularity_z1_theta1(1.0,1.,R);
  Angularity Angularity_z1_theta15( 1.5,1.,R);
  Angularity Angularity_z1_theta2(2.0,1.,R);
  Angularity Angularity_z1_theta3( 3.0,1.,R);

  Angularity Angularity_z2_theta1(1.0,2.,R);
  Angularity Angularity_z2_theta15( 1.5,2.,R);
  Angularity Angularity_z2_theta2(2.0,2.,R);
  Angularity Angularity_z2_theta3( 3.0,2.,R);  

  fastjet::contrib::OnePass_WTA_KT_Axes axes;
  fastjet::contrib::UnnormalizedMeasure unormbeta(1.0);
  fastjet::contrib::Nsubjettiness  nSub1_beta1(1, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub2_beta1(2, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub3_beta1(3, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub4_beta1(4, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub5_beta1(5, axes, unormbeta);
    
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number    
    iev++;
       
    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    
    //std::cout << "#merged: " << particlesMerged.size() << "  signal: " << particlesSig.size() << "  bkg: " << particlesBkg.size() << std::endl;

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(user_pt)))); // Inclusive jets to take a jets with pt over (pt_min)

    //calculate some angularities
    vector<double> z1_theta1;      z1_theta1.reserve(jetCollectionSig.getJet().size());
    vector<double> z1_theta15;     z1_theta15.reserve(jetCollectionSig.getJet().size());
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionSig.getJet().size());
    vector<double> z1_theta3;      z1_theta3.reserve(jetCollectionSig.getJet().size());

    vector<double> z2_theta1;      z2_theta1.reserve(jetCollectionSig.getJet().size());
    vector<double> z2_theta15;     z2_theta15.reserve(jetCollectionSig.getJet().size());
    vector<double> z2_theta2;      z2_theta2.reserve(jetCollectionSig.getJet().size());
    vector<double> z2_theta3;      z2_theta3.reserve(jetCollectionSig.getJet().size());  

    vector<double> antiKT_tau1;       antiKT_tau1.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau2;       antiKT_tau2.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau3;       antiKT_tau3.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau4;       antiKT_tau4.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau5;       antiKT_tau5.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau2tau1;   antiKT_tau2tau1.reserve(jetCollectionSig.getJet().size());
    vector<double> antiKT_tau3tau2;   antiKT_tau3tau2.reserve(jetCollectionSig.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      z1_theta15.push_back(Angularity_z1_theta15.result(jet));
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
      z1_theta3.push_back(Angularity_z1_theta3.result(jet));

      z2_theta1.push_back(Angularity_z2_theta1.result(jet));
      z2_theta15.push_back(Angularity_z2_theta15.result(jet));
      z2_theta2.push_back(Angularity_z2_theta2.result(jet));
      z2_theta3.push_back(Angularity_z2_theta3.result(jet));

      antiKT_tau1.push_back(nSub1_beta1(jet));
      antiKT_tau2.push_back(nSub2_beta1(jet));
      antiKT_tau3.push_back(nSub3_beta1(jet));
      antiKT_tau4.push_back(nSub4_beta1(jet));
      antiKT_tau5.push_back(nSub5_beta1(jet));

      if (nSub1_beta1(jet) != 0){
        antiKT_tau2tau1.push_back(nSub2_beta1(jet)/nSub1_beta1(jet));
      }
      if (nSub1_beta1(jet) == 0){
        //std::cout<<"Still zero tau1 "<<jet.constituents().size()<<std::endl;
        antiKT_tau2tau1.push_back(-999);
      }
      if (nSub2_beta1(jet) != 0){
        antiKT_tau3tau2.push_back(nSub3_beta1(jet)/nSub2_beta1(jet));
      }
      if (nSub2_beta1(jet) == 0){
        //std::cout<<"Still zero tau2 "<<jet.constituents().size()<<std::endl;
        antiKT_tau3tau2.push_back(-999);
      }
    }

    jetCollectionSig.addVector("z1_theta1", z1_theta1);
    jetCollectionSig.addVector("z1_theta15",z1_theta15);
    jetCollectionSig.addVector("z1_theta2", z1_theta2);
    jetCollectionSig.addVector("z1_theta3", z1_theta3);

    jetCollectionSig.addVector("z2_theta1", z2_theta1);
    jetCollectionSig.addVector("z2_theta15",z2_theta15);
    jetCollectionSig.addVector("z2_theta2", z2_theta2);
    jetCollectionSig.addVector("z2_theta3", z2_theta3);

    jetCollectionSig.addVector("tau1",  antiKT_tau1);
    jetCollectionSig.addVector("tau2",  antiKT_tau2);
    jetCollectionSig.addVector("tau3",  antiKT_tau3);
    jetCollectionSig.addVector("tau4",  antiKT_tau4);
    jetCollectionSig.addVector("tau5",  antiKT_tau5);
    jetCollectionSig.addVector("tau2tau1", antiKT_tau2tau1);
    jetCollectionSig.addVector("tau3tau2", antiKT_tau3tau2);
    
    //---------------------------------------------------------------------------
    //   Jet Charge
    //---------------------------------------------------------------------------

    vector<double> jetCharge;               jetCharge.reserve(jetCollectionSig.getJet().size());
    vector<double> jetChargeDynamical;      jetChargeDynamical.reserve(jetCollectionSig.getJet().size());

    JetCharge jetChargeFunction(0.5,-1); // kappa, ptmin
    JetChargeDynamical jetChargeDynamicalFunction(0.3,1.0,0.3,-1); // Xi, Kappa<, Kappa>, ptmin

    for(PseudoJet jet : jetCollectionSig.getJet()) {
      jetCharge.push_back(jetChargeFunction.result(jet));
      jetChargeDynamical.push_back(jetChargeDynamicalFunction.result(jet));
    }

    jetCollectionSig.addVector("jetCharge", jetCharge);
    jetCollectionSig.addVector("jetChargeDynamical", jetChargeDynamical);
    
    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the CS jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionCS_SD(sdgSigBeta00Z01.doGrooming(jetCollectionSig));

    jetCollectionCS_SD.addVector("SD_zg",    sdgSigBeta00Z01.getZgs());
    jetCollectionCS_SD.addVector("SD_ndrop", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionCS_SD.addVector("SD_dr12",  sdgSigBeta00Z01.getDR12());
    
    //calculate some angularities
    vector<double> SD_z1_theta1;      SD_z1_theta1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z1_theta15;     SD_z1_theta15.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z1_theta2;      SD_z1_theta2.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z1_theta3;      SD_z1_theta3.reserve(jetCollectionCS_SD.getJet().size());

    vector<double> SD_z2_theta1;      SD_z2_theta1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z2_theta15;     SD_z2_theta15.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z2_theta2;      SD_z2_theta2.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_z2_theta3;      SD_z2_theta3.reserve(jetCollectionCS_SD.getJet().size());

    vector<double> SD_tau1;       SD_tau1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau2;       SD_tau2.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau3;       SD_tau3.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau4;       SD_tau4.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau5;       SD_tau5.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau2tau1;   SD_tau2tau1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau3tau2;   SD_tau3tau2.reserve(jetCollectionCS_SD.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_SD.getJet()) {
      SD_z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      SD_z1_theta15.push_back(Angularity_z1_theta15.result(jet));
      SD_z1_theta2.push_back(Angularity_z1_theta2.result(jet));
      SD_z1_theta3.push_back(Angularity_z1_theta3.result(jet));

      SD_z2_theta1.push_back(Angularity_z2_theta1.result(jet));
      SD_z2_theta15.push_back(Angularity_z2_theta15.result(jet));
      SD_z2_theta2.push_back(Angularity_z2_theta2.result(jet));
      SD_z2_theta3.push_back(Angularity_z2_theta3.result(jet));

      SD_tau1.push_back(nSub1_beta1(jet));
      SD_tau2.push_back(nSub2_beta1(jet));
      SD_tau3.push_back(nSub3_beta1(jet));
      SD_tau4.push_back(nSub4_beta1(jet));
      SD_tau5.push_back(nSub5_beta1(jet));

      if (nSub1_beta1(jet) != 0){
        SD_tau2tau1.push_back(nSub2_beta1(jet)/nSub1_beta1(jet));
      }
      if (nSub1_beta1(jet) == 0){
        //std::cout<<"Still zero tau1 "<<jet.constituents().size()<<std::endl;
        SD_tau2tau1.push_back(-999);
      }
      if (nSub2_beta1(jet) != 0){
        SD_tau3tau2.push_back(nSub3_beta1(jet)/nSub2_beta1(jet));
      }
      if (nSub2_beta1(jet) == 0){
        //std::cout<<"Still zero tau2 "<<jet.constituents().size()<<std::endl;
        SD_tau3tau2.push_back(-999);
      }
    }

    jetCollectionCS_SD.addVector("SD_z1_theta1", SD_z1_theta1);
    jetCollectionCS_SD.addVector("SD_z1_theta15",SD_z1_theta15);
    jetCollectionCS_SD.addVector("SD_z1_theta2", SD_z1_theta2);
    jetCollectionCS_SD.addVector("SD_z1_theta3", SD_z1_theta3);

    jetCollectionCS_SD.addVector("SD_z2_theta1", SD_z2_theta1);
    jetCollectionCS_SD.addVector("SD_z2_theta15",SD_z2_theta15);
    jetCollectionCS_SD.addVector("SD_z2_theta2", SD_z2_theta2);
    jetCollectionCS_SD.addVector("SD_z2_theta3", SD_z2_theta3);

    jetCollectionCS_SD.addVector("SD_tau1", SD_tau1);
    jetCollectionCS_SD.addVector("SD_tau2", SD_tau2);
    jetCollectionCS_SD.addVector("SD_tau3", SD_tau3);
    jetCollectionCS_SD.addVector("SD_tau4", SD_tau4);
    jetCollectionCS_SD.addVector("SD_tau5", SD_tau5);
    jetCollectionCS_SD.addVector("SD_tau2tau1", SD_tau2tau1);
    jetCollectionCS_SD.addVector("SD_tau3tau2", SD_tau3tau2);

    //---------------------------------------------------------------------------
    //   SD Jet Charge
    //---------------------------------------------------------------------------

    vector<double> SDjetCharge;               SDjetCharge.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SDjetChargeDynamical;      SDjetChargeDynamical.reserve(jetCollectionCS_SD.getJet().size());

    for(PseudoJet jet : jetCollectionCS_SD.getJet()) {
      SDjetCharge.push_back(jetChargeFunction.result(jet));
      SDjetChargeDynamical.push_back(jetChargeDynamicalFunction.result(jet));
    }

    jetCollectionCS_SD.addVector("SD_jetCharge", SDjetCharge);
    jetCollectionCS_SD.addVector("SD_jetChargeDynamical", SDjetChargeDynamical);

    //---------------------------------------------------------------------------
    //   Dynamical grooming
    //---------------------------------------------------------------------------
    
    dyGroomer dygTDSig(2);
    jetCollection jetCollectionSigDYTD(dygTDSig.doGrooming(jetCollectionSig));
    trw.addCollection("kappa_TD",        dygTDSig.getKappas());
    trw.addCollection("zg_TD",        dygTDSig.getZgs());
    trw.addCollection("dR_TD",        dygTDSig.getDR12());
    
    dyGroomer dygKTDSig(1);
    jetCollection jetCollectionSigDYKTD(dygKTDSig.doGrooming(jetCollectionSig));
    trw.addCollection("kappa_KTD",        dygKTDSig.getKappas());
    trw.addCollection("zg_KTD",        dygKTDSig.getZgs());
    trw.addCollection("dR_KTD",        dygKTDSig.getDR12());
    
    dyGroomer dygzDSig(0.1);
    jetCollection jetCollectionSigDYzD(dygzDSig.doGrooming(jetCollectionSig));
    trw.addCollection("kappa_zD",        dygzDSig.getKappas());
    trw.addCollection("zg_zD",        dygzDSig.getZgs());
    trw.addCollection("dR_zD",        dygzDSig.getDR12());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    //trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("",        jetCollectionSig);
    trw.addCollection("SD_",      jetCollectionCS_SD);
    
  
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetNoBKG.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
