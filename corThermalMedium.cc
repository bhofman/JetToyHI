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
  double user_r = cmdline.value<double>("-r",0.2); 
  double user_alpha = cmdline.value<double>("-alpha",0.0); 

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

    // Jewel sub:
    fastjet::Selector dummy_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> particlesDummy = dummy_selector(particlesMergedAll);

    for(int i = 0; i < (int)particlesDummy.size(); i++)
    {
       if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000)
       {
          particlesDummy.erase(particlesDummy.begin() + i);
          i = i - 1;
       }
    }

    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);  // distance in eta-phi plane
    subtractor.set_max_distance(
        user_r);  // free parameter for the maximal allowed distance between particle i and ghost k
    subtractor.set_alpha(
        user_alpha);  // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
    subtractor.set_scale_fourmomentum();
    subtractor.set_remove_all_zero_pt_particles(true);

    std::vector<fastjet::PseudoJet> subtracted_particles = subtractor.do_subtraction(particlesSig, particlesDummy);

    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), subtracted_particles.begin(), subtracted_particles.end() );

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------
    
    fastjet::ClusterSequenceArea csSig(subtracted_particles, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(5.)))); // Inclusive jets to take a jets with pt over (pt_min)

    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {2.,2.} , {.2,0.05}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection jetCollectionCS_Sig(sorted_by_pt(jet_selector(fullSig.inclusive_jets(user_pt)))); 
    /*
    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(csFullJets);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(csFullJets);

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets.getJet()) {
      if(jet.has_constituents()){
        csFullJetsClean.push_back(jet);
      }
    }
    jetCollection jetCollectionCS_Sig(csFullJetsClean);
    */
    //calculate some angularities
    vector<double> z1_theta1;      z1_theta1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z1_theta15;     z1_theta15.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z1_theta3;      z1_theta3.reserve(jetCollectionCS_Sig.getJet().size());

    vector<double> z2_theta1;      z2_theta1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z2_theta15;     z2_theta15.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z2_theta2;      z2_theta2.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> z2_theta3;      z2_theta3.reserve(jetCollectionCS_Sig.getJet().size());  

    vector<double> tau1;       tau1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau2;       tau2.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau3;       tau3.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau4;       tau4.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau5;       tau5.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau2tau1;   tau2tau1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> tau3tau2;   tau3tau2.reserve(jetCollectionCS_Sig.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      z1_theta15.push_back(Angularity_z1_theta15.result(jet));
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
      z1_theta3.push_back(Angularity_z1_theta3.result(jet));

      z2_theta1.push_back(Angularity_z2_theta1.result(jet));
      z2_theta15.push_back(Angularity_z2_theta15.result(jet));
      z2_theta2.push_back(Angularity_z2_theta2.result(jet));
      z2_theta3.push_back(Angularity_z2_theta3.result(jet));

      tau1.push_back(nSub1_beta1(jet));
      tau2.push_back(nSub2_beta1(jet));
      tau3.push_back(nSub3_beta1(jet));
      tau4.push_back(nSub4_beta1(jet));
      tau5.push_back(nSub5_beta1(jet));

      if (nSub1_beta1(jet) != 0){
        tau2tau1.push_back(nSub2_beta1(jet)/nSub1_beta1(jet));
      }
      if (nSub1_beta1(jet) == 0){
        //std::cout<<"Still zero tau1 "<<jet.constituents().size()<<std::endl;
        tau2tau1.push_back(-999);
      }
      if (nSub2_beta1(jet) != 0){
        tau3tau2.push_back(nSub3_beta1(jet)/nSub2_beta1(jet));
      }
      if (nSub2_beta1(jet) == 0){
        //std::cout<<"Still zero tau2 "<<jet.constituents().size()<<std::endl;
        tau3tau2.push_back(-999);
      }
    }

    jetCollectionCS_Sig.addVector("z1_theta1", z1_theta1);
    jetCollectionCS_Sig.addVector("z1_theta15",z1_theta15);
    jetCollectionCS_Sig.addVector("z1_theta2", z1_theta2);
    jetCollectionCS_Sig.addVector("z1_theta3", z1_theta3);

    jetCollectionCS_Sig.addVector("z2_theta1", z2_theta1);
    jetCollectionCS_Sig.addVector("z2_theta15",z2_theta15);
    jetCollectionCS_Sig.addVector("z2_theta2", z2_theta2);
    jetCollectionCS_Sig.addVector("z2_theta3", z2_theta3);

    jetCollectionCS_Sig.addVector("tau1",  tau1);
    jetCollectionCS_Sig.addVector("tau2",  tau2);
    jetCollectionCS_Sig.addVector("tau3",  tau3);
    jetCollectionCS_Sig.addVector("tau4",  tau4);
    jetCollectionCS_Sig.addVector("tau5",  tau5);
    jetCollectionCS_Sig.addVector("tau2tau1", tau2tau1);
    jetCollectionCS_Sig.addVector("tau3tau2", tau3tau2);

    //---------------------------------------------------------------------------
    //   Jet Charge
    //---------------------------------------------------------------------------

    vector<double> jetCharge;               jetCharge.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> jetChargeDynamical;      jetChargeDynamical.reserve(jetCollectionCS_Sig.getJet().size());

    JetCharge jetChargeFunction(0.5,-1); // kappa, ptmin
    JetChargeDynamical jetChargeDynamicalFunction(0.3,1.0,0.3,-1); // Xi, Kappa<, Kappa>, ptmin

    for(PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      jetCharge.push_back(jetChargeFunction.result(jet));
      jetChargeDynamical.push_back(jetChargeDynamicalFunction.result(jet));
    }

    jetCollectionCS_Sig.addVector("jetCharge", jetCharge);
    jetCollectionCS_Sig.addVector("jetChargeDynamical", jetChargeDynamical);
    /*
    //---------------------------------------------------------------------------
    //   CS test statistics
    //---------------------------------------------------------------------------
    //Background densities used by constituent subtraction
    std::vector<double> rhoFull;
    std::vector<double> rhomFull;
    rhoFull.push_back(csSubFull.getRho());  
    rhomFull.push_back(csSubFull.getRhoM()); 
    
    std::vector<double> ptPull; ptPull.reserve(jetCollectionSig.getJet().size());
    std::vector<double> mPull; mPull.reserve(jetCollectionSig.getJet().size());
    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPull.push_back((csFullJets.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionSig.getJet()[i].pt()));
      mPull.push_back((csFullJets.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPull",        ptPull);
    trw.addCollection("mPull",        mPull);
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);
    */
    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the CS jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionCS_SD(sdgSigBeta00Z01.doGrooming(jetCollectionCS_Sig));

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
        SD_tau2tau1.push_back(-999);
      }
      if (nSub2_beta1(jet) != 0){
        SD_tau3tau2.push_back(nSub3_beta1(jet)/nSub2_beta1(jet));
      }
      if (nSub2_beta1(jet) == 0){
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
    jetCollection jetCollectionSigDYTD(dygTDSig.doGrooming(jetCollectionCS_Sig));
    trw.addCollection("kappa_TD",        dygTDSig.getKappas());
    trw.addCollection("zg_TD",        dygTDSig.getZgs());
    trw.addCollection("dR_TD",        dygTDSig.getDR12());
    
    dyGroomer dygKTDSig(1);
    jetCollection jetCollectionSigDYKTD(dygKTDSig.doGrooming(jetCollectionCS_Sig));
    trw.addCollection("kappa_KTD",        dygKTDSig.getKappas());
    trw.addCollection("zg_KTD",        dygKTDSig.getZgs());
    trw.addCollection("dR_KTD",        dygKTDSig.getDR12());
    
    dyGroomer dygzDSig(0.1);
    jetCollection jetCollectionSigDYzD(dygzDSig.doGrooming(jetCollectionCS_Sig));
    trw.addCollection("kappa_zD",        dygzDSig.getKappas());
    trw.addCollection("zg_zD",        dygzDSig.getZgs());
    trw.addCollection("dR_zD",        dygzDSig.getDR12());
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    //trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("",     jetCollectionCS_Sig);
    trw.addCollection("SD_",      jetCollectionCS_SD);
    
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetThermalBKG.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
