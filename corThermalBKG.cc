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

using namespace std;
using namespace fastjet;

//./corThermalBKG -hard samples/Jewel/30/jewelNR_30_vac-8.hepmc -hardtype HepMC2 -pileup samples/Thermal/ThermalEventsMult7000PtAv1.20_1.pu14 -nev 100000 -output samples/AnalysisSamples/30/Jewel_30_cs-8.root -pt 30+10

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
  //fastjet::JetDefinition jet_def_ca(cambridge_algorithm, 999.); // Should no longer be needed

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

  Angularity mr(1.,0.,R);
  Angularity mr2(2.,0.,R);
  Angularity r2z(2.,1.,R);

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

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(user_pt)))); // Inclusive jets to take a jets with pt over (pt_min)

    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {2.,2.} , {.1,0.075}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    csSubFull.setBackgroundGrid();
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection csFullJets(sorted_by_pt(jet_selector(fullSig.inclusive_jets(1.)))); 

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

    //calculate some angularities
    //std::cout << "calc angularities groomed jets" << std::endl;
    vector<double> antiKT_width;      antiKT_width.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_pTD;        antiKT_pTD.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_mr;         antiKT_mr.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_mr2;        antiKT_mr2.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_r2z;        antiKT_r2z.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau1;       antiKT_tau1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau2;       antiKT_tau2.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau3;       antiKT_tau3.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau4;       antiKT_tau4.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau5;       antiKT_tau5.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau2tau1;   antiKT_tau2tau1.reserve(jetCollectionCS_Sig.getJet().size());
    vector<double> antiKT_tau3tau2;   antiKT_tau3tau2.reserve(jetCollectionCS_Sig.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      antiKT_width.push_back(width.result(jet));
      antiKT_pTD.push_back(pTD.result(jet));
      antiKT_mr.push_back(mr.result(jet));
      antiKT_mr2.push_back(mr2.result(jet));
      antiKT_r2z.push_back(r2z.result(jet));
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

    jetCollectionCS_Sig.addVector("antiKT_width", antiKT_width);
    jetCollectionCS_Sig.addVector("antiKT_ptd",   antiKT_pTD);
    jetCollectionCS_Sig.addVector("antiKT_mr",    antiKT_mr);
    jetCollectionCS_Sig.addVector("antiKT_mr2",   antiKT_mr2);
    jetCollectionCS_Sig.addVector("antiKT_r2z",   antiKT_r2z);
    jetCollectionCS_Sig.addVector("antiKT_tau1",  antiKT_tau1);
    jetCollectionCS_Sig.addVector("antiKT_tau2",  antiKT_tau2);
    jetCollectionCS_Sig.addVector("antiKT_tau3",  antiKT_tau3);
    jetCollectionCS_Sig.addVector("antiKT_tau4",  antiKT_tau4);
    jetCollectionCS_Sig.addVector("antiKT_tau5",  antiKT_tau5);
    jetCollectionCS_Sig.addVector("antiKT_tau2tau1", antiKT_tau2tau1);
    jetCollectionCS_Sig.addVector("antiKT_tau3tau2", antiKT_tau3tau2);

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
      ptPull.push_back((csFullJets.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(csFullJets.getJet()[i].pt()+jetCollectionSig.getJet()[i].pt()));
      mPull.push_back((csFullJets.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(csFullJets.getJet()[i].m()+jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPull",        ptPull);
    trw.addCollection("mPull",        mPull);
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);

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
    //std::cout << "calc angularities groomed jets" << std::endl;
    vector<double> SD_width;    SD_width.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_pTD;      SD_pTD.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_mr;         SD_mr.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_mr2;        SD_mr2.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_r2z;        SD_r2z.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau1;       SD_tau1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau2;       SD_tau2.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau3;       SD_tau3.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau4;       SD_tau4.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau5;       SD_tau5.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau2tau1;   SD_tau2tau1.reserve(jetCollectionCS_SD.getJet().size());
    vector<double> SD_tau3tau2;   SD_tau3tau2.reserve(jetCollectionCS_SD.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_SD.getJet()) {
      SD_width.push_back(width.result(jet));
      SD_pTD.push_back(pTD.result(jet));
      SD_mr.push_back(mr.result(jet));
      SD_mr2.push_back(mr2.result(jet));
      SD_r2z.push_back(r2z.result(jet));
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

    jetCollectionCS_SD.addVector("SD_width", SD_width);
    jetCollectionCS_SD.addVector("SD_ptd", SD_pTD);
    jetCollectionCS_SD.addVector("SD_mr", SD_mr);
    jetCollectionCS_SD.addVector("SD_mr2", SD_mr2);
    jetCollectionCS_SD.addVector("SD_r2z", SD_r2z);
    jetCollectionCS_SD.addVector("SD_tau1", SD_tau1);
    jetCollectionCS_SD.addVector("SD_tau2", SD_tau2);
    jetCollectionCS_SD.addVector("SD_tau3", SD_tau3);
    jetCollectionCS_SD.addVector("SD_tau4", SD_tau4);
    jetCollectionCS_SD.addVector("SD_tau5", SD_tau5);
    jetCollectionCS_SD.addVector("SD_tau2tau1", SD_tau2tau1);
    jetCollectionCS_SD.addVector("SD_tau3tau2", SD_tau3tau2);

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
    trw.addCollection("jetCollectionSig_",        jetCollectionSig);
    trw.addCollection("antiKT_",        jetCollectionCS_Sig);
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
