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
#include "include/csSubFullEventIterative.hh"
#include "include/dyGroomer.hh"
//#include "include/jetCharge.hh"

using namespace std;
using namespace fastjet;

// ./runAnalysis -hard samples/PythiaEventsTune14PtHat120.pu14 -pileup samples/ThermalEventsMult12000PtAv0.70.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

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
    
    //std::cout << "#merged: " << particlesMerged.size() << "  signal: " << particlesSig.size() << "  bkg: " << particlesBkg.size() << std::endl;

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(25.)))); // Inclusive jets to take a jets with pt over (pt_min)
    /*
    //calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    vector<double> tau2Sig;  tau2Sig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
      tau2Sig.push_back(nSub2_beta1(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);
    jetCollectionSig.addVector("tau2Sig", tau2Sig);
    */
    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {0.,0.} , {.2,.2}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    csSubFull.setBackgroundGrid();
    //This line crashes:
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection csFullJets(sorted_by_pt(jet_selector(fullSig.inclusive_jets(25.)))); 

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets.getJet()) {
      if(jet.has_constituents())
        csFullJetsClean.push_back(jet);
    }

    jetCollection jetCollectionCSFull(csFullJetsClean);

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSFull);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCSFull);

    /*
    //Background densities used by constituent subtraction
    std::vector<double> rhoFull;
    std::vector<double> rhomFull;
    rhoFull.push_back(csSubFull.getRho());  
    rhomFull.push_back(csSubFull.getRhoM()); 
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);
    */

    //trw.addCollection("csFull",        jetCollectionCSFull);
    
    /*
    std::vector<double> ptPull; ptPull.reserve(jetCollectionSig.getJet().size());
    std::vector<double> rapPull; rapPull.reserve(jetCollectionSig.getJet().size());
    std::vector<double> phiPull; phiPull.reserve(jetCollectionSig.getJet().size());
    std::vector<double> mPull; mPull.reserve(jetCollectionSig.getJet().size());
    
    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPull.push_back((jetCollectionCSFull.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionCSFull.getJet()[i].pt()+jetCollectionSig.getJet()[i].pt()));
      rapPull.push_back((jetCollectionCSFull.getJet()[i].rap()-jetCollectionSig.getJet()[i].rap())/(jetCollectionCSFull.getJet()[i].rap()+jetCollectionSig.getJet()[i].rap()));
      phiPull.push_back((jetCollectionCSFull.getJet()[i].phi()-jetCollectionSig.getJet()[i].phi())/(jetCollectionCSFull.getJet()[i].phi()+jetCollectionSig.getJet()[i].phi()));
      mPull.push_back((jetCollectionCSFull.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionCSFull.getJet()[i].m()+jetCollectionSig.getJet()[i].m()));
      //if (jetCollectionCSFull.getJet()[i].m() == 0){
      //  cout<<"I Found pt=0 with: "<<endl<<jetCollectionCSFull.getJet()[i]<<" signal: "<<jetCollectionSig.getJet()[i]<<endl;
      //  cout<<"Area: "<<jetCollectionSig.getJet()[i].pt()<<endl;
      //}
    }

    trw.addCollection("ptPull",        ptPull);
    trw.addCollection("rapPull",        rapPull);
    trw.addCollection("phiPull",        phiPull);
    trw.addCollection("mPull",        mPull);
    */


    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the CS jets
    //---------------------------------------------------------------------------
    // jetCollectionSig
    // jetCollectionCSFull

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("SD_zg",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("SD_ndrop", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("SD_dr12",  sdgSigBeta00Z01.getDR12());
    
    //calculate some angularities
    //std::cout << "calc angularities groomed jets" << std::endl;
    vector<double> SD_width;    SD_width.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_pTD;      SD_pTD.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_mr;         SD_mr.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_mr2;        SD_mr2.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_r2z;        SD_r2z.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau1;       SD_tau1.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau2;       SD_tau2.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau3;       SD_tau3.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau4;       SD_tau4.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau5;       SD_tau5.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau2tau1;   SD_tau2tau1.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    vector<double> SD_tau3tau2;   SD_tau3tau2.reserve(jetCollectionSigSDBeta00Z01.getJet().size());
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionSigSDBeta00Z01.getJet()) {
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
      if (nSub2_beta1(jet) != 0){
          SD_tau3tau2.push_back(nSub3_beta1(jet)/nSub2_beta1(jet));
      }
    }

    jetCollectionSigSDBeta00Z01.addVector("SD_width", SD_width);
    jetCollectionSigSDBeta00Z01.addVector("SD_ptd", SD_pTD);
    jetCollectionSigSDBeta00Z01.addVector("SD_mr", SD_mr);
    jetCollectionSigSDBeta00Z01.addVector("SD_mr2", SD_mr2);
    jetCollectionSigSDBeta00Z01.addVector("SD_r2z", SD_r2z);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau1", SD_tau1);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau2", SD_tau2);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau3", SD_tau3);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau4", SD_tau4);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau5", SD_tau5);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau2tau1", SD_tau2tau1);
    jetCollectionSigSDBeta00Z01.addVector("SD_tau3tau2", SD_tau3tau2);

    //---------------------------------------------------------------------------
    //   SD Jet Charge  // Need some work on adding PDG code from pu14
    //---------------------------------------------------------------------------
    /*
    jetCharge charge10(jetCollectionSigSDBeta00Z01,1.0);
    trw.addCollection("charge10",        charge10.getCharge());
    */

    //---------------------------------------------------------------------------
    //   Dynamical grooming
    //---------------------------------------------------------------------------
    //std::cout << "do dynamical grooming signal jets" << std::endl;
    dyGroomer dygTDSig(2);
    jetCollection jetCollectionSigDYTD(dygTDSig.doGrooming(csFullJetsClean));
    trw.addCollection("kappa_TD",        dygTDSig.getKappas());
    trw.addCollection("zg_TD",        dygTDSig.getZgs());
    trw.addCollection("dR_TD",        dygTDSig.getDR12());
    
    dyGroomer dygKTDSig(1);
    jetCollection jetCollectionSigDYKTD(dygKTDSig.doGrooming(csFullJetsClean));
    trw.addCollection("kappa_KTD",        dygKTDSig.getKappas());
    trw.addCollection("zg_KTD",        dygKTDSig.getZgs());
    trw.addCollection("dR_KTD",        dygKTDSig.getDR12());
    
    dyGroomer dygzDSig(0.1);
    jetCollection jetCollectionSigDYzD(dygzDSig.doGrooming(csFullJetsClean));
    trw.addCollection("kappa_zD",        dygzDSig.getKappas());
    trw.addCollection("zg_zD",        dygzDSig.getZgs());
    trw.addCollection("dR_zD",        dygzDSig.getDR12());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    //trw.addCollection("eventWeight",   eventWeight);
    //trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("SD_",      jetCollectionSigSDBeta00Z01);
  
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetBKG.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}