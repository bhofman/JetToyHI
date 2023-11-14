#include <iostream>
#include <chrono>
#include "TFile.h"
#include "TTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
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

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);
  cout << "will run on " << nEvent << " events" << endl;

  TFile *fout = new TFile(cmdline.value<string>("-output", "PPMC.root").c_str(), "RECREATE");
  int user_pt = cmdline.value<int>("-pt",1); 

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

  double jetRapMax = 1.5;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);
  //Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  Angularity Angularity_z1_theta1(1.0,1.,R);
  Angularity Angularity_z1_theta2(2.0,1.,R);
  Angularity Angularity_z2_theta1(1.0,2.,R);
  Angularity Angularity_z2_theta2(2.0,2.,R);
  
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle((nEvent == -1 ? 7 : -1));

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && ( iev < nEvent || nEvent == -1 ) )
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
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.)))); // Inclusive jets to take a jets with pt over (pt_min)

        
    //calculate some angularities
    vector<double> z1_theta1;      z1_theta1.reserve(jetCollectionSig.getJet().size());
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionSig.getJet().size());
    vector<double> z2_theta1;      z2_theta1.reserve(jetCollectionSig.getJet().size());
    vector<double> z2_theta2;      z2_theta2.reserve(jetCollectionSig.getJet().size());

    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
      z2_theta1.push_back(Angularity_z2_theta1.result(jet));
      z2_theta2.push_back(Angularity_z2_theta2.result(jet));
    }

    jetCollectionSig.addVector("z1_theta1", z1_theta1);
    jetCollectionSig.addVector("z1_theta2", z1_theta2);
    jetCollectionSig.addVector("z2_theta1", z2_theta1);
    jetCollectionSig.addVector("z2_theta2", z2_theta2);

    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the signal jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollection_SD(sdgSigBeta00Z01.doGrooming(jetCollectionSig));

    jetCollection_SD.addVector("SD_zg",    sdgSigBeta00Z01.getZgs());
    jetCollection_SD.addVector("SD_ndrop", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollection_SD.addVector("SD_dr12",  sdgSigBeta00Z01.getDR12());

    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {2.,2.} , {.2,0.05}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection jetCollectionCS_Sig(sorted_by_pt(jet_selector(fullSig.inclusive_jets(10.)))); 
    
    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCS_Sig);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCS_Sig);

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> MatchedEvent;
    for(fastjet::PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      if(jet.has_constituents()){
        MatchedEvent.push_back(jet);
      }
    }
    jetCollection jetCollectionCS_Sig_Matched(MatchedEvent);
    
    //calculate some angularities
    vector<double> z1_theta1_CS;      z1_theta1_CS.reserve(jetCollectionCS_Sig_Matched.getJet().size());
    vector<double> z1_theta2_CS;      z1_theta2_CS.reserve(jetCollectionCS_Sig_Matched.getJet().size());
    vector<double> z2_theta1_CS;      z2_theta1_CS.reserve(jetCollectionCS_Sig_Matched.getJet().size());
    vector<double> z2_theta2_CS;      z2_theta2_CS.reserve(jetCollectionCS_Sig_Matched.getJet().size());

    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_Sig_Matched.getJet()) {
      z1_theta1_CS.push_back(Angularity_z1_theta1.result(jet));
      z1_theta2_CS.push_back(Angularity_z1_theta2.result(jet));
      z2_theta1_CS.push_back(Angularity_z2_theta1.result(jet));
      z2_theta2_CS.push_back(Angularity_z2_theta2.result(jet));
    }

    jetCollectionCS_Sig_Matched.addVector("CS_z1_theta1", z1_theta1_CS);
    jetCollectionCS_Sig_Matched.addVector("CS_z1_theta2", z1_theta2_CS);
    jetCollectionCS_Sig_Matched.addVector("CS_z2_theta1", z2_theta1_CS);
    jetCollectionCS_Sig_Matched.addVector("CS_z2_theta2", z2_theta2_CS);

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
      ptPull.push_back((jetCollectionCS_Sig_Matched.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionSig.getJet()[i].pt()));
      mPull.push_back((jetCollectionCS_Sig_Matched.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPull",        ptPull);
    trw.addCollection("mPull",        mPull);
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);
    
    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the CS jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01_CS(0.1, 0.0, R);
    jetCollection jetCollectionCS_SD(sdgSigBeta00Z01_CS.doGrooming(jetCollectionCS_Sig));

    jetCollectionCS_SD.addVector("SD_zg",    sdgSigBeta00Z01_CS.getZgs());
    jetCollectionCS_SD.addVector("SD_ndrop", sdgSigBeta00Z01_CS.getNDroppedSubjets());
    jetCollectionCS_SD.addVector("SD_dr12",  sdgSigBeta00Z01_CS.getDR12());
   
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("",     jetCollectionCS_Sig);
    trw.addCollection("SD_",      jetCollectionCS_SD);
    
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  fout->cd();
  TTree *trOut = trw.getTree();
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
