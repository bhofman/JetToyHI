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
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/AliceFastSim.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");
  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultSimpleJetAnalysis.root").c_str(), "RECREATE");

  double R = cmdline.value<double>("-R",0.);

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double ghostRapMax         = 1.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);

  double jetRapMax = 0.5;
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  AliceFastSim fastSim = AliceFastSim();

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number    
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesSig = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());

    //---------------------------------------------------------------------------
    //   fastsim
    //---------------------------------------------------------------------------
    fastSim.setInputEvent(particlesSig);
    vector<PseudoJet> truth = fastSim.AliceAcceptance();
    vector<PseudoJet> detector = fastSim.AliceDetector();

    //---------------------------------------------------------------------------
    //   jet clustering of small R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_smaller(antikt_algorithm, R);

    fastjet::ClusterSequenceArea sigTruth_smaller(truth, jet_def_smaller, area_def);
    jetCollection jetCollectionSig_Truth_smaller(sorted_by_pt(jet_selector(sigTruth_smaller.inclusive_jets(10.))));

    fastjet::ClusterSequenceArea sigDetector_smaller(detector, jet_def_smaller, area_def);
    jetCollection jetCollectionSig_Detector_smaller(sorted_by_pt(jet_selector(sigDetector_smaller.inclusive_jets(10.))));

    //match truth and detector jets
    jetMatcher jetMatch_smaller(R);
    jetMatch_smaller.setBaseJets(jetCollectionSig_Detector_smaller);
    jetMatch_smaller.setTagJets(jetCollectionSig_Truth_smaller);
    jetMatch_smaller.matchJets();
    jetMatch_smaller.reorderedToTag(jetCollectionSig_Detector_smaller);

    //---------------------------------------------------------------------------
    //   jet clustering of small R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_bigger(antikt_algorithm, R+0.05);

    fastjet::ClusterSequenceArea sigTruth_bigger(truth, jet_def_bigger, area_def);
    jetCollection jetCollectionSig_Truth_bigger(sorted_by_pt(jet_selector(sigTruth_bigger.inclusive_jets(10.))));

    fastjet::ClusterSequenceArea sigDetector_bigger(detector, jet_def_bigger, area_def);
    jetCollection jetCollectionSig_Detector_bigger(sorted_by_pt(jet_selector(sigDetector_bigger.inclusive_jets(10.))));

    //match bigger to smaller jets
    jetMatcher jetMatch_bigger_Truth(R+0.05);
    jetMatch_bigger_Truth.setBaseJets(jetCollectionSig_Truth_bigger);
    jetMatch_bigger_Truth.setTagJets(jetCollectionSig_Truth_smaller);
    jetMatch_bigger_Truth.matchJets();
    jetMatch_bigger_Truth.reorderedToTag(jetCollectionSig_Truth_bigger);

    jetMatcher jetMatch_bigger_Detector(R+0.05);
    jetMatch_bigger_Detector.setBaseJets(jetCollectionSig_Detector_bigger);
    jetMatch_bigger_Detector.setTagJets(jetCollectionSig_Detector_smaller);
    jetMatch_bigger_Detector.matchJets();
    jetMatch_bigger_Detector.reorderedToTag(jetCollectionSig_Detector_bigger);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    trw.addCollection("eventWeight",   eventWeight);

    trw.addCollection("sigJet_Truth_smaller",        jetCollectionSig_Truth_smaller);
    trw.addCollection("sigJet_Truth_bigger",        jetCollectionSig_Truth_bigger);
    trw.addCollection("sigJet_Detector_smaller",        jetCollectionSig_Detector_smaller);
    trw.addCollection("sigJet_Detector_bigger",        jetCollectionSig_Detector_bigger);

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
