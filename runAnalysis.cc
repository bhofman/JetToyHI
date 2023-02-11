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

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

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

    vector<PseudoJet> particlesSig = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());

    //---------------------------------------------------------------------------
    //   fastsim
    //---------------------------------------------------------------------------
    AliceFastSim fastSim(particlesSig);
    auto truth = fastSim.AliceAcceptance();
    auto detector = fastSim.AliceDetector();

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea sigTruth(truth, jet_def, area_def);
    jetCollection jetCollectionSig_Truth(sorted_by_pt(jet_selector(sigTruth.inclusive_jets(25.))));

    fastjet::ClusterSequenceArea sigDetector(detector, jet_def, area_def);
    jetCollection jetCollectionSig_Detector(sorted_by_pt(jet_selector(sigDetector.inclusive_jets(25.))));

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    trw.addCollection("eventWeight",   eventWeight);

    trw.addCollection("sigJet_Truth_",        jetCollectionSig_Truth);
    trw.addCollection("sigJet_Detector_",        jetCollectionSig_Detector);

    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultSimpleJetAnalysis.root").c_str(), "RECREATE");
  fout->cd();
  TTree *trOut = trw.getTree();
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
