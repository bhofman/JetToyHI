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
#include "include/AliceFastSim.hh"

#include "include/csSubtractorFullEvent.hh"

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
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);
    
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
    auto acceptance = fastSim.AliceAcceptance();
    auto detector = fastSim.AliceDetector();
    //fastSim.test();

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea sig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(sig.inclusive_jets(25.))));

    fastjet::ClusterSequenceArea sig_acc(acceptance, jet_def, area_def);
    jetCollection jetCollectionSig_acc(sorted_by_pt(jet_selector(sig_acc.inclusive_jets(25.))));

    fastjet::ClusterSequenceArea sig_det(detector, jet_def, area_def);
    jetCollection jetCollectionSig_det(sorted_by_pt(jet_selector(sig_det.inclusive_jets(25.))));

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported
    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("sigJet_sim",        jetCollectionSig_acc);
    trw.addCollection("sigJet_det",        jetCollectionSig_det);

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
