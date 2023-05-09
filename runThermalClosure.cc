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

#include "include/csSubFullEventIterative.hh"

#include "include/thermalEvent.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");
  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIThermalClosure.root").c_str(), "RECREATE");

  double R = cmdline.value<double>("-R",0.4);

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double ghostRapMax         = 1.0;
  double ghost_area          = 0.5;
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

    //---------------------------------------------------------------------------
    //   Embedding
    //---------------------------------------------------------------------------
    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    thermalEvent thrm(1000,0.7, -3.0, 3.0, 0.5);
    vector<PseudoJet> particlesBkg = thrm.createThermalEventAlice();

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //---------------------------------------------------------------------------
    //   jet clustering of TRUTH small R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_smaller(antikt_algorithm, R);

    fastjet::ClusterSequenceArea sigTruth_smaller(particlesSig, jet_def_smaller, area_def);
    jetCollection jetCollectionSig_Truth_smaller(sorted_by_pt(jet_selector(sigTruth_smaller.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   jet clustering of THERMAL small R
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubFullSmall( {0.} , {.25}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFullSmall.setInputParticles(particlesMerged);
    csSubFullSmall.setMaxEta(1.);
    csSubFullSmall.setBackgroundGrid();
    fastjet::ClusterSequenceArea sigThermal_smaller(csSubFullSmall.doSubtractionFullEvent(), jet_def_smaller, area_def);
    jetCollection csFullJets(sorted_by_pt(jet_selector(sigThermal_smaller.inclusive_jets(1.))));  

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets.getJet()) {
      if(jet.has_constituents())
        csFullJetsClean.push_back(jet);
    }

    jetCollection jetCollectionSig_Thermal_smaller(csFullJetsClean);

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(0.2);
    jmCSFull.setBaseJets(jetCollectionSig_Thermal_smaller);
    jmCSFull.setTagJets(jetCollectionSig_Truth_smaller);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionSig_Thermal_smaller);   

    //---------------------------------------------------------------------------
    //   jet clustering TRUTH of big R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_bigger(antikt_algorithm, R+0.05);

    fastjet::ClusterSequenceArea sigTruth_bigger(particlesSig, jet_def_bigger, area_def);
    jetCollection jetCollectionSig_Truth_bigger(sorted_by_pt(jet_selector(sigTruth_bigger.inclusive_jets(1.))));

    //match bigger to smaller jets
    jetMatcher jetMatch_bigger_Truth(0.2);
    jetMatch_bigger_Truth.setBaseJets(jetCollectionSig_Truth_bigger);
    jetMatch_bigger_Truth.setTagJets(jetCollectionSig_Truth_smaller);
    jetMatch_bigger_Truth.matchJets();
    jetMatch_bigger_Truth.reorderedToTag(jetCollectionSig_Truth_bigger);
    
    //---------------------------------------------------------------------------
    //   jet clustering of THERMAL bigger R
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubFullBig( {0.} , {.25}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFullBig.setInputParticles(particlesMerged);
    csSubFullBig.setMaxEta(1.);
    csSubFullBig.setBackgroundGrid();
    fastjet::ClusterSequenceArea sigThermal_bigger(csSubFullBig.doSubtractionFullEvent(), jet_def_bigger, area_def);
    jetCollection csFullJetsBig(sorted_by_pt(jet_selector(sigThermal_bigger.inclusive_jets(1.))));   

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsCleanBig;
    for(fastjet::PseudoJet jet : csFullJetsBig.getJet()) {
      if(jet.has_constituents())
        csFullJetsCleanBig.push_back(jet);
    }

    jetCollection jetCollectionSig_Thermal_bigger(csFullJetsCleanBig);

    //match CSFull jets to signal jets
    jetMatcher jetMatch_bigger_Thermal(0.2);
    jetMatch_bigger_Thermal.setBaseJets(jetCollectionSig_Thermal_bigger);
    jetMatch_bigger_Thermal.setTagJets(jetCollectionSig_Thermal_smaller);
    jetMatch_bigger_Thermal.matchJets();
    jetMatch_bigger_Thermal.reorderedToTag(jetCollectionSig_Thermal_bigger);  

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    trw.addCollection("sigJet_Truth_smaller",        jetCollectionSig_Truth_smaller);
    trw.addCollection("sigJet_Truth_bigger",        jetCollectionSig_Truth_bigger);
    trw.addCollection("sigJet_Thermal_smaller",        jetCollectionSig_Thermal_smaller);
    trw.addCollection("sigJet_Thermal_bigger",        jetCollectionSig_Thermal_bigger);

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
