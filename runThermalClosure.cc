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

  //double jetRapMax = 0.5;
  Selector jet_selector_truth = SelectorAbsEtaMax(1.0);
  Selector jet_selector_detector = SelectorAbsEtaMax(0.5);

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
    //   Defining data samples
    //---------------------------------------------------------------------------
    vector<PseudoJet> particlesMergedAll = mixer.particles();

    // Truth event
    fastjet::Selector hard_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesHard = hard_selector(particlesMergedAll);

    // Detector event
    fastjet::Selector reco_selector = SelectorVertexNumber(99);
    vector<PseudoJet> particlesReco = reco_selector(particlesMergedAll);

    // Pileup event
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesEmbedded = particlesBkg;
    particlesEmbedded.insert( particlesEmbedded.end(), particlesReco.begin(), particlesReco.end() );

    //---------------------------------------------------------------------------
    //   jet clustering of 3 samples
    //---------------------------------------------------------------------------
    /*
    JetDefinition jet_def(antikt_algorithm, R);
    fastjet::ClusterSequenceArea jets_Truth(particlesHard, jet_def, area_def);
    jetCollection jetCollection_Truth(sorted_by_pt(jet_selector(jets_Truth.inclusive_jets(10.))));
    trw.addCollection("jetCollection_Truth_",        jetCollection_Truth);

    fastjet::ClusterSequenceArea jets_Reco(particlesReco, jet_def, area_def);
    jetCollection jetCollection_Reco(sorted_by_pt(jet_selector(jets_Reco.inclusive_jets(10.))));
    trw.addCollection("jetCollection_Reco_",        jetCollection_Reco);

    // Randomly reject for PbPb kinematic efficiency

    fastjet::ClusterSequenceArea jets_MB(particlesBkg, jet_def, area_def);
    jetCollection jetCollection_MB(sorted_by_pt(jet_selector(jets_MB.inclusive_jets(10.))));
    trw.addCollection("jetCollection_MB_",        jetCollection_MB);

    fastjet::ClusterSequenceArea jets_Embedded(particlesEmbedded, jet_def, area_def);
    jetCollection jetCollection_Embedded(sorted_by_pt(jet_selector(jets_Embedded.inclusive_jets(10.))));
    trw.addCollection("jetCollection_Embedded_",        jetCollection_Embedded);
    */
    //---------------------------------------------------------------------------
    //   jet clustering of TRUTH small R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_smaller(antikt_algorithm, R);
    fastjet::ClusterSequenceArea sigTruth_smaller(particlesHard, jet_def_smaller, area_def);
    jetCollection jetCollectionSig_Truth_smaller(sorted_by_pt(jet_selector_truth(sigTruth_smaller.inclusive_jets(8.))));

    vector<double> TrackOver8GeV_Truth_smaller;      TrackOver8GeV_Truth_smaller.reserve(jetCollectionSig_Truth_smaller.getJet().size());
    double found;
    for(fastjet::PseudoJet jet : jetCollectionSig_Truth_smaller.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        TrackOver8GeV_Truth_smaller.push_back(found);
      }
    }
    jetCollectionSig_Truth_smaller.addVector("TrackOver8GeV_Truth_smaller", TrackOver8GeV_Truth_smaller);
    
    trw.addCollection("sigJet_Truth_smaller",        jetCollectionSig_Truth_smaller);

    //---------------------------------------------------------------------------
    //   jet clustering of THERMAL small R
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubFullSmall( {0.} , {.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFullSmall.setInputParticles(particlesEmbedded);
    csSubFullSmall.setMaxEta(1.);
    csSubFullSmall.setBackgroundGrid();
    fastjet::ClusterSequenceArea sigThermal_smaller(csSubFullSmall.doSubtractionFullEvent(), jet_def_smaller, area_def);
    jetCollection csFullJets(sorted_by_pt(jet_selector_detector(sigThermal_smaller.inclusive_jets(1.))));  
    

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets.getJet()) {
      if(jet.has_constituents())
        csFullJetsClean.push_back(jet);
    }

    jetCollection jetCollectionSig_Thermal_smaller(csFullJetsClean);

    trw.addCollection("sigJet_Thermal_smaller_beforeMatching",        jetCollectionSig_Thermal_smaller);

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(0.2);
    jmCSFull.setBaseJets(jetCollectionSig_Thermal_smaller);
    jmCSFull.setTagJets(jetCollectionSig_Truth_smaller);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionSig_Thermal_smaller);   

    vector<double> TrackOver8GeV_Thermal_smaller;      TrackOver8GeV_Thermal_smaller.reserve(jetCollectionSig_Thermal_smaller.getJet().size());
    double found;
    for(fastjet::PseudoJet jet : jetCollectionSig_Thermal_smaller.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        TrackOver8GeV_Thermal_smaller.push_back(found);
      }
    }
    jetCollectionSig_Thermal_smaller.addVector("TrackOver8GeV_Thermal_smaller", TrackOver8GeV_Thermal_smaller);

    trw.addCollection("sigJet_Thermal_smaller",        jetCollectionSig_Thermal_smaller);

    //---------------------------------------------------------------------------
    //   jet clustering TRUTH of big R
    //---------------------------------------------------------------------------
    JetDefinition jet_def_bigger(antikt_algorithm, R+0.05);

    fastjet::ClusterSequenceArea sigTruth_bigger(particlesHard, jet_def_bigger, area_def);
    jetCollection jetCollectionSig_Truth_bigger(sorted_by_pt(jet_selector_truth(sigTruth_bigger.inclusive_jets(8.))));

    trw.addCollection("sigJet_Truth_bigger_beforeMatching",        jetCollectionSig_Truth_bigger);

    //match bigger to smaller jets
    jetMatcher jetMatch_bigger_Truth(0.2);
    jetMatch_bigger_Truth.setBaseJets(jetCollectionSig_Truth_bigger);
    jetMatch_bigger_Truth.setTagJets(jetCollectionSig_Truth_smaller);
    jetMatch_bigger_Truth.matchJets();
    jetMatch_bigger_Truth.reorderedToTag(jetCollectionSig_Truth_bigger);

    vector<double> TrackOver8GeV_Truth_bigger;      TrackOver8GeV_Truth_bigger.reserve(jetCollectionSig_Truth_bigger.getJet().size());
    double found;
    for(fastjet::PseudoJet jet : jetCollectionSig_Truth_bigger.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        TrackOver8GeV_Truth_bigger.push_back(found);
      }
    }
    jetCollectionSig_Truth_bigger.addVector("TrackOver8GeV_Truth_bigger", TrackOver8GeV_Truth_bigger);

    trw.addCollection("sigJet_Truth_bigger",        jetCollectionSig_Truth_bigger);
    
    //---------------------------------------------------------------------------
    //   jet clustering of THERMAL bigger R
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubFullBig( {0.} , {.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFullBig.setInputParticles(particlesMerged);
    csSubFullBig.setMaxEta(1.);
    csSubFullBig.setBackgroundGrid();
    fastjet::ClusterSequenceArea sigThermal_bigger(csSubFullBig.doSubtractionFullEvent(), jet_def_bigger, area_def);
    jetCollection csFullJetsBig(sorted_by_pt(jet_selector_detector(sigThermal_bigger.inclusive_jets(8.))));   

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsCleanBig;
    for(fastjet::PseudoJet jet : csFullJetsBig.getJet()) {
      if(jet.has_constituents())
        csFullJetsCleanBig.push_back(jet);
    }

    jetCollection jetCollectionSig_Thermal_bigger(csFullJetsCleanBig);

    trw.addCollection("sigJet_Thermal_bigger_beforeMatching",        jetCollectionSig_Thermal_bigger);

    //match CSFull jets to signal jets
    jetMatcher jetMatch_bigger_Thermal(0.2);
    jetMatch_bigger_Thermal.setBaseJets(jetCollectionSig_Thermal_bigger);
    jetMatch_bigger_Thermal.setTagJets(jetCollectionSig_Thermal_smaller);
    jetMatch_bigger_Thermal.matchJets();
    jetMatch_bigger_Thermal.reorderedToTag(jetCollectionSig_Thermal_bigger);  

    vector<double> TrackOver8GeV_Thermal_bigger;      TrackOver8GeV_Thermal_bigger.reserve(jetCollectionSig_Thermal_bigger.getJet().size());
    double found;
    for(fastjet::PseudoJet jet : jetCollectionSig_Thermal_bigger.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        TrackOver8GeV_Thermal_bigger.push_back(found);
      }
    }
    jetCollectionSig_Thermal_bigger.addVector("TrackOver8GeV_Thermal_bigger", TrackOver8GeV_Thermal_bigger);

    trw.addCollection("sigJet_Thermal_bigger",        jetCollectionSig_Thermal_bigger);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------

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
