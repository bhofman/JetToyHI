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

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);
  cout << "will run on " << nEvent << " events" << endl;

  TFile *fout = new TFile(cmdline.value<string>("-output", "PPMC.root").c_str(), "RECREATE");

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.2;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;     
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 1.5;
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  Angularity Angularity_z1_theta1(1.0,1.,R);
  Angularity Angularity_z1_theta2(2.0,1.,R);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle((nEvent == -1 ? 7 : -1));

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  AliceFastSim fastSim = AliceFastSim();

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

    vector<PseudoJet> particlesSig = mixer.particles();

    //---------------------------------------------------------------------------
    //   fastsim
    //---------------------------------------------------------------------------
    fastSim.setInputEvent(particlesSig);
    vector<PseudoJet> particlesTruth = fastSim.AliceAcceptance();
    vector<PseudoJet> particlesDetector = fastSim.AliceDetector();

    //---------------------------------------------------------------------------
    //   jet clustering of Truth jets
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csTruth(particlesTruth, jet_def, area_def);
    jetCollection jetCollectionTruth(sorted_by_pt(jet_selector(csTruth.inclusive_jets(20.))));

    // Angularities of truth matched jets
    vector<double> Truth_z1_theta1;      Truth_z1_theta1.reserve(jetCollectionTruth.getJet().size()); 
    vector<double> Truth_z1_theta2;      Truth_z1_theta2.reserve(jetCollectionTruth.getJet().size()); 
    for(PseudoJet jet : jetCollectionTruth.getJet()) {
      if (!jet.has_constituents()) { 
        Truth_z1_theta1.push_back(0);
        Truth_z1_theta2.push_back(0);
        } 
      else { 
        Truth_z1_theta1.push_back(Angularity_z1_theta1.result(jet));
        Truth_z1_theta2.push_back(Angularity_z1_theta2.result(jet));
        }
    }
    jetCollectionTruth.addVector("Truth_z1_theta1", Truth_z1_theta1);
    jetCollectionTruth.addVector("Truth_z1_theta2", Truth_z1_theta2);
    
    trw.addCollection("Truth_",      jetCollectionTruth);

    //              ## SD Groom jets ## 
    //SoftDrop grooming alice (zcut=0.2, beta=0)
    softDropGroomer Truth_SDGroomer(0.2, 0.0, R);
    jetCollection Truth_SD_(Truth_SDGroomer.doGrooming(jetCollectionTruth));
    Truth_SD_.addVector("Truth_SD_zg",    Truth_SDGroomer.getZgs());
    Truth_SD_.addVector("Truth_SD_ndrop", Truth_SDGroomer.getNDroppedSubjets());
    Truth_SD_.addVector("Truth_SD_dr12",  Truth_SDGroomer.getDR12());
    trw.addCollection("Truth_SD_",      Truth_SD_);

    //---------------------------------------------------------------------------
    //   Detector level jets before any matching
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csDetector(particlesDetector, jet_def, area_def);
    jetCollection jetCollectionDetector(sorted_by_pt(jet_selector(csDetector.inclusive_jets(15.)))); 

    //---------------------------------------------------------------------------
    //   Matching of detector jets to truth jets
    //---------------------------------------------------------------------------
    jetCollection jetCollectionDetectorMatched(jetCollectionDetector);  
    
    //match CSFull jets to signal jets
    jetMatcher jmDetector(0.6*R);
    jmDetector.setBaseJets(jetCollectionTruth);
    jmDetector.setTagJets(jetCollectionDetectorMatched);
    jmDetector.matchJets();
    jmDetector.reorderedToBase(jetCollectionDetectorMatched);

    // 100 GeV track finder
    vector<double> Detector_Matched_TrackOver100GeV; double found;
    Detector_Matched_TrackOver100GeV.reserve(jetCollectionDetectorMatched.getJet().size());
    for(fastjet::PseudoJet jet : jetCollectionDetectorMatched.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        Detector_Matched_TrackOver100GeV.push_back(found);
      }
      else { Detector_Matched_TrackOver100GeV.push_back(0); }
    } 
    jetCollectionDetectorMatched.addVector("Detector_Matched_TrackOver100GeV", Detector_Matched_TrackOver100GeV); 
        
    // Angularities
    vector<double> Detector_Matched_z1_theta1;      Detector_Matched_z1_theta1.reserve(jetCollectionDetectorMatched.getJet().size()); 
    vector<double> Detector_Matched_z1_theta2;      Detector_Matched_z1_theta2.reserve(jetCollectionDetectorMatched.getJet().size()); 
    for(PseudoJet jet : jetCollectionDetectorMatched.getJet()) {
      if (!jet.has_constituents()) { 
        Detector_Matched_z1_theta1.push_back(0);
        Detector_Matched_z1_theta2.push_back(0);
        } 
      else { 
        Detector_Matched_z1_theta1.push_back(Angularity_z1_theta1.result(jet));
        Detector_Matched_z1_theta2.push_back(Angularity_z1_theta2.result(jet));
        }
    }
    jetCollectionDetectorMatched.addVector("Detector_Matched_z1_theta1", Detector_Matched_z1_theta1);
    jetCollectionDetectorMatched.addVector("Detector_Matched_z1_theta2", Detector_Matched_z1_theta2);

    trw.addCollection("Detector_Matched_",     jetCollectionDetectorMatched);

    //              ## SD Groom jets ## 
    //SoftDrop grooming alice (zcut=0.2, beta=0)
    softDropGroomer Detector_Matched_SDGroomer(0.2, 0.0, R);
    jetCollection Detector_Matched_SD_(Detector_Matched_SDGroomer.doGrooming(jetCollectionDetectorMatched));
    Detector_Matched_SD_.addVector("Detector_Matched_SD_zg",    Detector_Matched_SDGroomer.getZgs());
    Detector_Matched_SD_.addVector("Detector_Matched_SD_ndrop", Detector_Matched_SDGroomer.getNDroppedSubjets());
    Detector_Matched_SD_.addVector("Detector_Matched_SD_dr12",  Detector_Matched_SDGroomer.getDR12());
    trw.addCollection("Detector_Matched_SD_",      Detector_Matched_SD_);
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported
    trw.addCollection("eventWeight",   eventWeight);
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
