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

  double jetRapMax = 0.5;
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
    
    trw.addCollection("eventWeight",   eventWeight);
    trw.fillTree();

    // run as --hard DETECTOR
    fastjet::Selector detector_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesDetector = detector_selector(particlesMergedAll);

    // run as --pileup TRUTH
    fastjet::Selector truth_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesTruth = truth_selector(particlesMergedAll);
    
    //---------------------------------------------------------------------------
    //   jet clustering of Detector jets
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csDetector(particlesDetector, jet_def, area_def);
    jetCollection jetCollectionDetector(sorted_by_pt(jet_selector(csDetector.inclusive_jets(user_pt)))); // Inclusive jets to take a jets with pt over (pt_min)

    //calculate some angularities
    vector<double> z1_theta1;      z1_theta1.reserve(jetCollectionDetector.getJet().size());
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionDetector.getJet().size());
    vector<double> z2_theta1;      z2_theta1.reserve(jetCollectionDetector.getJet().size());
    vector<double> z2_theta2;      z2_theta2.reserve(jetCollectionDetector.getJet().size());  

    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionDetector.getJet()) {
      z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
      z2_theta1.push_back(Angularity_z2_theta1.result(jet));
      z2_theta2.push_back(Angularity_z2_theta2.result(jet));
    }

    jetCollectionDetector.addVector("Det_z1_theta1", z1_theta1);
    jetCollectionDetector.addVector("Det_z1_theta2", z1_theta2);
    jetCollectionDetector.addVector("Det_z2_theta1", z2_theta1);
    jetCollectionDetector.addVector("Det_z2_theta2", z2_theta2);

    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the Detector jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0) // zcut=0.2
    softDropGroomer sdgSigBeta00Z01_Detector(0.2, 0.0, R);
    jetCollection jetCollectionDetector_SD(sdgSigBeta00Z01_Detector.doGrooming(jetCollectionDetector));

    jetCollectionDetector_SD.addVector("Det_SD_zg",    sdgSigBeta00Z01_Detector.getZgs());
    jetCollectionDetector_SD.addVector("Det_SD_ndrop", sdgSigBeta00Z01_Detector.getNDroppedSubjets());
    jetCollectionDetector_SD.addVector("Det_SD_dr12",  sdgSigBeta00Z01_Detector.getDR12());

    //---------------------------------------------------------------------------
    //   jet clustering of Truth jets
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csTruth(particlesTruth, jet_def, area_def);
    jetCollection jetCollectionTruth(sorted_by_pt(jet_selector(csTruth.inclusive_jets(1.)))); // Inclusive jets to take a jets with pt over (pt_min)

    trw.addCollection("Truth_",      jetCollectionTruth);

    //match truth jets to detector jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionTruth);
    jmCSFull.setTagJets(jetCollectionDetector);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionTruth);

    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> MatchedEvent;
    for(fastjet::PseudoJet jet : jetCollectionTruth.getJet()) {
      if(jet.has_constituents()){
        MatchedEvent.push_back(jet);
      }
    }
    jetCollection jetCollectionTruthMatched(MatchedEvent);    

    //calculate some angularities
    vector<double> z1_theta1_truth;      z1_theta1_truth.reserve(jetCollectionTruthMatched.getJet().size());
    vector<double> z1_theta2_truth;      z1_theta2_truth.reserve(jetCollectionTruthMatched.getJet().size());
    vector<double> z2_theta1_truth;      z2_theta1_truth.reserve(jetCollectionTruthMatched.getJet().size());
    vector<double> z2_theta2_truth;      z2_theta2_truth.reserve(jetCollectionTruthMatched.getJet().size());  

    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionTruthMatched.getJet()) {
      z1_theta1_truth.push_back(Angularity_z1_theta1.result(jet));
      z1_theta2_truth.push_back(Angularity_z1_theta2.result(jet));
      z2_theta1_truth.push_back(Angularity_z2_theta1.result(jet));
      z2_theta2_truth.push_back(Angularity_z2_theta2.result(jet));
    }

    jetCollectionTruthMatched.addVector("Truth_z1_theta1", z1_theta1_truth);
    jetCollectionTruthMatched.addVector("Truth_z1_theta2", z1_theta2_truth);
    jetCollectionTruthMatched.addVector("Truth_z2_theta1", z2_theta1_truth);
    jetCollectionTruthMatched.addVector("Truth_z2_theta2", z2_theta2_truth);

    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the Truth jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0) // zcut=0.2
    softDropGroomer sdgSigBeta00Z01_Truth(0.2, 0.0, R);
    jetCollection jetCollectionTruthMatched_SD(sdgSigBeta00Z01_Truth.doGrooming(jetCollectionTruthMatched));

    jetCollectionTruthMatched_SD.addVector("Truth_Matched_SD_zg",    sdgSigBeta00Z01_Truth.getZgs());
    jetCollectionTruthMatched_SD.addVector("Truth_Matched_SD_ndrop", sdgSigBeta00Z01_Truth.getNDroppedSubjets());
    jetCollectionTruthMatched_SD.addVector("Truth_Matched_SD_dr12",  sdgSigBeta00Z01_Truth.getDR12());
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("Det_",        jetCollectionDetector);
    trw.addCollection("Det_SD_",      jetCollectionDetector_SD);
    trw.addCollection("Truth_Matched_",        jetCollectionTruthMatched);
    trw.addCollection("Truth_Matched_SD_",      jetCollectionTruthMatched_SD);
    
  
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
