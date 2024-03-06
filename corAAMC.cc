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
#include "include/Angularity.hh"
#include "include/csSubFullEventIterative.hh"

#include <TRandom.h>

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
  double ghostRapMax         = 1.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;     
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 1.2;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  Angularity Angularity_z1_theta2(2.0,1.,R);
  
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

    fastjet::Selector hard_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesTruth = hard_selector(particlesMergedAll);
    
    fastjet::Selector reco_selector = SelectorVertexNumber(99);
    vector<PseudoJet> particlesReco = reco_selector(particlesMergedAll);

    // Randomly reject 2% of tracks
    TRandom randomGenerator;
    vector<PseudoJet> particlesReducedTracking;
    for(fastjet::PseudoJet particle : particlesReco) {
        double randomNumber = randomGenerator.Rndm();
        double eff = 0.98;
        if (randomNumber < eff) {
            particlesReducedTracking.push_back(particle);
        }
    }

    fastjet::Selector pileup_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesPileup = pileup_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesReducedTracking;
    particlesMerged.insert( particlesMerged.end(), particlesPileup.begin(), particlesPileup.end() );
    
    //---------------------------------------------------------------------------
    //   Truth level jets before any matching
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csTruth(particlesTruth, jet_def, area_def);
    jetCollection jetCollectionTruth(sorted_by_pt(jet_selector(csTruth.inclusive_jets(15.)))); 

    // Angularities of truth matched jets
    vector<double> z1_theta2_truth;      z1_theta2_truth.reserve(jetCollectionTruth.getJet().size()); 
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionTruth.getJet()) {
      if (!jet.has_constituents()) {
        z1_theta2_truth.push_back(0);
      } else {
        z1_theta2_truth.push_back(Angularity_z1_theta2.result(jet));
      }
    }
    jetCollectionTruth.addVector("Truth_z1_theta2", z1_theta2_truth);
    trw.addCollection("Truth_",     jetCollectionTruth);
    //---------------------------------------------------------------------------
    //   Detector level jets before any matching
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csDetector(particlesReducedTracking, jet_def, area_def);
    jetCollection jetCollectionDetector(sorted_by_pt(jet_selector(csDetector.inclusive_jets(15.)))); 

    // 100 GeV track finder
    double found;
    vector<double> Detector_TrackOver100GeV;      Detector_TrackOver100GeV.reserve(jetCollectionDetector.getJet().size()); 
    for(fastjet::PseudoJet jet : jetCollectionDetector.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        Detector_TrackOver100GeV.push_back(found);
      }
    } 
    jetCollectionDetector.addVector("Detector_TrackOver100GeV", Detector_TrackOver100GeV); 
 
    //calculate some angularities
    vector<double> z1_theta2_detector;      z1_theta2_detector.reserve(jetCollectionDetector.getJet().size()); 
    for(PseudoJet jet : jetCollectionDetector.getJet()) {
      if (!jet.has_constituents()) {
        z1_theta2_detector.push_back(0);
      } else {
        z1_theta2_detector.push_back(Angularity_z1_theta2.result(jet));
      }
    }
    jetCollectionDetector.addVector("Detector_z1_theta2", z1_theta2_detector);

    trw.addCollection("Detector_",     jetCollectionDetector);
    
    //---------------------------------------------------------------------------
    //   Embedded level jets before any matching
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubEmbedded( {0.0} , {0.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubEmbedded.setInputParticles(particlesMerged);
    csSubEmbedded.setMaxEta(1.0);
    fastjet::ClusterSequenceArea csEmbedded(csSubEmbedded.Subtract(), jet_def, area_def);
    jetCollection jetCollectionEmbedded(sorted_by_pt(jet_selector(csEmbedded.inclusive_jets(15.)))); 

    // 100 GeV track finder
    vector<double> Embedded_TrackOver100GeV;      Embedded_TrackOver100GeV.reserve(jetCollectionEmbedded.getJet().size());
    for(fastjet::PseudoJet jet : jetCollectionEmbedded.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        Embedded_TrackOver100GeV.push_back(found);
      }
    } 
    jetCollectionEmbedded.addVector("Embedded_TrackOver100GeV", Embedded_TrackOver100GeV); 
  
    //calculate some angularities
    vector<double> Embedded_z1_theta2;      Embedded_z1_theta2.reserve(jetCollectionEmbedded.getJet().size());
    for(PseudoJet jet : jetCollectionEmbedded.getJet()) {
      Embedded_z1_theta2.push_back(Angularity_z1_theta2.result(jet));
    }
    jetCollectionEmbedded.addVector("Embedded_z1_theta2", Embedded_z1_theta2);

    trw.addCollection("Embedded_",     jetCollectionEmbedded);

    //---------------------------------------------------------------------------
    //   Matching of detector jets to truth jets
    //---------------------------------------------------------------------------
    //match CSFull jets to signal jets
    jetMatcher jmTruth(0.6*R);
    jmTruth.setBaseJets(jetCollectionTruth);
    jmTruth.setTagJets(jetCollectionDetector);
    jmTruth.matchJets();
    jmTruth.reorderedToBase(jetCollectionDetector);
    
    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> MatchedEventDetector;
    for(fastjet::PseudoJet jet : jetCollectionDetector.getJet()) {
      if(jet.has_constituents()){
        MatchedEventDetector.push_back(jet);
      }
    }
    jetCollection jetCollectionDetectorMatched(MatchedEventDetector);

    // 100 GeV track finder
    vector<double> Detector_Matched_TrackOver100GeV;
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
    } 
    jetCollectionDetectorMatched.addVector("Detector_Matched_TrackOver100GeV", Detector_Matched_TrackOver100GeV); 
        
    vector<double> z1_theta2_detector_matched;      z1_theta2_detector_matched.reserve(jetCollectionDetectorMatched.getJet().size()); 
    for(PseudoJet jet : jetCollectionDetectorMatched.getJet()) {
        z1_theta2_detector_matched.push_back(Angularity_z1_theta2.result(jet));
    }
    jetCollectionDetectorMatched.addVector("Detector_Matched_z1_theta2", z1_theta2_detector_matched);
    
    trw.addCollection("Detector_Matched_",     jetCollectionDetectorMatched);

    //---------------------------------------------------------------------------
    //   Matching of embedded jets to detector jets
    //---------------------------------------------------------------------------
    //match CSFull jets to signal jets
    jetMatcher jmEmbedded(0.6*R);
    jmEmbedded.setBaseJets(jetCollectionDetector);
    jmEmbedded.setTagJets(jetCollectionEmbedded);
    jmEmbedded.matchJets();
    jmEmbedded.reorderedToBase(jetCollectionEmbedded);
    
    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> MatchedEventEmbedded;
    for(fastjet::PseudoJet jet : jetCollectionEmbedded.getJet()) {
      if(jet.has_constituents()){
        MatchedEventEmbedded.push_back(jet);
      }
    }
    jetCollection jetCollectionEmbeddedMatched(MatchedEventEmbedded);

    // 100 GeV track finder
    vector<double> Embedded_Matched_TrackOver100GeV;
    Embedded_Matched_TrackOver100GeV.reserve(jetCollectionEmbeddedMatched.getJet().size());
    for(fastjet::PseudoJet jet : jetCollectionEmbeddedMatched.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        Embedded_Matched_TrackOver100GeV.push_back(found);
      }
    } 
    jetCollectionEmbeddedMatched.addVector("Embedded_Matched_TrackOver100GeV", Embedded_Matched_TrackOver100GeV); 
    
    vector<double> z1_theta2_embedded_matched;      z1_theta2_embedded_matched.reserve(jetCollectionEmbeddedMatched.getJet().size()); 
    for(PseudoJet jet : jetCollectionEmbeddedMatched.getJet()) {
        z1_theta2_embedded_matched.push_back(Angularity_z1_theta2.result(jet));
    }
    jetCollectionEmbeddedMatched.addVector("Embedded_Matched_z1_theta2", z1_theta2_embedded_matched);
    
    trw.addCollection("Embedded_Matched_",     jetCollectionEmbeddedMatched);
    
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
