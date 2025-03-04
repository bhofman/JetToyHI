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
#include "include/softDropGroomer.hh"
#include <TRandom.h>

#include "include/AliceFastSim.hh"
#include "include/thermalAlice.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);

  cout << "will run on " << nEvent << " events" << endl;

  TFile *fout = new TFile(cmdline.value<string>("-output", "test.root").c_str(), "RECREATE");

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

  double jetRapMax = 1.0;
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle((nEvent == -1 ? 7 : -1));

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

AliceFastSim fastSim = AliceFastSim();
  thermalAlice thrmEvent;

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

    fastSim.setInputEvent(particlesMergedAll);
    vector<PseudoJet> particlesTruth = fastSim.AliceAcceptance();
    vector<PseudoJet> particlesDetector = fastSim.AliceDetector();
    
    vector<fastjet::PseudoJet> particlesPileup = thrmEvent.createThermalEventAlice();
    //std::cout << "Pileup particles: " << particlesPileup.size() << std::endl;
    //trw.addCollection("particlesPileup_",      particlesPileup);

    //std::vector<double> Npileup;
    //Npileup.push_back(particlesPileup.size()); 
    //trw.addCollection("Npileup",         Npileup);
  
    vector<PseudoJet> particlesMerged = particlesDetector;
    particlesMerged.insert( particlesMerged.end(), particlesPileup.begin(), particlesPileup.end() );
    
    //---------------------------------------------------------------------------
    //   Truth level jets before any matching
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csTruth(particlesTruth, jet_def, area_def);
    jetCollection jetCollectionTruth(sorted_by_pt(jet_selector(csTruth.inclusive_jets(15.)))); 
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
    //   Embedded level jets before any matching
    //---------------------------------------------------------------------------
    csSubFullEventIterative csSubEmbedded( {0.0} , {0.1}, 0.005,ghostRapMax); // alpha, rParam, ghA, ghRapMax
    csSubEmbedded.setInputParticles(particlesMerged);
    csSubEmbedded.setMaxEta(1.0);
    fastjet::ClusterSequenceArea csEmbedded(csSubEmbedded.Subtract(), jet_def, area_def);
    jetCollection jetCollectionEmbedded(sorted_by_pt(jet_selector(csEmbedded.inclusive_jets(15.)))); 
    trw.addCollection("Embedded_",      jetCollectionEmbedded);
 
    std::vector<double> rho;
    rho.push_back(csSubEmbedded.getRho()); 
    trw.addCollection("rho",         rho);

    //              ## SD Groom jets ## 
    //SoftDrop grooming alice (zcut=0.2, beta=0)
    softDropGroomer Embedded_SDGroomer(0.2, 0.0, R);
    jetCollection Embedded_SD_(Embedded_SDGroomer.doGrooming(jetCollectionEmbedded));
    Embedded_SD_.addVector("Embedded_SD_zg",    Embedded_SDGroomer.getZgs());
    Embedded_SD_.addVector("Embedded_SD_ndrop", Embedded_SDGroomer.getNDroppedSubjets());
    Embedded_SD_.addVector("Embedded_SD_dr12",  Embedded_SDGroomer.getDR12());
    trw.addCollection("Embedded_SD_",      Embedded_SD_);    
    
    //---------------------------------------------------------------------------
    //   Matching of embedded jets to matched detector jets
    //---------------------------------------------------------------------------
    jetCollection jetCollectionEmbeddedMatched(jetCollectionEmbedded); 

    //match CSFull jets to signal jets
    jetMatcher jmEmbedded(0.6*R);
    jmEmbedded.setBaseJets(jetCollectionTruth);
    jmEmbedded.setTagJets(jetCollectionEmbeddedMatched);
    jmEmbedded.matchJets();
    jmEmbedded.reorderedToBase(jetCollectionEmbeddedMatched);

    // 100 GeV track finder
    vector<double> Embedded_Matched_TrackOver100GeV; int foundEmbeddedMatched;
    Embedded_Matched_TrackOver100GeV.reserve(jetCollectionEmbeddedMatched.getJet().size());
    for(fastjet::PseudoJet jet : jetCollectionEmbeddedMatched.getJet()) {
      if(jet.has_constituents()) {
        foundEmbeddedMatched = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {
            if (constituent.perp() > 100.) {
                foundEmbeddedMatched = constituent.perp();
                break;
            }
        }
        Embedded_Matched_TrackOver100GeV.push_back(foundEmbeddedMatched);
      }
      else { Embedded_Matched_TrackOver100GeV.push_back(0); }
    } 
    jetCollectionEmbeddedMatched.addVector("Embedded_Matched_TrackOver100GeV: ", Embedded_Matched_TrackOver100GeV); 
    trw.addCollection("Embedded_Matched_",     jetCollectionEmbeddedMatched);

    //              ## SD Groom jets ## 
    //SoftDrop grooming alice (zcut=0.2, beta=0)
    softDropGroomer Embedded_Matched_SDGroomer(0.2, 0.0, R);
    jetCollection Embedded_Matched_SD_(Embedded_Matched_SDGroomer.doGrooming(jetCollectionEmbeddedMatched));
    Embedded_Matched_SD_.addVector("Embedded_Matched_SD_zg",    Embedded_Matched_SDGroomer.getZgs());
    Embedded_Matched_SD_.addVector("Embedded_Matched_SD_ndrop", Embedded_Matched_SDGroomer.getNDroppedSubjets());
    Embedded_Matched_SD_.addVector("Embedded_Matched_SD_dr12",  Embedded_Matched_SDGroomer.getDR12());
    trw.addCollection("Embedded_Matched_SD_",      Embedded_Matched_SD_);

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
