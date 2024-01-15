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
  double R                   = 0.2;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;     
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 0.7;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);
  //Selector jet_selector = SelectorAbsEtaMax(jetRapMax);

  Angularity Angularity_z1_theta2(2.0,1.0,R);
  
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle((nEvent == -1 ? 8 : -1));

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

    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {0.0} , {0.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMergedAll);
    csSubFull.setMaxEta(1.5);
    fastjet::ClusterSequenceArea fullSig(csSubFull.Subtract(), jet_def, area_def);
    jetCollection jetCollectionCS_Sig(sorted_by_pt(jet_selector(fullSig.inclusive_jets(10.)))); 

    vector<double> TrackOver100GeV;      TrackOver100GeV.reserve(jetCollectionCS_Sig.getJet().size());
    double found;
    for(fastjet::PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      if(jet.has_constituents()) {
        found = 0;
        for(fastjet::PseudoJet constituent : jet.constituents()) {

    if (constituent.perp() > 100.) {
                found = constituent.perp();
                break;
            }
        }
        TrackOver100GeV.push_back(found);
      }
    } 
    jetCollectionCS_Sig.addVector("TrackOver100GeV", TrackOver100GeV);  

    //calculate some angularities
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionCS_Sig.getJet().size());
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
    }
    jetCollectionCS_Sig.addVector("z1_theta2", z1_theta2);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("",     jetCollectionCS_Sig);
    
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
