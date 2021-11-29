#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/extraInfo.hh"
#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
//#include "include/Angularity.hh"
//#include "include/dyGroomer.hh"

//#include "include/csSubtractor.hh"
//#include "include/csSubFullEventIterative.hh"

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
  //fastjet::JetDefinition jet_def_ca(cambridge_algorithm, 999.); // Should no longer be needed

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);
    
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

    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    
    //std::cout << "#merged: " << particlesMerged.size() << "  signal: " << particlesSig.size() << "  bkg: " << particlesBkg.size() << std::endl;

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(130.)))); // Inclusive jets to take a jets with pt over (pt_min)
    
    std::vector<double>  jetMultiplicity;
    std::vector<double>  constituentPt;
    std::vector<double>  constituentPhi;
    std::vector<double>  constituentRap;

    std::vector<double>  constituentDeltaR;
    std::vector<double>  constituentEEoverQ2;

    for (long unsigned int i = 0; i < jetCollectionSig.getJet().size(); ++i) {
      fastjet::PseudoJet jet = jetCollectionSig.getJet()[i];
      jetMultiplicity.push_back(jet.constituents().size());

      if(jet.has_constituents()) {
        for (long unsigned int i = 0; i < jet.constituents().size(); ++i) {
          fastjet::PseudoJet constituent = jet.constituents()[i];
          constituentPt.push_back(constituent.perp());
          constituentPhi.push_back(constituent.phi());
          constituentRap.push_back(constituent.rap());
            for (long unsigned int j = i+1; j < jet.constituents().size(); ++j) {
              fastjet::PseudoJet second_constituent = jet.constituents()[j];
              constituentDeltaR.push_back(constituent.delta_R(second_constituent));
              constituentEEoverQ2.push_back((constituent.perp()*second_constituent.perp())/(jet.perp()*jet.perp()));
            }
        }
      }
    }

    trw.addCollection("sigJet",        jetCollectionSig);

    trw.addCollection("jetMultiplicity", jetMultiplicity);

    trw.addCollection("constituentPt",        constituentPt);
    trw.addCollection("constituentRap",        constituentRap);
    trw.addCollection("constituentPhi",        constituentPhi);

    trw.addCollection("constituentDeltaR",        constituentDeltaR);
    trw.addCollection("constituentEEoverQ2",        constituentEEoverQ2);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetEEC.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}

// 