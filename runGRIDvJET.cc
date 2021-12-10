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
#include "include/Angularity.hh"
#include "include/dyGroomer.hh"

#include "include/csSubtractor.hh"
#include "include/csSubFullEventIterative.hh"

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
    jetCollection jetCollection_Sig(sorted_by_pt(jet_selector(csSig.inclusive_jets(60.)))); // Inclusive jets to take a jets with pt over (pt_min)

    //---------------------------------------------------------------------------
    //   background subtraction 
    //---------------------------------------------------------------------------
    
    //GRID rho
    csSubFullEventIterative csSubFull_GRID( {2.} , {0.1}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull_GRID.setInputParticles(particlesMerged);
    csSubFull_GRID.setMaxEta(jetRapMax);
    csSubFull_GRID.setBackgroundGrid();
    fastjet::ClusterSequenceArea fullSig_GRID(csSubFull_GRID.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection csFullJets_GRID(sorted_by_pt(jet_selector(fullSig_GRID.inclusive_jets(0.)))); 

    std::vector<double> rhoFull_GRID, rhomFull_GRID;
    double rho_GRID = csSubFull_GRID.getRho();
    double rhom_GRID = csSubFull_GRID.getRhoM();

    rhoFull_GRID.push_back(rho_GRID);  
    rhomFull_GRID.push_back(rhom_GRID); 
    trw.addCollection("rhoFull_GRID",         rhoFull_GRID);
    trw.addCollection("rhomFull_GRID",        rhomFull_GRID);  

    //JET rho
    csSubtractor csSubFull_Jet(R, 1., -1, 0.005,ghostRapMax,jetRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull_Jet.setInputParticles(particlesMerged);
    fastjet::ClusterSequenceArea fullSig_JET(csSubFull_Jet.doSubtraction(), jet_def, area_def);
    jetCollection csFullJets_JET(sorted_by_pt(jet_selector(fullSig_JET.inclusive_jets(0.)))); 

    std::vector<double> rhoFull_Jet, rhomFull_Jet;
    double rho_Jet = csSubFull_Jet.getRho();
    double rhom_Jet = csSubFull_Jet.getRhoM();

    rhoFull_Jet.push_back(rho_Jet);  
    rhomFull_Jet.push_back(rhom_Jet); 
    trw.addCollection("rhoFull_Jet",          rhoFull_Jet);
    trw.addCollection("rhomFull_Jet",        rhomFull_Jet);  

    std::vector<double> rhoDIFF;
    std::vector<double> rhomDIFF;
    rhoDIFF.push_back((rho_GRID-rho_Jet)/rho_Jet);
    rhomDIFF.push_back((rhom_GRID-rhom_Jet)/rhom_Jet);

    trw.addCollection("rhoDIFF",         rhoDIFF);
    trw.addCollection("rhomDIFF",        rhomDIFF); 
    
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

  TFile *fout = new TFile(cmdline.value<string>("-output", "GRIDvJETcs.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}

// 
