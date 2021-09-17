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
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/csSubFullEventIterative.hh"

using namespace std;
using namespace fastjet;

// ./runAnalysis -hard samples/PythiaEventsTune14PtHat120.pu14 -pileup samples/ThermalEventsMult12000PtAv0.70.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  // first argument: command line option; second argument: default value
  int nEvent = cmdline.value<int>("-nev",1); 
  cout << "will run on " << nEvent <<"events"<<endl<<endl; 
  //bool verbose = cmdline.present("-verbose");

  bool JETJET = cmdline.present("-jet");
  bool FULL_EVENT = cmdline.present("-full");
  bool FULL_EVENT_ITERATIVE = cmdline.present("-iter"); 
  if(JETJET) {  cout << "Using method: Jet by Jet" << endl<<endl; }
  if(FULL_EVENT) {  cout << "Using method: Full event" << endl<<endl; }
  if(FULL_EVENT_ITERATIVE) {  cout << "Using method: Full Iterative" << endl<<endl; }
  
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
    //cout << iev << endl;
       
    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract hard partons that initiated the jets
    fastjet::Selector parton_selector = SelectorVertexNumber(-1);  
    vector<PseudoJet> partons = parton_selector(particlesMergedAll);
    
    // select final state particles from hard event only
    //vector<PseudoJet> particlesBkg, particlesSig;
    //SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event

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
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(40.))));
    trw.addCollection("sigJet",        jetCollectionSig);

    //---------------------------------------------------------------------------
    //   jet clustering of signal+background jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csRaw(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionRaw(sorted_by_pt(jet_selector(csRaw.inclusive_jets(0.))));

    //match Raw(=unsubtracted) jets to signal jets
    jetMatcher jmRaw(R);
    jmRaw.setBaseJets(jetCollectionRaw);
    jmRaw.setTagJets(jetCollectionSig);
    jmRaw.matchJets();
    jmRaw.reorderedToTag(jetCollectionRaw);

    trw.addCollection("rawJet",        jetCollectionRaw);
    
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01Sig(0.1, 0.0, R);

    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01Sig.doGrooming(jetCollectionSig));
    trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);

    //---------------------------------------------------------------------------
    //   constituents
    //---------------------------------------------------------------------------
    
    std::vector<double>  cons_sig_pt, cons_SD_pt, cons_raw_pt ;
    std::vector<double>  cons_sig_dr, cons_SD_dr, cons_raw_dr ;

    for(fastjet::PseudoJet jet : jetCollectionSig.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_sig_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_sig_dr.push_back(DeltaR);
        }
      }
    } 

    for(fastjet::PseudoJet jet : jetCollectionSigSDBeta00Z01.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_SD_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_SD_dr.push_back(DeltaR);
        }
      }
    }

    for(fastjet::PseudoJet jet : jetCollectionRaw.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_raw_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_raw_dr.push_back(DeltaR);
        }
      }
    } 

    trw.addCollection("cons_sig_pt",        cons_sig_pt);
    trw.addCollection("cons_SD_pt",        cons_SD_pt);
    trw.addCollection("cons_raw_pt",        cons_raw_pt);

    trw.addCollection("cons_sig_dr",        cons_sig_dr);
    trw.addCollection("cons_SD_dr",        cons_SD_dr);
    trw.addCollection("cons_raw_dr",        cons_raw_dr);

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

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyBasis.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}