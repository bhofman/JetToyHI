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
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(25.)))); // Inclusive jets to take a jets with pt over (pt_min)

    /*
    //calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);
    */
    
    //---------------------------------------------------------------------------
    //   jet clustering of signal+background jets
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csRaw(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionRaw(sorted_by_pt(jet_selector(csRaw.inclusive_jets(25.))));

    //calculate some angularities
    vector<double> widthRaw; widthRaw.reserve(jetCollectionRaw.getJet().size());
    vector<double> pTDRaw;   pTDRaw.reserve(jetCollectionRaw.getJet().size());
    for(PseudoJet jet : jetCollectionRaw.getJet()) {
      widthRaw.push_back(width.result(jet));
      pTDRaw.push_back(pTD.result(jet));
    }

    //match Raw(=unsubtracted) jets to signal jets
    jetMatcher jmRaw(R);
    jmRaw.setBaseJets(jetCollectionRaw);
    jmRaw.setTagJets(jetCollectionSig);
    jmRaw.matchJets();
    jmRaw.reorderedToTag(jetCollectionRaw);

    jetCollectionRaw.addVector("widthRaw", widthRaw);
    jetCollectionRaw.addVector("pTDRaw", pTDRaw);
    trw.addCollection("rawJet",       jetCollectionRaw);

    //---------------------------------------------------------------------------
    //   background subtraction Jet-by-Jet
    //---------------------------------------------------------------------------
    if(JETJET) {
    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    csSubtractor csSub(R, 0., -1, 0.005,ghostRapMax,jetRapMax);  // Rjet, alpha, rParam, ghA, ghostRapMax, jetRapMax
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());
    
    //Background densities used by constituent subtraction
    std::vector<double> rho;
    std::vector<double> rhom;
    rho.push_back(csSub.getRho());    
    rhom.push_back(csSub.getRhoM());  

    //match CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();
    jmCS.reorderedToTag(jetCollectionCS);
    
    trw.addCollection("csJet",        jetCollectionCS);
    trw.addCollection("csJetRho",         rho);
    trw.addCollection("csJetRhom",        rhom);

    std::vector<double> ptPullJet; ptPullJet.reserve(jetCollectionSig.getJet().size());
    std::vector<double> mPullJet; mPullJet.reserve(jetCollectionSig.getJet().size());

    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPullJet.push_back((jetCollectionCS.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionSig.getJet()[i].pt()));
      mPullJet.push_back((jetCollectionCS.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPullJet",        ptPullJet);
    trw.addCollection("mPullJet",        mPullJet);
    }
    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT
    //---------------------------------------------------------------------------
    if(FULL_EVENT) {
    //We want to substract for full event instead:
    csSubtractorFullEvent csSubFull( 0., .4, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection jetCollectionCSFull(sorted_by_pt(jet_selector(fullSig.inclusive_jets(25.)))); 

    //Background densities used by constituent subtraction
    std::vector<double> rhoFull;
    std::vector<double> rhomFull;
    rhoFull.push_back(csSubFull.getRho());  
    rhomFull.push_back(csSubFull.getRhoM()); 

    //match CS FULL jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSFull);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCSFull);

    trw.addCollection("csFull",        jetCollectionCSFull);
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);

    std::vector<double> ptPullFull; ptPullFull.reserve(jetCollectionSig.getJet().size());
    std::vector<double> mPullFull; mPullFull.reserve(jetCollectionSig.getJet().size());

    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPullFull.push_back((jetCollectionCSFull.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionSig.getJet()[i].pt()));
      mPullFull.push_back((jetCollectionCSFull.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPullFull",        ptPullFull);
    trw.addCollection("mPullFull",        mPullFull);
    }
    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT ITERATIVE
    //---------------------------------------------------------------------------
    if(FULL_EVENT_ITERATIVE) {
    //We want to substract for full event instead:
    csSubFullEventIterative csSubFull( {0.,0.} , {.2,.2}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    csSubFull.setBackground();
    //This line crashes:
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection jetCollectionCSIter(sorted_by_pt(jet_selector(fullSig.inclusive_jets(25.)))); 

    //match CS FULL jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSIter);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCSIter);

    //Background densities used by constituent subtraction
    std::vector<double> rhoIter;
    std::vector<double> rhomIter;
    rhoIter.push_back(csSubFull.getRho());  
    rhomIter.push_back(csSubFull.getRhoM()); 

    trw.addCollection("csIter",        jetCollectionCSIter);
    trw.addCollection("csIterRho",         rhoIter);
    trw.addCollection("csIterRhom",        rhomIter);

    std::vector<double> ptPullIter; ptPullIter.reserve(jetCollectionSig.getJet().size());
    std::vector<double> mPullIter; mPullIter.reserve(jetCollectionSig.getJet().size());

    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPullIter.push_back((jetCollectionCSIter.getJet()[i].pt()-jetCollectionSig.getJet()[i].pt())/(jetCollectionSig.getJet()[i].pt()));
      mPullIter.push_back((jetCollectionCSIter.getJet()[i].m()-jetCollectionSig.getJet()[i].m())/(jetCollectionSig.getJet()[i].m()));
    }

    trw.addCollection("ptPullIter",        ptPullIter);
    trw.addCollection("mPullIter",        mPullIter);
    }
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    //trw.addPartonCollection("partons",       partons);
    trw.addCollection("sigJet",        jetCollectionSig);
    //trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);
 
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultConstituentSubtraction.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}