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

#include "include/csSubtractorFullEvent.hh"

using namespace std;
using namespace fastjet;

// ./runSimpleJetAnalysis -hard /Users/mverweij/mnt/eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

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

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesMerged = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract hard partons that initiated the jets
    fastjet::Selector parton_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> partons = parton_selector(particlesMerged);
    
    // select final state particles from hard event only
    vector<PseudoJet> particlesBkg, particlesSig;
    SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event

    //---------------------------------------------------------------------------
    //   jet clustering --- SIGNAL ONLY
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(25.))));

    //calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);

    //---------------------------------------------------------------------------
    //   jet clustering --- SIG + BKG 
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csRaw(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionRaw(sorted_by_pt(jet_selector(csRaw.inclusive_jets(0.))));

    //match Raw(=unsubtracted) jets to signal jets, so jet samples correspond to each other
    //Otherwise many more jets in sig + bkg
    jetMatcher jmRaw(R);
    jmRaw.setBaseJets(jetCollectionRaw);
    jmRaw.setTagJets(jetCollectionSig);
    jmRaw.matchJets();
    jmRaw.reorderedToTag(jetCollectionRaw);

    //---------------------------------------------------------------------------
    //   background subtraction FULL EVENT
    //---------------------------------------------------------------------------

    csSubtractorFullEvent csSubFull( 0., .25, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection jetCollectionCS(sorted_by_pt(jet_selector(fullSig.inclusive_jets(0.)))); 

    //match background subtracted jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCS);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCS);

    // Remove jets which are empty after BKG subtraction
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : jetCollectionCS.getJet()) {
      if(jet.has_constituents()){
        csFullJetsClean.push_back(jet);
      }
    }
    jetCollection jetCollectionCS_noEmpty(csFullJetsClean);

    //calculate some angularities
    vector<double> widthSigCS; widthSigCS.reserve(jetCollectionCS_noEmpty.getJet().size());
    vector<double> pTDSigCS;   pTDSigCS.reserve(jetCollectionCS_noEmpty.getJet().size());
    for(PseudoJet jet : jetCollectionCS_noEmpty.getJet()) {
      widthSigCS.push_back(width.result(jet));
      pTDSigCS.push_back(pTD.result(jet));
    }
    jetCollectionCS_noEmpty.addVector("widthSigCS", widthSigCS);
    jetCollectionCS_noEmpty.addVector("pTDSigCS", pTDSigCS);

    //---------------------------------------------------------------------------
    //   SoftDrop Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());

    //SoftDrop grooming classic for signal + bkg subtracted jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01CS(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01CS(sdgSigBeta00Z01CS.doGrooming(jetCollectionCS_noEmpty));
    jetCollectionSigSDBeta00Z01CS.addVector("zgSigSDBeta00Z01CS",    sdgSigBeta00Z01CS.getZgs());
    jetCollectionSigSDBeta00Z01CS.addVector("ndropSigSDBeta00Z01CS", sdgSigBeta00Z01CS.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01CS.addVector("dr12SigSDBeta00Z01CS",  sdgSigBeta00Z01CS.getDR12());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addPartonCollection("partons",       partons);

    // Signal jets:
    trw.addCollection("sigJet",        jetCollectionSig);
    // Sig + BKG jets:
    trw.addCollection("rawJet",       jetCollectionRaw);
    // BKG subtracted jets:
    trw.addCollection("csFull",        jetCollectionCS_noEmpty);
    
    // SoftDropped signal jets
    trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);
    // SoftDropped BKG subtracted jets
    trw.addCollection("sigJetSDBeta00Z01CS",      jetCollectionSigSDBeta00Z01CS);

    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultSimpleJetAnalysis.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
