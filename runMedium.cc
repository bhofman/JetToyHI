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

//hsdjdfkjfdwkjlsdafklnjsdfakljsdfakjlsdsdffsdsfdsfd
using namespace std;
using namespace fastjet;

// ./runAnalysis -hard samples/PythiaEventsTune14PtHat120.pu14 -pileup samples/ThermalEventsMult12000PtAv0.70.pu14 -nev 10

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  int user_pt = cmdline.value<int>("-pt",10); 

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

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

  Angularity mr(1.,0.,R);
  Angularity mr2(2.,0.,R);
  Angularity r2z(2.,1.,R);

  fastjet::contrib::OnePass_WTA_KT_Axes axes;
  fastjet::contrib::UnnormalizedMeasure unormbeta(1.0);
  fastjet::contrib::Nsubjettiness  nSub1_beta1(1, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub2_beta1(2, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub3_beta1(3, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub4_beta1(4, axes, unormbeta);
  fastjet::contrib::Nsubjettiness  nSub5_beta1(5, axes, unormbeta);
    
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

    // Jewel sub:
    fastjet::Selector dummy_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> particlesDummy = dummy_selector(particlesMergedAll);

    for(int i = 0; i < (int)particlesDummy.size(); i++)
    {
       if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000)
       {
          particlesDummy.erase(particlesDummy.begin() + i);
          i = i - 1;
       }
    }

    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);  // distance in eta-phi plane
    subtractor.set_max_distance(
        0.1);  // free parameter for the maximal allowed distance between particle i and ghost k
    subtractor.set_alpha(
        2.);  // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore
    subtractor.set_do_mass_subtraction();
    subtractor.set_remove_all_zero_pt_particles(true);

    std::vector<fastjet::PseudoJet> subtracted_particles = subtractor.do_subtraction(particlesMerged, particlesDummy);

    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------
    ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(user_pt)))); // te vrpef, tel dummy pt niet mee.
    
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    jetCollectionSigSDBeta00Z01.addVector("sig_SD_zg",    sdgSigBeta00Z01.getZgs());
    jetCollectionSigSDBeta00Z01.addVector("sig_SD_ndrop", sdgSigBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01.addVector("sig_SD_dr12",  sdgSigBeta00Z01.getDR12());
    trw.addCollection("sig_SD_",        jetCollectionSigSDBeta00Z01);
    
    //---------------------------------------------------------------------------
    //  
    //---------------------------------------------------------------------------
    
    jetCollection jetCollectionSigJewelSub(GetCorrectedJets(jetCollectionSig.getJet(), particlesDummy));

    softDropGroomer sdgSigSubBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollectionSigSDSubBeta00Z01(sdgSigSubBeta00Z01.doGroomingWithJewelSub(jetCollectionSig,particlesDummy));
    jetCollectionSigSDSubBeta00Z01.addVector("jewelsub_SD_zg",    sdgSigSubBeta00Z01.getZgs());
    jetCollectionSigSDSubBeta00Z01.addVector("jewelsub_SD_ndrop", sdgSigSubBeta00Z01.getNDroppedSubjets());
    jetCollectionSigSDSubBeta00Z01.addVector("jewelsub_SD_dr12",  sdgSigSubBeta00Z01.getDR12());
    trw.addCollection("jewelsub_SD_",      jetCollectionSigSDSubBeta00Z01);
    
    //---------------------------------------------------------------------------
    //   
    //---------------------------------------------------------------------------
    fastjet::ClusterSequenceArea csSigSub(subtracted_particles, jet_def, area_def);
    jetCollection jetCollectionCSSub(sorted_by_pt(jet_selector(csSigSub.inclusive_jets(1.)))); // Inclusive jets to take a jets with pt over (pt_min) 

    //match CSFull jets to signal jets
    jetMatcher jmCSFullSD(R);
    jmCSFullSD.setBaseJets(jetCollectionCSSub);
    jmCSFullSD.setTagJets(jetCollectionSig);
    jmCSFullSD.matchJets();
    jmCSFullSD.reorderedToTag(jetCollectionCSSub);
    
    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : jetCollectionCSSub.getJet()) {
      if(jet.has_constituents()){
        csFullJetsClean.push_back(jet);
      }
    }
    
    jetCollection jetCollectionCS_Sig(csFullJetsClean);

    softDropGroomer sdgSigBeta00Z01_cs(0.1, 0.0, R);
    jetCollection jetCollectionSigSDBeta00Z01_cs(sdgSigBeta00Z01_cs.doGrooming(jetCollectionCS_Sig));
    jetCollectionSigSDBeta00Z01_cs.addVector("cssub_SD_zg",    sdgSigBeta00Z01_cs.getZgs());
    jetCollectionSigSDBeta00Z01_cs.addVector("cssub_SD_ndrop", sdgSigBeta00Z01_cs.getNDroppedSubjets());
    jetCollectionSigSDBeta00Z01_cs.addVector("cssub_SD_dr12",  sdgSigBeta00Z01_cs.getDR12());
    trw.addCollection("cssub_SD_",      jetCollectionSigSDBeta00Z01_cs);

    //---------------------------------------------------------------------------
    //   
    //---------------------------------------------------------------------------
    
    std::vector<double> ptPull; ptPull.reserve(jetCollectionSig.getJet().size());
    for (unsigned int i = 0; i < jetCollectionSig.getJet().size(); i++) {
      ptPull.push_back((jetCollectionSigJewelSub.getJet()[i].pt()-jetCollectionCSSub.getJet()[i].pt())/jetCollectionSigJewelSub.getJet()[i].pt());
    }
    trw.addCollection("ptPull",        ptPull);
    
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported
    trw.addCollection("sig",        jetCollectionSig);
    trw.addCollection("cssub",      jetCollectionCSSub);
    trw.addCollection("jewelsub",      jetCollectionSigJewelSub);
  
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetBKG.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}