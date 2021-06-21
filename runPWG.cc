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
    jetCollection jetCollection_Sig(sorted_by_pt(jet_selector(csSig.inclusive_jets(130.)))); // Inclusive jets to take a jets with pt over (pt_min)

    //---------------------------------------------------------------------------
    //   background subtraction 
    //---------------------------------------------------------------------------
    /*
    //run jet-by-jet constituent subtraction on mixed (hard+UE) event
    csSubtractor csSub(R, 1., -1, 0.005,ghostRapMax,jetRapMax);  // Rjet, alpha, rParam, ghA, ghostRapMax, jetRapMax
    csSub.setInputParticles(particlesMerged);
    jetCollection csFullJets(csSub.doSubtraction());
    */
    
    //We want to substract for full event iterative:
    csSubFullEventIterative csSubFull( {2.,2.} , {0.1,.075}, 0.005,ghostRapMax);  // alpha, rParam, ghA, ghRapMax
    csSubFull.setInputParticles(particlesMerged);
    csSubFull.setMaxEta(3.);
    csSubFull.setBackgroundGrid();
    fastjet::ClusterSequenceArea fullSig(csSubFull.doSubtractionFullEvent(), jet_def, area_def);
    jetCollection csFullJets(sorted_by_pt(jet_selector(fullSig.inclusive_jets(0.)))); 
    
    // Make sure our groomed jets have constituents
    std::vector<fastjet::PseudoJet> csFullJetsClean;
    for(fastjet::PseudoJet jet : csFullJets.getJet()) {
      if(jet.has_constituents()){
        csFullJetsClean.push_back(jet);
      }
    }
    jetCollection jetCollectionCS_Sig(csFullJetsClean);
    
    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCS_Sig);
    jmCSFull.setTagJets(jetCollection_Sig);
    jmCSFull.matchJets();
    jmCSFull.reorderedToTag(jetCollectionCS_Sig);
    
    /*
    //Background densities used by constituent subtraction
    std::vector<double> rhoFull;
    std::vector<double> rhomFull;
    rhoFull.push_back(csSubFull.getRho());  
    rhomFull.push_back(csSubFull.getRhoM()); 
    trw.addCollection("csFullRho",         rhoFull);
    trw.addCollection("csFullRhom",        rhomFull);  
    */
    std::vector<double> ptPull; ptPull.reserve(jetCollection_Sig.getJet().size());
    std::vector<double> mPull; mPull.reserve(jetCollection_Sig.getJet().size());
    for (unsigned int i = 0; i < jetCollection_Sig.getJet().size(); i++) {
      ptPull.push_back((jetCollectionCS_Sig.getJet()[i].pt()-jetCollection_Sig.getJet()[i].pt())/(jetCollectionCS_Sig.getJet()[i].pt()+jetCollection_Sig.getJet()[i].pt()));
      mPull.push_back((jetCollectionCS_Sig.getJet()[i].m()-jetCollection_Sig.getJet()[i].m())/(jetCollectionCS_Sig.getJet()[i].m()+jetCollection_Sig.getJet()[i].m()));
    }

    //---------------------------------------------------------------------------
    //   SOFTDROP Groom the CS jets
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSigBeta00Z01_CS(0.1, 0.0, R);
    jetCollection jetCollectionCS_SD(sdgSigBeta00Z01_CS.doGrooming(csFullJetsClean));
    softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    jetCollection jetCollection_SD(sdgSigBeta00Z01.doGrooming(jetCollection_Sig));

    //match CSFull jets to signal jets
    jetMatcher jmCSFullSD(R);
    jmCSFullSD.setBaseJets(jetCollectionCS_SD);
    jmCSFullSD.setTagJets(jetCollection_SD);
    jmCSFullSD.matchJets();
    jmCSFullSD.reorderedToTag(jetCollectionCS_SD);
    
    std::vector<double> ptPull_SD; ptPull_SD.reserve(jetCollection_SD.getJet().size());
    std::vector<double> mPull_SD; mPull_SD.reserve(jetCollection_SD.getJet().size());
    for (unsigned int i = 0; i < jetCollection_SD.getJet().size(); i++) {
      ptPull_SD.push_back((jetCollectionCS_SD.getJet()[i].pt()-jetCollection_SD.getJet()[i].pt())/(jetCollectionCS_SD.getJet()[i].pt()+jetCollection_SD.getJet()[i].pt()));
      mPull_SD.push_back((jetCollectionCS_SD.getJet()[i].m()-jetCollection_SD.getJet()[i].m())/(jetCollectionCS_SD.getJet()[i].m()+jetCollection_SD.getJet()[i].m()));
    }
    
    trw.addCollection("SIG_",        jetCollection_Sig);
    trw.addCollection("csSIG_",        jetCollectionCS_Sig);

    trw.addCollection("SD_",      jetCollection_SD);
    trw.addCollection("csSD_",      jetCollectionCS_SD);

    trw.addCollection("ptPull",        ptPull);
    trw.addCollection("mPull",        mPull);

    trw.addCollection("ptPull_SD",        ptPull_SD);
    trw.addCollection("mPull_SD",        mPull_SD);
    
    //---------------------------------------------------------------------------
    //   pt of constituents
    //---------------------------------------------------------------------------
    
    std::vector<double>  cons_sig_pt, cons_sigCS_pt, cons_SD_pt, cons_SDCS_pt;
    std::vector<double>  cons_sig_dr, cons_sigCS_dr, cons_SD_dr, cons_SDCS_dr;

    for(fastjet::PseudoJet jet : jetCollection_Sig.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_sig_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_sig_dr.push_back(DeltaR);
        }
      }
    } 
    
    for(fastjet::PseudoJet jet : jetCollectionCS_Sig.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_sigCS_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_sigCS_dr.push_back(DeltaR);
        }
      }
    }
    
    for(fastjet::PseudoJet jet : jetCollection_SD.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_SD_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_SD_dr.push_back(DeltaR);
        }
      }
    }

    for(fastjet::PseudoJet jet : jetCollectionCS_SD.getJet()) {
      if(jet.has_constituents()) {
        for(fastjet::PseudoJet constituent : jet.constituents()) {
          cons_SDCS_pt.push_back(constituent.perp());
          double DeltaR = std::sqrt(constituent.squared_distance(jet));
          cons_SDCS_dr.push_back(DeltaR);
        }
      }
    }

    trw.addCollection("cons_sig_pt",        cons_sig_pt);
    trw.addCollection("cons_sigCS_pt",        cons_sigCS_pt);
    trw.addCollection("cons_SD_pt",        cons_SD_pt);
    trw.addCollection("cons_SDCS_pt",        cons_SDCS_pt);

    trw.addCollection("cons_sig_dr",        cons_sig_dr);
    trw.addCollection("cons_sigCS_dr",        cons_sigCS_dr);
    trw.addCollection("cons_SD_dr",        cons_SD_dr);
    trw.addCollection("cons_SDCS_dr",        cons_SDCS_dr);    
    
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

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetCHECK.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}

// 
