#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiaEvent.hh"
#include "include/extraInfo.hh"

#include "PU14/CmdLine.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv)
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);
 
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  unsigned int nEvent = cmdline.value<unsigned int>("-nev",1);  // first argument: command line option; second argument: default value

  //event generator settings
  double       ptHat = cmdline.value<double>("-pthat",20);
  unsigned int tune  = cmdline.value<int>("-tune",14);
  unsigned int job   = cmdline.value<int>("-job",0);
  unsigned int seed  = cmdline.value<int>("-seed",0);

  std::cout << "generating " << nEvent << " events with pthat = " << ptHat << " and tune = " << tune << std::endl;  

  pythiaEvent pyt(ptHat, tune, -0.95, 0.95, false, false, seed);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");
  TString outFileName = Form("%s/PythiaEventsTune%dPtHat%.0f_%i.pu14",dir,tune,ptHat,job);
  
  fout.open(outFileName.Data());
  
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------

    fout << "# event " << ie << "\n";

    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    std::vector<fastjet::PseudoJet> partons = pyt.getPartonList();
    for(fastjet::PseudoJet p : partons) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
   
    for(fastjet::PseudoJet p : particlesSig) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
    fout << "end\n";
  }

  pyt.getStat();
  
  fout.close();

  std::cout << "\n Finished generating PYTHIA events" << std::endl;
    
}
