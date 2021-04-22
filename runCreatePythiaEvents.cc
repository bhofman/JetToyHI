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
 
  // Number of events, generated and listed ones.
  //unsigned int nEvent    = 10000;

  //event generator settings
  double       ptHat = cmdline.value<double>("-pthat",120);//120.;
  unsigned int tune  = cmdline.value<int>("-tune",14);

  std::cout << "generating " << nEvent << " events with pthat = " << ptHat << " and tune = " << tune << std::endl;  

  pythiaEvent pyt(ptHat, tune, -3.0, 3.0);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString outFileName = Form("%s/PythiaEventsTune%dPtHat%.0f.pu14",dir,tune,ptHat);
  
  fout.open(outFileName.Data());
  
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);

    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------

    fout << "# event " << ie << "\n";

    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    std::vector<fastjet::PseudoJet> partons = pyt.getPartonList();
    int moeder_gluon = 0 ; // Hierin houden we bij of moeder een gluon is
    int dochter_cb = 0 ; //   Is dochter een charm of beauty??

  for(fastjet::PseudoJet p : partons) {

    const int & pdgid = p.user_info<extraInfo>().pdg_id();
    const int & vtx   = p.user_info<extraInfo>().vertex_number();
    std::cout<<vtx<<std::endl;
    if (vtx == -1) {// We zijn een moeder tegengekomen
      if (pdgid == 21){ // Moeder is ook een gluon 
       moeder_gluon = 1 ; } // We hebben een moeder gluon gezien
      else {
        moeder_gluon = 0 ; } // moeder is geen gluon
    }
    if (vtx == -2) { // We zijn een dochter tegengekomen
      if (abs(pdgid) == 4 || abs(pdgid) == 5){ // dochter is charm of beauty
        dochter_cb = 1 ; }
      else{
        dochter_cb = 0 ; }
      if (moeder_gluon == 1 && dochter_cb == 1) { continue; }// we hebben wat we willen
      }
    }

    if(moeder_gluon==1 && dochter_cb==1){
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
  }
    fout << "end\n";
  }

  fout.close();

  std::cout << "\n Finished generating PYTHIA events" << std::endl;
    
}
