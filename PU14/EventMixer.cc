//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Yi Chen (FHead) https://github.com/FHead
// https://github.com/FHead/JetToyHI/blob/49d264cc304602341e56a315f0a9dbd768016f57/PU14/EventMixer.cc
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "PU14.hh"
#include "helpers.hh"
#include <cstdlib>
#include <iostream>

using namespace std;

//----------------------------------------------------------------------
EventMixer::EventMixer(CmdLine * cmdline) : _cmdline(cmdline) {

  _hard_name   = _cmdline->value<string>("-hard");
  _hard_type   = _cmdline->value<string>("-hardtype", "PU14");
  _hard_varname   = _cmdline->value<string>("-hardvarname", "particle_gen");
  _hard_treename   = _cmdline->value<string>("-hardtreename", "AliAnalysisTaskTrackSkim_tree");

  _reco_name   = _cmdline->value<string>("-reco","");
  _reco_type   = _cmdline->value<string>("-recotype", "PU14");
  _reco_varname   = _cmdline->value<string>("-recovarname", "particle_data");
  _reco_treename   = _cmdline->value<string>("-recotreename", "AliAnalysisTaskTrackSkim_tree");

  _pileup_name = _cmdline->value<string>("-pileup", "");
  _pileup_type = _cmdline->value<string>("-pileuptype", "PU14");
  _pileup_varname = _cmdline->value<string>("-pileupvarname", "particle_");
  _pileup_treename = _cmdline->value<string>("-pileuptreename", "AliAnalysisTaskTrackSkim_tree");

  // setting the multiplicity of pileup events (background HI)
  //
  //  -npu <npu>  : fixed <npu> number of PU vertices - default to 1

  // fixed (the default)
  _npu = _cmdline->value("-npu", 1);
  if(_npu > 1)
  {
     cerr << "WARNING: number of background event requested = " << _npu << endl;
     cerr << "   make sure you actually want that!" << endl;
  }

  _massless = _cmdline->present("-massless");

  if (_cmdline->present("-chs")) {
    set_chs_rescaling_factor(1e-60);
  } else {
    set_chs_rescaling_factor(1.0); // this effectively turns off CHS
  }

  _hard  .reset(new EventSource(_hard_name  , _hard_type , _hard_varname , _hard_treename));
  //_hard->Recycle = false;

  if (_reco_name.empty()){
    cerr << "INFO: no detector level requested" << endl;
    _reco.reset();
  } else {
    _reco.reset(new EventSource(_reco_name, _reco_type , _reco_varname , _reco_treename));
    //_reco->Recycle = false;
  }

  if (_pileup_name.empty()){
    cerr << "INFO: no background requested" << endl;
    _pileup.reset();
    _npu=0;
  } else {
    _pileup.reset(new EventSource(_pileup_name, _pileup_type , _pileup_varname , _pileup_treename, true));
    //_pileup->Recycle = true;
  }
}

//----------------------------------------------------------------------
bool EventMixer::next_event() {
  _particles.resize(0);
  _particles.clear();
  _hard_event_weight = 1;
  _pu_event_weight = 1;
  
  // first get the hard event ( truth level )
  if (! _hard->append_next_event(_particles,_hard_event_weight,0)) return false;

  unsigned hard_size = _particles.size();

  // add detector level event if available
  if (_reco.get()){
      if (! _reco->append_next_event(_particles,_pu_event_weight,99)) std::cout << "Empty reco event" << std::endl; //return false;
  }
  
  // add pileup if available
  if (_pileup.get()){
    for (int i = 1; i <= _npu; i++) {
      if (! _pileup->append_next_event(_particles,_pu_event_weight,i,true)) std::cout << "Empty pileup event" << std::endl; //return false;
    }
  }

  // make particles massless if requested
  if (_massless){
    _particles = MasslessTransformer()(_particles);
  }

  // apply CHS rescaling factor if requested
  if (chs_rescaling_factor() != 1.0) {
    for (unsigned i = hard_size; i < _particles.size(); i++) {
      if (_particles[i].user_info<PU14>().charge() != 0) _particles[i] *= _chs_rescaling_factor;
    }
  }

  return true;
}

//----------------------------------------------------------------------
string EventMixer::description() const {
  ostringstream ostr;
  ostr << "Event mixer using hard events from " << _hard_name;

  if (_npu > 0) {
    ostr << " and " << _npu << " pileup events from " << _pileup_name;
  
    if (chs_rescaling_factor() != 1.0) {
    ostr << " with CHS (rescaling charged PU by a factor "
         << chs_rescaling_factor() << ")";
    }
  }
  if (_massless){
    ostr << " and massless particles";
  }
  return ostr.str();
}
