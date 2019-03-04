////////////////////////////////////////////////////////////////////////
// Class:       TriggerValidation
// Plugin Type: analyzer (art v3_01_02)
// File:        TriggerValidation_module.cc
//
// Generated at Sun Mar  3 21:45:37 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "TTree.h"

class TriggerValidation;


class TriggerValidation : public art::EDFilter {
public:
  explicit TriggerValidation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerValidation(TriggerValidation const&) = delete;
  TriggerValidation(TriggerValidation&&) = delete;
  TriggerValidation& operator=(TriggerValidation const&) = delete;
  TriggerValidation& operator=(TriggerValidation&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  int _run, _sub, _evt;
  unsigned int _time_low;
  unsigned int _time_high;

  std::string fOpHitProducer;

  // TTree storing pulse-by-pulse information
  TTree* _pulse_tree;
  int   _ch;
  float _time;
  float _pe;
  float _pe_amp_raw;
  float _pe_amp_gain;
  float _pe_area_raw;
  float _pe_area_gain;
  float _gain_area;
  float _gain_amp;
  float _largest_pe; // largest PE in event

};


TriggerValidation::TriggerValidation(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  fOpHitProducer = p.get< std::string >("OpHitProducer", "ophitBeam");

}

bool TriggerValidation::filter(art::Event& e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  _time_low  = e.time().timeLow();
  _time_high = e.time().timeHigh();

  // grab PMT calibrations if needed
  art::ServiceHandle<geo::Geometry> geo;
  const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();

  std::vector<float> chgain_amp_v = std::vector<float>{};
  std::vector<float> chgain_area_v = std::vector<float>{};
  for (size_t ch=0; ch < 32; ch++) {
    //std::cout << "querying gain at ch " << ch;
    chgain_amp_v.push_back( gain_provider.ExtraInfo(ch).GetFloatData("amplitude_gain") );
    chgain_area_v.push_back( gain_provider.Gain(ch) );
    //std::cout << "\t amp : " << chgain_amp_v.at(chgain_amp_v.size()-1)
    //	      << "\t area : " << chgain_area_v.at(chgain_area_v.size()-1) << std::endl;
  }

  // load optical flashes
  art::Handle<std::vector<recob::OpHit> > ophit_h;
  e.getByLabel( fOpHitProducer, ophit_h );

  _largest_pe = 0;
  for (size_t h=0; h < ophit_h->size(); h++) {
    auto const& ophit = ophit_h->at(h);
    if (ophit.PE() > _largest_pe) { _largest_pe = ophit.PE(); }
  }

  for (size_t h=0; h < ophit_h->size(); h++) {

    auto const& ophit = ophit_h->at(h);
    
    _ch         = ophit.OpChannel();

    if (_ch >= 32) continue;

    _time       = ophit.PeakTime();

    _gain_area = chgain_area_v[_ch];
    _gain_amp  = chgain_amp_v[_ch];

    _pe = ophit.PE();

    _pe_amp_raw   = ophit.Amplitude();
    _pe_area_raw  = ophit.Area();
    _pe_amp_gain  = _pe_amp_raw  * 20.  / _gain_amp;
    _pe_area_gain = _pe_area_raw * 130. / _gain_area;

    _pulse_tree->Fill();

  }// for all optical hits
    
  return true;
}

void TriggerValidation::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  _pulse_tree = tfs->make<TTree>("pulse_tree","One entry per recorder pulse");
  _pulse_tree->Branch ("_run",&_run,"run/I");
  _pulse_tree->Branch("_sub",&_sub,"sub/I");
  _pulse_tree->Branch("_evt",&_evt,"evt/I");
  _pulse_tree->Branch("_time_low" ,&_time_low ,"time_low/i" );
  _pulse_tree->Branch("_time_high",&_time_high,"time_high/i");
  _pulse_tree->Branch("_ch", &_ch, "ch/I" );
  _pulse_tree->Branch("_time",        &_time,        "time/F"        );
  _pulse_tree->Branch("_gain_amp",    &_gain_amp,    "gain_amp/F"    );
  _pulse_tree->Branch("_gain_area",   &_gain_area,   "gain_area/F"   );
  _pulse_tree->Branch("_pe",          &_pe,          "pe/F"          );
  _pulse_tree->Branch("_pe_amp_raw",  &_pe_amp_raw,  "pe_amp_raw/F"  );
  _pulse_tree->Branch("_pe_amp_gain", &_pe_amp_gain, "pe_amp_gain/F" );
  _pulse_tree->Branch("_pe_area_raw", &_pe_area_raw, "pe_area_raw/F" );
  _pulse_tree->Branch("_pe_area_gain",&_pe_area_gain,"pe_area_gain/F");
  _pulse_tree->Branch("_largest_pe",  &_largest_pe,  "largest_pe/F"  );

}

void TriggerValidation::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerValidation)
