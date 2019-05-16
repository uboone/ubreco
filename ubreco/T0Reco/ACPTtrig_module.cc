////////////////////////////////////////////////////////////////////////
// Class:       ACPTtrig
// Plugin Type: producer (art v3_00_00)
// File:        ACPTtrig_module.cc
//
// Generated at Wed Jan  2 13:44:06 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
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
#include "lardata/Utilities/FindManyInChainP.h"

#include <memory>

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// save info associated to common optical filter
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <TTree.h>

class ACPTtrig;


class ACPTtrig : public art::EDFilter {
public:
  explicit ACPTtrig(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ACPTtrig(ACPTtrig const&) = delete;
  ACPTtrig(ACPTtrig&&) = delete;
  ACPTtrig& operator=(ACPTtrig const&) = delete;
  ACPTtrig& operator=(ACPTtrig&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  double fTrackLenMin;
  double fCathodeMin, fCathodeMax;
  double fAnodeMin, fAnodeMax;
  double fYMin, fYMax;
  double fBeamSpillStart, fBeamSpillEnd;
  std::string fTrackProducer, fFlashProducer, fOpDetWfmProducer, fCRTTagProducer, fCaloProducer;
  bool fSaveTree, fSaveWf;

  TTree* _tree;
  float _len;
  float _trk_beg_x, _trk_beg_y, _trk_beg_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  float _trk_beg_x_off, _trk_beg_y_off, _trk_beg_z_off;
  float _trk_end_x_off, _trk_end_y_off, _trk_end_z_off;
  float _trk_len; //corrected by SCE

  float _calorimetry_energy_plane_0; float _calorimetry_energy_plane_1; float _calorimetry_energy_plane_2; // float _residualrange; float _dedx;
  //float _dqdx; float _total_charge;
  std::vector<float> _residualrange_vector_plane_0;
  std::vector<float> _residualrange_vector_plane_1;
  std::vector<float> _residualrange_vector_plane_2;
  std::vector<float> _dedx_vector_plane_0;
  std::vector<float> _dedx_vector_plane_1;
  std::vector<float> _dedx_vector_plane_2;

  float _trk_charge_0;
  float _trk_charge_1;
  float _trk_charge_2;
  
  float _flash_t, _flash_pe, _flash_yw, _flash_yc, _flash_zw, _flash_zc;
  float  _opfilter_pe_beam, _opfilter_pe_beam_tot, _opfilter_pe_veto, _opfilter_pe_veto_tot;
  int _crttag;
  float _crtt0;
  int    _nflash;
  std::vector<double> _flash_pe_v;
  std::vector<std::vector<short>> _wf_v;

  int _run;
  int _subrun;
  int _event;

  void ResetTTree();

};


ACPTtrig::ACPTtrig(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{

  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::Track, anab::T0> >();
  produces< art::Assns <recob::Track, recob::OpFlash> >();

  fTrackProducer    = p.get<std::string>("TrackProducer");
  fFlashProducer    = p.get<std::string>("FlashProducer");
  fOpDetWfmProducer = p.get<std::string>("OpDetWfmProducer");
  fCRTTagProducer   = p.get<std::string>("CRTTagProducer","");
  fCaloProducer     = p.get<std::string>("CaloProducer");
  fTrackLenMin      = p.get<double>("TrackLenMin");
  fCathodeMin       = p.get<double>("CathodeMin");
  fCathodeMax       = p.get<double>("CathodeMax");
  fAnodeMin         = p.get<double>("AnodeMin");
  fAnodeMax         = p.get<double>("AnodeMax");
  fYMin             = p.get<double>("YMin");
  fYMax             = p.get<double>("YMax");
  fBeamSpillStart   = p.get<double>("BeamSpillStart");
  fBeamSpillEnd     = p.get<double>("BeamSpillEnd");
  fSaveTree         = p.get<bool>("SaveTree",false);
  fSaveWf           = p.get<bool>("SaveWf",false);

  _wf_v.resize(32);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","ACPT trigger tree");
  
  _tree->Branch("_run",&_run, "run/I");
  _tree->Branch("_subrun",&_subrun, "subrun/I");
  _tree->Branch("_event", &_event, "event/I");

  _tree->Branch("_trk_beg_x",&_trk_beg_x,"trk_beg_x/F");
  _tree->Branch("_trk_beg_y",&_trk_beg_y,"trk_beg_y/F");
  _tree->Branch("_trk_beg_z",&_trk_beg_z,"trk_beg_z/F");
  _tree->Branch("_trk_end_x",&_trk_end_x,"trk_end_x/F");
  _tree->Branch("_trk_end_y",&_trk_end_y,"trk_end_y/F");
  _tree->Branch("_trk_end_z",&_trk_end_z,"trk_end_z/F");
  _tree->Branch("_len",&_len,"len/F");
  
  _tree->Branch("_trk_beg_x_off",&_trk_beg_x_off,"trk_beg_x_off/F");
  _tree->Branch("_trk_beg_y_off",&_trk_beg_y_off,"trk_beg_y_off/F");
  _tree->Branch("_trk_beg_z_off",&_trk_beg_z_off,"trk_beg_z_off/F");
  _tree->Branch("_trk_end_x_off",&_trk_end_x_off,"trk_end_x_off/F");
  _tree->Branch("_trk_end_y_off",&_trk_end_y_off,"trk_end_y_off/F");
  _tree->Branch("_trk_end_z_off",&_trk_end_z_off,"trk_end_z_off/F");
  _tree->Branch("_trk_len",&_trk_len,"trk_len/F");

  _tree->Branch("_calorimetry_energy_plane_0",&_calorimetry_energy_plane_0,"calorimetry_energy_plane_0/F");
  _tree->Branch("_calorimetry_energy_plane_1",&_calorimetry_energy_plane_1,"calorimetry_energy_plane_1/F");
  _tree->Branch("_calorimetry_energy_plane_2",&_calorimetry_energy_plane_2,"calorimetry_energy_plane_2/F");

  // _tree->Branch("_residualrange",&_residualrange,"residualrange/F");
  // _tree->Branch("_dedx",&_dedx,"dedx/F");
  _tree->Branch("_residualrange_vector_plane_0","std::vector<float>",&_residualrange_vector_plane_0);
  _tree->Branch("_residualrange_vector_plane_1","std::vector<float>",&_residualrange_vector_plane_1);
  _tree->Branch("_residualrange_vector_plane_2","std::vector<float>",&_residualrange_vector_plane_2);
  _tree->Branch("_dedx_vector_plane_0","std::vector<float>",&_dedx_vector_plane_0);
  _tree->Branch("_dedx_vector_plane_1","std::vector<float>",&_dedx_vector_plane_1);
  _tree->Branch("_dedx_vector_plane_2","std::vector<float>",&_dedx_vector_plane_2);
//_tree->Branch("_dqdx",&_dqdx,"dqdx/F");
  //_tree->Branch("_total_charge",&_total_charge,"total_charge/F");

  _tree->Branch("_trk_charge_0",&_trk_charge_0,"trk_charge_0/F");
  _tree->Branch("_trk_charge_1",&_trk_charge_1,"trk_charge_1/F");
  _tree->Branch("_trk_charge_2",&_trk_charge_2,"trk_charge_2/F");

  _tree->Branch("_nflash",&_nflash,"nflash/I");
  _tree->Branch("_flash_t",&_flash_t,"flash_t/F");
  _tree->Branch("_flash_pe",&_flash_pe,"flash_pe/F");
  _tree->Branch("_flash_pe_v","std::vector<double>",&_flash_pe_v);
  _tree->Branch("_flash_yw",&_flash_yw,"flash_yw/F");
  _tree->Branch("_flash_zw",&_flash_zw,"flash_zw/F");
  _tree->Branch("_flash_yc",&_flash_yc,"flash_yc/F");
  _tree->Branch("_flash_zc",&_flash_zc,"flash_zc/F");
  _tree->Branch("_opfilter_pe_beam",&_opfilter_pe_beam,"opfilter_pe_beam/F");
  _tree->Branch("_opfilter_pe_beam_tot",&_opfilter_pe_beam_tot,"opfilter_pe_beam_tot/F");
  _tree->Branch("_opfilter_pe_veto",&_opfilter_pe_veto,"opfilter_pe_veto/F");
  _tree->Branch("_opfilter_pe_veto_tot",&_opfilter_pe_veto_tot,"opfilter_pe_veto_tot/F");

  _tree->Branch("_crttag",&_crttag,"crttag/I");
  _tree->Branch("_crtt0" ,&_crtt0 ,"crtt0/F");

  // save waveform
  if (fSaveWf) {
    _tree->Branch( "wf_00", "std::vector<short>", &(_wf_v[0])  );
    _tree->Branch( "wf_01", "std::vector<short>", &(_wf_v[1])  );
    _tree->Branch( "wf_02", "std::vector<short>", &(_wf_v[2])  );
    _tree->Branch( "wf_03", "std::vector<short>", &(_wf_v[3])  );
    _tree->Branch( "wf_04", "std::vector<short>", &(_wf_v[4])  );
    _tree->Branch( "wf_05", "std::vector<short>", &(_wf_v[5])  );
    _tree->Branch( "wf_06", "std::vector<short>", &(_wf_v[6])  );
    _tree->Branch( "wf_07", "std::vector<short>", &(_wf_v[7])  );
    _tree->Branch( "wf_08", "std::vector<short>", &(_wf_v[8])  );
    _tree->Branch( "wf_09", "std::vector<short>", &(_wf_v[9])  );
    _tree->Branch( "wf_10", "std::vector<short>", &(_wf_v[10])  );
    _tree->Branch( "wf_11", "std::vector<short>", &(_wf_v[11])  );
    _tree->Branch( "wf_12", "std::vector<short>", &(_wf_v[12])  );
    _tree->Branch( "wf_13", "std::vector<short>", &(_wf_v[13])  );
    _tree->Branch( "wf_14", "std::vector<short>", &(_wf_v[14])  );
    _tree->Branch( "wf_15", "std::vector<short>", &(_wf_v[15])  );
    _tree->Branch( "wf_16", "std::vector<short>", &(_wf_v[16])  );
    _tree->Branch( "wf_17", "std::vector<short>", &(_wf_v[17])  );
    _tree->Branch( "wf_18", "std::vector<short>", &(_wf_v[18])  );
    _tree->Branch( "wf_19", "std::vector<short>", &(_wf_v[19])  );
    _tree->Branch( "wf_20", "std::vector<short>", &(_wf_v[20])  );
    _tree->Branch( "wf_21", "std::vector<short>", &(_wf_v[21])  );
    _tree->Branch( "wf_22", "std::vector<short>", &(_wf_v[22])  );
    _tree->Branch( "wf_23", "std::vector<short>", &(_wf_v[23])  );
    _tree->Branch( "wf_24", "std::vector<short>", &(_wf_v[24])  );
    _tree->Branch( "wf_25", "std::vector<short>", &(_wf_v[25])  );
    _tree->Branch( "wf_26", "std::vector<short>", &(_wf_v[26])  );
    _tree->Branch( "wf_27", "std::vector<short>", &(_wf_v[27])  );
    _tree->Branch( "wf_28", "std::vector<short>", &(_wf_v[28])  );
    _tree->Branch( "wf_29", "std::vector<short>", &(_wf_v[29])  );
    _tree->Branch( "wf_30", "std::vector<short>", &(_wf_v[30])  );
    _tree->Branch( "wf_31", "std::vector<short>", &(_wf_v[31])  );
  }
}

bool ACPTtrig::filter(art::Event& e)
{

  ResetTTree();

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::Track, anab::T0> >       trk_t0_assn_v   ( new art::Assns<recob::Track, anab::T0>       );
  std::unique_ptr< art::Assns <recob::Track, recob::OpFlash> > trk_flash_assn_v( new art::Assns<recob::Track, recob::OpFlash> );

  // load Flash
  art::Handle<std::vector<recob::OpFlash> > flash_h;
  e.getByLabel(fFlashProducer,flash_h);

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    throw std::exception();
  }

  // load commont-optical-filter output
  art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
  art::InputTag fCommonOpFiltTag("opfiltercommonext");
  e.getByLabel(fCommonOpFiltTag, CommonOpticalFilter_h);

  _opfilter_pe_beam     = CommonOpticalFilter_h->PE_Beam();
  _opfilter_pe_beam_tot = CommonOpticalFilter_h->PE_Beam_Total();
  _opfilter_pe_veto     = CommonOpticalFilter_h->PE_Veto();
  _opfilter_pe_veto_tot = CommonOpticalFilter_h->PE_Veto_Total();

  // SCE service
  ///auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  //Save run, subrun, event in tree
  _event  = e.id().event();
  _run    = e.run();
  _subrun = e.subRun();

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    throw std::exception();
  }

  //  select largest flash within beam-spill window, if available
  size_t flash_ctr = 0;
  double PEmax = 0;
  size_t flash_idx = 0;
  for (size_t f=0; f < flash_h->size(); f++) {
    auto const& flash = flash_h->at(f);
    if ( (flash.Time() > fBeamSpillStart) && (flash.Time() < fBeamSpillEnd) ) {
      flash_ctr += 1;
      if (flash.TotalPE() > PEmax) {
	PEmax = flash.TotalPE();
	flash_idx = f;
      }// if largest flash so far
    }// if within beam spill window
  }// for all flashes

  // load waveform
  art::Handle<std::vector<raw::OpDetWaveform> > wf_h;
  e.getByLabel(fOpDetWfmProducer,wf_h);

  // make sure waveforms look good
  if(!wf_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate OpDetWaveform!"<<std::endl;
    throw std::exception();
  }

  _wf_v.resize(32);

  for (size_t w=0; w < wf_h->size(); w++) {
    
    auto const& wf = wf_h->at(w);
    auto ch = wf.ChannelNumber();
    if (ch >= 32) continue;
    _wf_v[ch] = std::vector<short>(800,0);
    if (wf.size() > 800) {
      for (size_t tick=0; tick < 800; tick++) 
	_wf_v[ch][tick] = wf[tick];
    }// if long enough waveform
  }// for all waveforms

  // load tracks previously created for which T0 reconstruction should occur
  art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(fTrackProducer,track_h);

  // make sure tracks look good
  if(!track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Track!"<<std::endl;
    throw std::exception();
  }


  // ADDITION FROM PETRILLO
  e.getValidHandle<std::vector<recob::PFParticle>>(fTrackProducer);
  // grab the hits associated to the tracks
  auto trk_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::PFParticle>::find(track_h, e, fTrackProducer);

  // grab calo associated to tracks
  art::FindMany<anab::Calorimetry> trk_calo_assn_v(track_h, e, fCaloProducer); 

  // if we are using CRT tagged tracks
  // create map linking track index to CRT hit time
  std::map<unsigned int, double> CRThitMap;
  if (fCRTTagProducer != "") {
    art::FindMany<anab::T0> trk_crttag_assn_v(track_h, e, fCRTTagProducer);
    for (size_t a=0; a < trk_crttag_assn_v.size(); a++) {
      auto trk_crttag_assn = trk_crttag_assn_v.at(a);
      if (trk_crttag_assn.size() == 1) 
	CRThitMap[ a ] = trk_crttag_assn.at(0)->Time();
    }// for all track -> crt t0 associations
  }// if we are loading CRT tags in the first place

  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, track_h);
  
  // loop through tracks
  size_t trkctr = 0;

  for (auto& track : TrkVec){

    trkctr += 1;

    if (track->Length() < fTrackLenMin) continue;
    
    auto const& beg = track->Vertex();
    auto const& end = track->End();


    // is this track associated to a CRT hit?
    if (CRThitMap.find(trkctr - 1) != CRThitMap.end()) {
      _crttag = 1;
      _crtt0 = CRThitMap[ trkctr - 1 ];
    }// found CRT hit
    else { _crttag = 0; }

    // filter out events that do not have a CRT tag if a CRT producer was specified  
    if ( (fCRTTagProducer != "") && (_crttag == 0) )
      continue;
    // if CRT tag out of time, skip
    if ( ( _crtt0 < (fBeamSpillStart - 5.) ) || ( _crtt0 > (fBeamSpillEnd + 5.) ) )
      continue;
    
    bool tagged = false;
    
    // cathode-side tracks
    // for tracks that enter from the cathode:
    if ( (beg.X() > end.X()) && (beg.Y() < fYMax) && (end.Y() < fYMin) ) {
      if ( (beg.X() > fCathodeMin) && (beg.X() < fCathodeMax) ) {
	tagged = true;
      }// if within cathode X constraint
    }// for tracks that enter from the side
    // for tracks that exit from the cathode:
    if ( (beg.X() < end.X()) && (beg.Y() > fYMax) && (end.Y() > fYMin) ) {
      if ( (end.X() > fCathodeMin) && (end.X() < fCathodeMax) ) {
	tagged = true;
      }// if within cathode X constraint
    }// for tracks that exit from the side

    // anode-side tracks
    // for tracks that enter from the side:
    if ( (beg.X() < end.X()) && (beg.Y() < fYMax) && (end.Y() < fYMin) ) {
      if ( (beg.X() > fAnodeMin) && (beg.X() < fAnodeMax) ) {
	tagged = true;
      }// if within anode X constraint
    }// for tracks that enter from the side
    // for tracks that exit from the side:
    if ( (beg.X() > end.X()) && (beg.Y() > fYMax) && (end.Y() > fYMin) ) {
      if ( (end.X() > fAnodeMin) && (end.X() < fAnodeMax) ) {
	tagged = true;
      }// if within anode X constraint
    }// for tracks that exit from the side
    
    // has the track been tagged?
    if (tagged == true) {

      // associated hits
      const std::vector< art::Ptr<recob::Hit> >& hit_v = trk_hit_assn_v.at(trkctr - 1);
      _trk_charge_0 = 0.;
      _trk_charge_1 = 0.;
      _trk_charge_2 = 0.;
      for (auto const& hit : hit_v) {
	  if (hit->WireID().Plane == 0 )
	    _trk_charge_0 += hit->Integral();
	  if (hit->WireID().Plane == 1 )
	    _trk_charge_1 += hit->Integral();
	  if (hit->WireID().Plane == 2 )
	    _trk_charge_2 += hit->Integral();
      }/// for all hits

      //Using Calorimetry information
      const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(trkctr - 1);
      
      _residualrange_vector_plane_0.clear();
      _dedx_vector_plane_0.clear();
      _residualrange_vector_plane_1.clear();
      _dedx_vector_plane_1.clear();
      _residualrange_vector_plane_2.clear();
      _dedx_vector_plane_2.clear();

      for (size_t pl=0; pl < Calo_v.size(); pl++)
	
	{
	  auto const& calo = Calo_v.at(pl);
	  auto const& plane = calo->PlaneID().Plane;
	  //	auto const& dqdx  = calo->dQdx();
	  auto const& dedx  = calo->dEdx();
	  auto const& rr    = calo->ResidualRange();
	  //auto const& xyz   = calo->XYZ();
	  //	fillCalorimetry(plane, dqdx, rr);//, xyz);
	  //_residualrange = rr;
	  //_dedx = dedx;
	  //auto const& Charge = calo->charge();
	  // auto const& dqdx = calo->dQdx();

	  //Will find the calorimetry energy for all three planes separately to see what we get depending on the plane (induction vs collection)
	  //Induction 1
	  if(plane==0){
	    Int_t N_plane_0 = 0;
	    N_plane_0 = rr.size();
	    _calorimetry_energy_plane_0 = 0.;
	    for(int i = 0; i < (N_plane_0 - 1); i++)
	      {
		_calorimetry_energy_plane_0 +=dedx.at(i)*(rr.at(i+1) - rr.at(i));
		_residualrange_vector_plane_0.push_back(rr.at(i));
		_dedx_vector_plane_0.push_back(dedx.at(i));
	      }
	  }
	  //induction 2
	  if(plane==1){
	    Int_t N_plane_1 = 0;
	    N_plane_1 = rr.size();
          _calorimetry_energy_plane_1 = 0.;
          for(int i = 0; i < (N_plane_1 - 1); i++)
            {
              _calorimetry_energy_plane_1 +=dedx.at(i)*(rr.at(i+1) - rr.at(i));
	      _residualrange_vector_plane_1.push_back(rr.at(i));
	      _dedx_vector_plane_1.push_back(dedx.at(i));
	    }
	  }
	  //collection
	  if(plane==2){
	    Int_t N = 0;
	    N = rr.size();
	    _calorimetry_energy_plane_2 = 0.;
	    for(int i = 0; i < (N-1); i++){
	      _calorimetry_energy_plane_2 += dedx.at(i)*(rr.at(i+1) - rr.at(i));
	      _residualrange_vector_plane_2.push_back(rr.at(i));
	      _dedx_vector_plane_2.push_back(dedx.at(i));
	    }
	  }
	  //	_residualrange_vector.push_back(rr);
	  //	_dedx_vector.push_back(calo->dEdx());
	  //_residualrange = rr.at(pl);
	  
	}// for all calorimetry objects associated to this track


      double time = 0.;
      // associate  flash if a flash exits
      if (flash_ctr > 0) {
	time = flash_h->at(flash_idx).Time();
	// get pointer to associated flash
	const art::Ptr<recob::OpFlash> flash_ptr(flash_h, flash_idx);
	trk_flash_assn_v->addSingle(track, flash_ptr);
      }// if there is an associated flash

      // create T0 object with this information!
      anab::T0 t0(time, 0, 0);
      T0_v->emplace_back(t0);
      util::CreateAssn(*this, e, *T0_v, track, *trk_t0_assn_v);

      // save TTree
      if (fSaveTree) {

	_len = track->Length();

	_trk_beg_x = beg.X();
	_trk_end_x = end.X();
	_trk_beg_y = beg.y();
	_trk_end_y = end.y();
	_trk_beg_z = beg.Z();
	_trk_end_z = end.Z();
	
	//Space Charge Effects Corrections implementation
	//	geo::Point_t posOffsets;
	geo::Point_t begOffsets;
	geo::Point_t endOffsets;
	auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
	//	if(sce->EnableCalSpatialSCE()){
	//  posOffsets = sce->GetCalPosOffsets(beg);
	//	}
	
	if(sce->EnableCalSpatialSCE()) {
	  begOffsets = sce->GetCalPosOffsets(beg);
	  _trk_beg_x_off = _trk_beg_x - begOffsets.X();
	  _trk_beg_y_off = _trk_beg_y + begOffsets.Y();
	  _trk_beg_z_off = _trk_beg_z + begOffsets.Z();
	}

	if(sce->EnableCalSpatialSCE()) {
	  endOffsets = sce->GetCalPosOffsets(end);
	  _trk_end_x_off = _trk_end_x - endOffsets.X();
	  _trk_end_y_off = _trk_end_y + endOffsets.Y();
	  _trk_end_z_off = _trk_end_z + endOffsets.Z();
	}
	
	_trk_len = sqrt((_trk_beg_x_off-_trk_end_x_off)*(_trk_beg_x_off-_trk_end_x_off)+(_trk_beg_y_off-_trk_end_y_off)*(_trk_beg_y_off-_trk_end_y_off)+(_trk_beg_z_off-_trk_end_z_off)*(_trk_beg_z_off-_trk_end_z_off));

	_nflash    = flash_ctr;
	_flash_t   = time;

	if (flash_ctr > 0) {
	  auto const& flash = flash_h->at(flash_idx);
	  _flash_pe   = flash.TotalPE();
	  _flash_pe_v = flash.PEs();
	  _flash_yc   = flash.YCenter();
	  _flash_yw   = flash.YWidth();
	  _flash_zc   = flash.ZCenter();
	  _flash_zw   = flash.ZWidth();
	}// if a flash was found
	_tree->Fill();
      }// if we are saving to the TTree

    }// if the track was tagged

  }// for all tracks

  bool selected = false;
  if (T0_v->size() > 0)
    selected = true;

  e.put(std::move(T0_v));
  e.put(std::move(trk_t0_assn_v));
  e.put(std::move(trk_flash_assn_v));

  return selected;
}

void ACPTtrig::ResetTTree() {

  _len = 0;
  _trk_beg_x = 0;
  _trk_beg_y = 0;
  _trk_beg_z = 0;
  _trk_end_x = 0;
  _trk_end_y = 0;
  _trk_end_z = 0;
  _trk_beg_x_off = 0;
  _trk_beg_y_off = 0;
  _trk_beg_z_off = 0;
  _trk_end_x_off = 0;
  _trk_end_y_off = 0;
  _trk_end_z_off = 0;
  _trk_len = 0;
  _flash_t = 0;
  _flash_pe = 0;
  _flash_yw = 0;
  _flash_yc = 0;
  _flash_zw = 0;
  _flash_zc = 0;
  _opfilter_pe_beam = 0;
  _opfilter_pe_beam_tot = 0;
  _opfilter_pe_veto = 0;
  _opfilter_pe_veto_tot = 0;
  _crttag = -1;
  _crtt0  = -1;
  _nflash = 0;
  _flash_pe_v.clear();
}

void ACPTtrig::beginJob()
{
  // Implementation of optional member function here.
}

void ACPTtrig::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ACPTtrig)
