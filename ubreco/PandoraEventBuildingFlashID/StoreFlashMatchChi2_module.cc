////////////////////////////////////////////////////////////////////////
// Class:       StoreFlashMatchChi2
// Plugin Type: producer (art v3_01_02)
// File:        StoreFlashMatchChi2_module.cc
//
// Generated at Mon May 20 07:11:36 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"

//#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderFMWKInterface.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

#include "larcore/Geometry/Geometry.h"

#include "TTree.h"

#include <memory>

class StoreFlashMatchChi2;


class StoreFlashMatchChi2 : public art::EDProducer {
public:
  explicit StoreFlashMatchChi2(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoreFlashMatchChi2(StoreFlashMatchChi2 const&) = delete;
  StoreFlashMatchChi2(StoreFlashMatchChi2&&) = delete;
  StoreFlashMatchChi2& operator=(StoreFlashMatchChi2 const&) = delete;
  StoreFlashMatchChi2& operator=(StoreFlashMatchChi2&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  ::flashana::FlashMatchManager m_flashMatchManager; ///< The flash match manager
  art::InputTag fFlashProducer;
  art::InputTag fT0Producer; // producer for ACPT in-time anab::T0 <-> recob::Track assocaition
  std::string fPandoraProducer, fSpacePointProducer;
  float fBeamWindowEnd, fBeamWindowStart;
  float fMaxTotalPE;
  float fChargeToNPhotonsShower, fChargeToNPhotonsTrack;
  std::vector<float> fPMTChannelCorrection;

  ::flashana::Flash_t GetFlashPESpectrum(const recob::OpFlash& opflash);

  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap, 
				    const art::Ptr<recob::PFParticle> &particle,
				    lar_pandora::PFParticleVector &downstreamPFParticles) const;

  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap, 
				    const lar_pandora::PFParticleVector &parentPFParticles,
				    lar_pandora::PFParticleVector &downstreamPFParticles) const;

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  
		    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, 
		    std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

  TTree* _flashmatch_acpt_tree;
  TTree* _flashmatch_nuslice_tree;
  std::vector<float> _pe_reco_v, _pe_hypo_v;
  float _trk_vtx_x, _trk_vtx_y, _trk_vtx_z, _trk_end_x, _trk_end_y, _trk_end_z;
  float _nuvtx_x, _nuvtx_y, _nuvtx_z;
  int _evt, _run, _sub;
  float _flashtime;
  float _flashpe;

  // PFP map
  std::map<unsigned int, unsigned int> _pfpmap;

};


StoreFlashMatchChi2::StoreFlashMatchChi2(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{

  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::PFParticle, anab::T0> >();

  fFlashProducer      = p.get<art::InputTag>("FlashProducer"   );
  fPandoraProducer    = p.get<std::string>("PandoraProducer"   );
  fT0Producer         = p.get<std::string>("T0Producer"        );
  fSpacePointProducer = p.get<std::string>("SpacePointProducer");
  fBeamWindowStart = p.get<float>("BeamWindowStart");
  fBeamWindowEnd   = p.get<float>("BeamWindowEnd");
  fMaxTotalPE      = p.get<float>("MaxTotalPE");
  fChargeToNPhotonsShower   = p.get<float>("ChargeToNPhotonsShower");
  fChargeToNPhotonsTrack    = p.get<float>("ChargeToNPhotonsTrack");
  fPMTChannelCorrection         = p.get<std::vector<float>>("PMTChannelCorrection");

  m_flashMatchManager.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  art::ServiceHandle<art::TFileService> tfs;

  // Tree to store ACPT track flash-matching information
  _flashmatch_acpt_tree = tfs->make<TTree>("ACPTFMtree","ACPT FlashMatch tree");
  _flashmatch_acpt_tree->Branch("evt",&_evt,"evt/I");
  _flashmatch_acpt_tree->Branch("run",&_run,"run/I");
  _flashmatch_acpt_tree->Branch("sub",&_sub,"sub/I");
  _flashmatch_acpt_tree->Branch("flashtime",&_flashtime,"flashtime/F");
  _flashmatch_acpt_tree->Branch("flashpe"  ,&_flashpe  ,"flashpe/F"  );
  _flashmatch_acpt_tree->Branch("pe_reco_v","std::vector<float>",&_pe_reco_v);
  _flashmatch_acpt_tree->Branch("pe_hypo_v","std::vector<float>",&_pe_hypo_v);
  _flashmatch_acpt_tree->Branch("trk_beg_x",&_trk_vtx_x,"trk_beg_x/F");
  _flashmatch_acpt_tree->Branch("trk_beg_y",&_trk_vtx_y,"trk_beg_y/F");
  _flashmatch_acpt_tree->Branch("trk_beg_z",&_trk_vtx_z,"trk_beg_z/F");
  _flashmatch_acpt_tree->Branch("trk_end_x",&_trk_end_x,"trk_end_x/F");
  _flashmatch_acpt_tree->Branch("trk_end_y",&_trk_end_y,"trk_end_y/F");
  _flashmatch_acpt_tree->Branch("trk_end_z",&_trk_end_z,"trk_end_z/F");


  // Tree to store neutrino flash-matching information
  _flashmatch_nuslice_tree = tfs->make<TTree>("nuslicetree","ACPT FlashMatch tree");
  _flashmatch_nuslice_tree->Branch("evt",&_evt,"evt/I");
  _flashmatch_nuslice_tree->Branch("run",&_run,"run/I");
  _flashmatch_nuslice_tree->Branch("sub",&_sub,"sub/I");
  _flashmatch_nuslice_tree->Branch("flashtime",&_flashtime,"flashtime/F");
  _flashmatch_nuslice_tree->Branch("flashpe"  ,&_flashpe  ,"flashpe/F"  );
  _flashmatch_nuslice_tree->Branch("pe_reco_v","std::vector<float>",&_pe_reco_v);
  _flashmatch_nuslice_tree->Branch("pe_hypo_v","std::vector<float>",&_pe_hypo_v);
  _flashmatch_nuslice_tree->Branch("nuvtx_x",&_nuvtx_x,"nuvtx_x/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_y",&_nuvtx_y,"nuvtx_y/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_z",&_nuvtx_z,"nuvtx_z/F");
  

  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void StoreFlashMatchChi2::produce(art::Event& e)
{

  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::PFParticle, anab::T0> > pfp_t0_assn_v( new art::Assns<recob::PFParticle, anab::T0>  );

  // reset TTree variables
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  _flashtime = -9999.;
  _flashpe   = -9999.;

  //  prepare flash object
  ::flashana::Flash_t beamflash;
  beamflash.time = 0.;
  float maxEventPE = 0.;

  // Collect all flashes from the event
  const auto flashes(*e.getValidHandle< std::vector<recob::OpFlash> >(fFlashProducer));


  for (const auto &flash : flashes) {

    // get claibrated PE spectrum
    ::flashana::Flash_t thisflash = GetFlashPESpectrum(flash);

    if ( (thisflash.time < fBeamWindowStart) || (thisflash.time > fBeamWindowEnd) )
      continue;
    
    auto totalPE = std::accumulate(thisflash.pe_v.begin(), thisflash.pe_v.end(), 0);

    if (totalPE < fMaxTotalPE)
      continue;

    if (totalPE < maxEventPE) 
      continue;

    // made it this far, save as beam candidate    
    maxEventPE = totalPE;
    beamflash = thisflash;
  }// for all flashes


  if (beamflash.time == 0) {
    // quit here!
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  _flashtime = beamflash.time;
  _flashpe   = std::accumulate( beamflash.pe_v.begin(), beamflash.pe_v.end(), 0);

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);
  
  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPandoraProducer);
  // grab tracks associated with PFParticles
  art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, fPandoraProducer);
  // grab vertices associated with PFParticles
  art::FindManyP<recob::Vertex> pfp_vertex_assn_v(pfp_h, e, fPandoraProducer);

  // grab tracks in the event
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track> >(fPandoraProducer);
  // grab anab::T0 for ACPT in-time associated to Tracks
  art::FindManyP<anab::T0> trk_t0_assn_v(trk_h, e, fT0Producer);
  
  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPandoraProducer);    
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointProducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointProducer);

  _pfpmap.clear();

  for (unsigned int p=0; p < pfp_h->size(); p++)
    _pfpmap[pfp_h->at(p).Self()] = p;
  
  for (unsigned int p=0; p < pfp_h->size(); p++){
    
    auto const& pfp = pfp_h->at(p);

    ::flashana::Flash_t beamflashcopy = beamflash;
    m_flashMatchManager.Reset();
    flashana::QCluster_t lightCluster;
    
    // start from primary PFParticles
    if (pfp.IsPrimary() == false) continue;
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    
    // now build vectors of PFParticles, space-points, and hits for this slice
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint>>> spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit>>> hit_v_v;
    
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

    // go through these pfparticles and fill info needed for matching
    for (size_t i=0; i < pfp_ptr_v.size(); i++) {
      
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      
      pfp_v.push_back(pfp);
      
      //auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      
      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {

	auto SP = spacepoint_ptr_v[sp];

	auto const& spkey = SP.key();
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	for (size_t h=0; h < this_hit_ptr_v.size(); h++) {

	  auto hit = this_hit_ptr_v.at( h );

	  // Only use hits from the collection plane
	  if (hit->View() != geo::kZ)
	    continue;
	  
	  // Add the charged point to the vector
	  const auto &position(SP->XYZ());
	  const auto charge(hit->Integral());
	  lightCluster.emplace_back(position[0], position[1], position[2], charge * (lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr) ? fChargeToNPhotonsTrack : fChargeToNPhotonsShower));

	}// for all hits associated to this spacepoint
      }// fpr all spacepoints
      
    }// for all pfp pointers
    
    // Perform the match
    m_flashMatchManager.Emplace(std::move(beamflashcopy));
    m_flashMatchManager.Emplace(std::move(lightCluster));


    const auto matches(m_flashMatchManager.Match());
    
    float FMscore = 9999.;
    
    _pe_hypo_v.clear();
    if (matches.size() != 0) {
      FMscore = matches.back().score;
      for (auto hypo_pe : matches.back().hypothesis)
        _pe_hypo_v.push_back(static_cast<float>(hypo_pe));
    }

    // is this the neutrino slice?
    if ( (pfp.PdgCode() == 12) || (pfp.PdgCode() == 14) ) {
      // grab vertex
      auto const& vtx_assn_v = pfp_vertex_assn_v.at(p);
      if (vtx_assn_v.size() == 1) {
	art::Ptr<recob::Vertex> vtx = vtx_assn_v.at(0);
	Double_t xyz[3] = {};
	vtx->XYZ(xyz);
	_nuvtx_x = xyz[0];
	_nuvtx_y = xyz[1];
	_nuvtx_z = xyz[2];
	_flashmatch_nuslice_tree->Fill();
      }// if vertex is found
    }// if neutrino reco slice
    
    // load associated tracks
    auto const& pfp_track_assn = pfp_track_assn_v.at(p);
    if ( pfp_track_assn.size() == 1 ) {
      art::Ptr<recob::Track> trk = pfp_track_assn.at(0);
      // grab track key to find anab::T0 association, if present
      auto const T0_v = trk_t0_assn_v.at( trk.key() );
      if (T0_v.size() == 1) {
	_trk_vtx_x = trk->Vertex().X();
	_trk_vtx_y = trk->Vertex().Y();
	_trk_vtx_z = trk->Vertex().Z();
	_trk_end_x = trk->End().X();
	_trk_end_y = trk->End().Y();
	_trk_end_z = trk->End().Z();
	_flashmatch_acpt_tree->Fill();
      }// if this track is ACPT tagged
    }// if there is a track associated to this primary pfparticle

    // create T0 object with this information!
    anab::T0 t0(beamflash.time, 0, 0, 0, FMscore);
    T0_v->emplace_back(t0);
    util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
    
  }//  for all PFParticles

  e.put(std::move(T0_v));
  e.put(std::move(pfp_t0_assn_v));

}

::flashana::Flash_t StoreFlashMatchChi2::GetFlashPESpectrum(const recob::OpFlash& opflash) {
  
  // prepare conainer to store flash
  ::flashana::Flash_t flash;
  flash.time = opflash.Time();

  
  // geometry service
  const art::ServiceHandle<geo::Geometry> geometry;
  uint nOpDets(geometry->NOpDets());
  std::vector<float> PEspectrum;
  PEspectrum.resize(nOpDets);
  
  // apply gain to OpDets
  for (uint OpChannel = 0; OpChannel < nOpDets; ++OpChannel)
    {
      uint opdet = geometry->OpDetFromOpChannel(OpChannel);
      PEspectrum[opdet] = opflash.PEs().at(OpChannel);
    }

  _pe_reco_v = PEspectrum;
  
  // Reset variables
  flash.x = flash.y = flash.z = 0;
  flash.x_err = flash.y_err = flash.z_err = 0;
  float totalPE = 0.;
  float sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
  
  for (unsigned int opdet = 0; opdet < PEspectrum.size(); opdet++)
    {
      double PMTxyz[3];
      geometry->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
      // Add up the position, weighting with PEs
      sumy += PEspectrum[opdet] * PMTxyz[1];
      sumy2 += PEspectrum[opdet] * PMTxyz[1] * PMTxyz[1];
      sumz += PEspectrum[opdet] * PMTxyz[2];
      sumz2 += PEspectrum[opdet] * PMTxyz[2] * PMTxyz[2];
      totalPE += PEspectrum[opdet];
    }
  flash.y = sumy / totalPE;
  flash.z = sumz / totalPE;
  // This is just sqrt(<x^2> - <x>^2)
  if ((sumy2 * totalPE - sumy * sumy) > 0.)
    flash.y_err = std::sqrt(sumy2 * totalPE - sumy * sumy) / totalPE;
  
  if ((sumz2 * totalPE - sumz * sumz) > 0.)
    flash.z_err = std::sqrt(sumz2 * totalPE - sumz * sumz) / totalPE;
  
  // Set the flash properties
  flash.pe_v.resize(nOpDets);
  flash.pe_err_v.resize(nOpDets);
  
  // Fill the flash with the PE spectrum
  for (unsigned int i = 0; i < nOpDets; ++i)
    {
      const auto PE(PEspectrum.at(i));
      flash.pe_v.at(i) = PE;
      flash.pe_err_v.at(i) = std::sqrt(PE);
    }
  
  if (flash.pe_v.size() != nOpDets)
    throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;

  return flash;
}


void StoreFlashMatchChi2::CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap, 
						       const art::Ptr<recob::PFParticle> &particle,
						       lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
    downstreamPFParticles.push_back(particle);
  
  for (const auto &daughterId : particle->Daughters())
    {
      const auto iter(pfParticleMap.find(daughterId));
      if (iter == pfParticleMap.end())
	throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;
      
      this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
    }
}


void StoreFlashMatchChi2::CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap, 
						       const lar_pandora::PFParticleVector &parentPFParticles,
						       lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  for (const auto &particle : parentPFParticles)
    this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
}

void StoreFlashMatchChi2::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
  auto daughters = pfp_ptr->Daughters();
  
  pfp_v.push_back(pfp_ptr);
  
  for(auto const& daughterid : daughters) {
    
    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
    AddDaughters(pfp_ptr, pfp_h, pfp_v);
    
  }// for all daughters
  
  return;
}


void StoreFlashMatchChi2::beginJob()
{
  // Implementation of optional member function here.
}

void StoreFlashMatchChi2::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(StoreFlashMatchChi2)
