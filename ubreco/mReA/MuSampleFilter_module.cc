////////////////////////////////////////////////////////////////////////
// Class:       MuSampleFilter
// Plugin Type: filter (art v3_01_02)
// File:        MuSampleFilter_module.cc
//
// Generated at Mon Feb 24 15:58:08 2020 by Matthew Rosenberg using cetskelgen
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

#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "Pandora/PdgTable.h"

#include <memory>

class MuSampleFilter;


class MuSampleFilter : public art::EDFilter {
public:
  explicit MuSampleFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuSampleFilter(MuSampleFilter const&) = delete;
  MuSampleFilter(MuSampleFilter&&) = delete;
  MuSampleFilter& operator=(MuSampleFilter const&) = delete;
  MuSampleFilter& operator=(MuSampleFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  std::string m_pfp_producer;
  std::string m_muon_producer;
  std::string m_pfpToVtx_producer;
  std::string m_pfpToTrk_producer;
  std::string m_mct_producer;

  bool m_applyTruthCuts;
  bool m_applyStrictSelection;

  float m_fidVolXStart;
  float m_fidVolXEnd;
  float m_fidVolYStart;
  float m_fidVolYEnd;
  float m_fidVolZStart;
  float m_fidVolZEnd;
  float m_cntBuffer; //min distance from track end to detector edge for track to be considered contained

  float m_nuScoreCut;
  float m_muTrackScoreCut;

  bool PrintAndReturn(bool evtPassed);
  bool isFiducial(const double x[3]) const;
  bool isContained(const double x[3]) const;

  lar_pandora::LArPandoraHelper larpandora;

};


MuSampleFilter::MuSampleFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
  m_pfp_producer(p.get<std::string>("PFParticleLabel","pandora")),
  m_muon_producer(p.get<std::string>("MuonFilterLabel","NuCCproducer")),
  m_pfpToVtx_producer(p.get<std::string>("PFPartToVertexLabel","pandora")),
  m_pfpToTrk_producer(p.get<std::string>("PFPartToTrackLabel","pandora")),
  m_mct_producer(p.get<std::string>("MCTruthLabel","generator")),
  m_applyTruthCuts(p.get<bool>("ApplyTruthCuts",false)),
  m_applyStrictSelection(p.get<bool>("ApplyStrictSelection",true)),
  m_fidVolXStart(p.get<float>("FidVolXStart",20.)),
  m_fidVolXEnd(p.get<float>("FidVolXEnd",20.)),
  m_fidVolYStart(p.get<float>("FidVolYStart",20.)),
  m_fidVolYEnd(p.get<float>("FidVolYEnd",20.)),
  m_fidVolZStart(p.get<float>("FidVolZStart",20.)),
  m_fidVolZEnd(p.get<float>("FidVolZEnd",60.)),
  m_cntBuffer(p.get<float>("ContainmentBuffer", 10.)),
  m_nuScoreCut(p.get<float>("NuScoreCut", 0.2)),
  m_muTrackScoreCut(p.get<float>("MuTrackScoreCut", 0.98))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


bool MuSampleFilter::filter(art::Event& e)
{
  // Implementation of required member function here.

  //Apply truth cuts if requested
  if(m_applyTruthCuts){

    //Verify we have vaild MCTruth object
    art::Handle< std::vector<simb::MCTruth> > mcTruthBlocks;
    e.getByLabel(m_mct_producer, mcTruthBlocks);
    if(!mcTruthBlocks.isValid()) return PrintAndReturn(false);
    if(mcTruthBlocks->size() != 1) return PrintAndReturn(false);
    const art::Ptr<simb::MCTruth> mcTruth(mcTruthBlocks, 0);
    if(mcTruth->Origin() != simb::kBeamNeutrino) return PrintAndReturn(false);

    //Only look at true nu_mu CC interaction in fiducial volume
    const simb::MCNeutrino& mcNu = mcTruth->GetNeutrino();
    const int neutrinoPDG = mcNu.Nu().PdgCode();
    const int leptonPDG = mcNu.Lepton().PdgCode();
    if(std::abs(neutrinoPDG) != pandora::NU_MU) return PrintAndReturn(false);
    if(std::abs(leptonPDG) != pandora::MU_MINUS) return PrintAndReturn(false);
    const double mcNuVtx[3] = {mcNu.Nu().Vx(), mcNu.Nu().Vy(), mcNu.Nu().Vz()};
    if(!isFiducial(mcNuVtx)) return PrintAndReturn(false);

  }

  //Find reconstructed muon candidate
  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  e.getByLabel(m_pfp_producer, pfParticleHandle);
  art::FindManyP<anab::T0> pfp_muon_assn(pfParticleHandle, e, m_muon_producer);
  std::map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;
  std::vector<size_t> muonParticleKey;

  for(unsigned int iPF = 0; iPF < pfParticleHandle->size(); ++iPF){
    const art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, iPF);
    if(!pfParticleIdMap.emplace(pfParticle->Self(), pfParticle).second)
      throw cet::exception("MuSampleFilter") << "Repeated PFParticles in input!" << std::endl;
    const std::vector< art::Ptr<anab::T0> > T0_muon = pfp_muon_assn.at(pfParticle.key());
    if(!T0_muon.empty()) muonParticleKey.push_back(pfParticle->Self());
  }

  if(muonParticleKey.empty()) throw cet::exception("MuSampleFilter") <<
    "No tagged muon candidate. Run only over filtered output (NuCCfilter or similar)" << std::endl;

  if(muonParticleKey.size() > 1) throw cet::exception("MuSampleFilter") <<
    "Multiple tagged muon candidates. Run only over filtered output (NuCCfilter or similar)" << std::endl;

  const art::Ptr<recob::PFParticle>& recoMuon = pfParticleIdMap[ muonParticleKey[0] ];

  //Get reconstructed neutrino vertex
  const auto nuIterator = pfParticleIdMap.find(recoMuon->Parent());
  if(nuIterator == pfParticleIdMap.end()) throw cet::exception("MuSampleFilter")
    << "reco muon's parent is not in the input collection (should be neutrino)" << std::endl;
  const auto recoNu = nuIterator->second;
  if(!recoNu->IsPrimary() || std::abs(recoNu->PdgCode()) != pandora::NU_MU)
    throw cet::exception("MuSampleFilter") << "reco muon's parent is not a primary nu_mu" << std::endl;

  art::FindManyP<recob::Vertex> pfPartToVertexAssoc(pfParticleHandle, e, m_pfpToVtx_producer);
  const auto associatedVertices = pfPartToVertexAssoc.at(recoNu.key());
  if(associatedVertices.size() != 1) throw cet::exception("MuSampleFilter") <<
    "Number of neutrino vertices != 1. Is this supposed to happen?" << std::endl;
  recob::tracking::Point_t nuVtxPos = associatedVertices[0]->position();

  //apply stricter topological and track score cuts if requested
  if(m_applyStrictSelection){
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticlesToMetadata particlesToMetadata;
    larpandora.CollectPFParticleMetadata(e, m_pfp_producer, pfparticles, particlesToMetadata);
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties =
     particlesToMetadata.at(recoNu).front()->GetPropertiesMap();
    const larpandoraobj::PFParticleMetadata::PropertiesMap &muon_properties = 
     particlesToMetadata.at(recoMuon).front()->GetPropertiesMap();
    float nuScore = neutrino_properties.at("NuScore");
    float muTrackScore = muon_properties.at("TrackScore");
    if(nuScore < m_nuScoreCut || muTrackScore < m_muTrackScoreCut) return PrintAndReturn(false);
  }

  //apply fiducial volume cut
  const double nuVtx[3] = {nuVtxPos.X(), nuVtxPos.Y(), nuVtxPos.Z()};
  if(!isFiducial(nuVtx)) return PrintAndReturn(false);

  //Get reconstructed muon track end position
  art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, e, m_pfpToTrk_producer);
  const auto associatedTracks = pfPartToTrackAssoc.at(recoMuon.key());
  if(associatedTracks.size() != 1) throw cet::exception("MuSampleFilter")
    << "reconstructed muon does not have a single associated track" << std::endl;
  recob::tracking::Point_t muEndPos = associatedTracks[0]->End();

  //apply muon containment cut
  const double muEnd[3] = {muEndPos.X(), muEndPos.Y(), muEndPos.Z()};
  return PrintAndReturn( isContained(muEnd) );

}


bool MuSampleFilter::PrintAndReturn(bool evtPassed){
  if(evtPassed) mf::LogDebug("MuSampleFilter") << "event passed\n";
  else mf::LogDebug("MuSampleFilter") << "event failed\n";
  return evtPassed;
}


//copied from NuMuSelection::isFiducial from elee analysis
bool MuSampleFilter::isFiducial(const double x[3]) const
{      
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2. * geo->DetHalfWidth(),
    -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};
        
  bool is_x = x[0] > (bnd[0] + m_fidVolXStart) && x[0] < (bnd[1] - m_fidVolXEnd);
  bool is_y = x[1] > (bnd[2] + m_fidVolYStart) && x[1] < (bnd[3] - m_fidVolYEnd);
  bool is_z = x[2] > (bnd[4] + m_fidVolZStart) && x[2] < (bnd[5] - m_fidVolZEnd);
        
  return is_x && is_y && is_z;
}


bool MuSampleFilter::isContained(const double x[3]) const
{      
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2. * geo->DetHalfWidth(),
    -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};
        
  bool is_x = x[0] > (bnd[0] + m_cntBuffer) && x[0] < (bnd[1] - m_cntBuffer);
  bool is_y = x[1] > (bnd[2] + m_cntBuffer) && x[1] < (bnd[3] - m_cntBuffer);
  bool is_z = x[2] > (bnd[4] + m_cntBuffer) && x[2] < (bnd[5] - m_cntBuffer);
        
  return is_x && is_y && is_z;
}


DEFINE_ART_MODULE(MuSampleFilter)
