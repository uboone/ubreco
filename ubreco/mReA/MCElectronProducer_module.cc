////////////////////////////////////////////////////////////////////////
// Class:       MCElectronProducer
// Plugin Type: producer (art v3_01_02)
// File:        MCElectronProducer_module.cc
//
// Generated at Thu Feb 13 11:53:09 2020 by Matthew Rosenberg using cetskelgen
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

#include "canvas/Persistency/Common/FindManyP.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "Pandora/PdgTable.h"

#include "TLorentzVector.h"

#include <memory>
#include <limits>

class MCElectronProducer;


class MCElectronProducer : public art::EDProducer {
public:
  explicit MCElectronProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCElectronProducer(MCElectronProducer const&) = delete;
  MCElectronProducer(MCElectronProducer&&) = delete;
  MCElectronProducer& operator=(MCElectronProducer const&) = delete;
  MCElectronProducer& operator=(MCElectronProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Optional functions.
  void endJob() override;

private:

  // Declare member data here.
  std::string m_pfp_producer;
  std::string m_muon_producer;
  std::string m_trk_producer;
  std::string m_pfpToVtx_producer;
  std::string m_pfpToTrk_producer;
  std::string m_flash_producer;
  std::string m_trkToOpf_producer;
  std::string m_mcp_producer;

  bool m_truthOverride;
  bool m_runMuonTest;

  std::string m_beamName;
  std::string m_triggerName;
  double m_beamLow;
  double m_beamHigh;
  double m_triggerTime;

  unsigned int m_noFlash_counter;
  unsigned int m_adjFlash_counter;

  trkf::TrackMomentumCalculator tmc;
  pandora::PdgTable partInfo;
  double m_muMass;
  double m_eMass;

  //Set Beam gate in ns, geant coordinates (based on table 2 in docdb 12290)
  void SetBeamGate(const std::string& beamName);
  //Set trigger time in ns, geant coordinates (based on table 3 in docdb 12290)
  void SetTriggerTime(const std::string& triggerName);
  //Find earliest flash time in collection and return value in ns, geant coordinates
  //double GetFirstFlashTime(const std::vector< art::Ptr<recob::OpFlash> >& flashVec);
  double FlashTime(const art::Ptr<recob::OpFlash>& flash){ return flash->Time(); }
  double FlashTime(const recob::OpFlash& flash){ return flash.Time(); }
  template <typename T>
  double GetFirstFlashTime(const std::vector<T>& flashVec);

};


void MCElectronProducer::SetBeamGate(const std::string& beamName){
  if(beamName == "BNB"){
    m_beamLow = 3125.0;
    m_beamHigh = 4725.0;
  }
  else if(beamName == "NUMI"){
    m_beamLow = 4687.5;
    m_beamHigh = 14287.5;
  }
  else{
    throw cet::exception("MCElectronProducer") << "Unknown beam name (" <<
     beamName << ") supplied. Use either BNB or NUMI" << std::endl;
  }
}


void MCElectronProducer::SetTriggerTime(const std::string& triggerName){
  if(triggerName == "none") m_triggerTime = 0.;
  else if(triggerName == "BNB") m_triggerTime = 31.25;
  else if(triggerName == "NUMI") m_triggerTime = 54.8675;
  else if(triggerName == "EXT") m_triggerTime = 414.0625;
  else if(triggerName == "MUCS") m_triggerTime = 405.25;
  else throw cet::exception("MCElectronProducer") << "Unknown trigger name (" <<
   triggerName << ") supplied. Use either none, BNB, NUMI, EXT, or MUCS" << std::endl;
}


MCElectronProducer::MCElectronProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  m_pfp_producer(p.get<std::string>("PFParticleLabel","pandora")),
  m_muon_producer(p.get<std::string>("MuonFilterLabel","NuCCproducer")),
  m_trk_producer(p.get<std::string>("TrackLabel","pandora")),
  m_pfpToVtx_producer(p.get<std::string>("PFPartToVertexLabel","pandora")),
  m_pfpToTrk_producer(p.get<std::string>("PFPartToTrackLabel","pandora")),
  m_flash_producer(p.get<std::string>("OpticalFlashLabel", "simpleFlashBeam")),
  m_trkToOpf_producer(p.get<std::string>("TrackToFlashLabel","acpttrigtagger")),
  m_mcp_producer(p.get<std::string>("MCParticleLabel","largeant")),
  m_truthOverride(p.get<bool>("TruthOverride", false)),
  m_runMuonTest(p.get<bool>("RunMuonTest", false)),
  m_beamName(p.get<std::string>("BeamName","BNB")),
  m_triggerName(p.get<std::string>("TriggerName","BNB")),
  m_noFlash_counter(0),
  m_adjFlash_counter(0)
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<simb::MCTruth> >();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  SetBeamGate(m_beamName);
  SetTriggerTime(m_triggerName);
  m_muMass = partInfo.GetParticleMass(pandora::MU_MINUS);
  m_eMass = partInfo.GetParticleMass(pandora::E_MINUS);
}


//double MCElectronProducer::GetFirstFlashTime(const std::vector< art::Ptr<recob::OpFlash> >& flashVec){
template <typename T>
double MCElectronProducer::GetFirstFlashTime(const std::vector<T>& flashVec){
  //default to middle of beam gate
  if(flashVec.size() < 1){
    ++m_noFlash_counter;
    mf::LogWarning("MCElectronProducer") << "No optical flashes found in " << m_flash_producer
     << " collection. Setting MC electron time to the middle of the beam gate..." << std::endl;
    return (m_beamLow + m_beamHigh)/2.0;
  }
  //find earliest flash time
  double simTime = std::numeric_limits<double>::max();
  for(const auto flash : flashVec){
    if(FlashTime(flash) < simTime) simTime = FlashTime(flash);
  }
  //convert to geant coords (us -> ns, subtract trigger time and 60ns PMT readout delay)
  simTime = simTime*1000. - m_triggerTime - 60.;
  //check if flash time is in the beam window
  if(simTime < m_beamLow){
    ++m_adjFlash_counter;
    mf::LogWarning("MCElectronProducer") << "Optical flash used for nu interaction is outside " <<
     "the beam window (t = " << simTime <<"ns, geant coords). Setting MC electron time to " <<
     "lower bound of beam gate (" << m_beamLow << " ns)" << std::endl;
    simTime = m_beamLow;
  }
  if(simTime > m_beamHigh){
    ++m_adjFlash_counter;
    mf::LogWarning("MCElectronProducer") << "Optical flash used for nu interaction is outside " <<
     "the beam window (t = " << simTime <<"ns, geant coords). Setting MC electron time to " <<
     "upper bound of beam gate (" << m_beamHigh << " ns)" << std::endl;
    simTime = m_beamHigh;
  }
  return simTime;
}


void MCElectronProducer::produce(art::Event& e)
{
  // Implementation of required member function here.
  mf::LogDebug("MCElectronProducer") << "Creating an MCTruth object containing a single electron MCParticle"
   << " with the kinematics of the removed muon\n";

  //Create MCParticle for simulated electron
  int simEtrackId = 4;
  int simEPDG = pandora::E_MINUS;
  double simEmass = m_eMass;
  if(m_runMuonTest){
    simEPDG = pandora::MU_MINUS;
    simEmass = m_muMass;
  }
  simb::MCParticle simElec(simEtrackId, simEPDG, "primary", -1, simEmass, 1);

  //Use truth info if truth override was requested
  if(m_truthOverride){

    //Get simulated MC muon
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    e.getByLabel(m_mcp_producer, mcParticleHandle);

    for (unsigned int i = 0; i < mcParticleHandle->size(); ++i){

      const art::Ptr<simb::MCParticle> mcMuon(mcParticleHandle, i);

      if(std::abs(mcMuon->PdgCode()) == pandora::MU_MINUS && mcMuon->Mother() == 0){

        //MC muon found. Calculate energy and assign start position, 4-momentum
        double leptonE = std::sqrt( std::pow(mcMuon->Px(),2.) + std::pow(mcMuon->Py(),2.) +
         std::pow(mcMuon->Pz(),2.) + std::pow(m_eMass,2.) );
        if(m_runMuonTest) leptonE = mcMuon->E();

        TLorentzVector vtx(mcMuon->Vx(), mcMuon->Vy(), mcMuon->Vz(), mcMuon->T());
        TLorentzVector mom(mcMuon->Px(), mcMuon->Py(), mcMuon->Pz(), leptonE);
        simElec.AddTrajectoryPoint(vtx, mom);

        //Create MCTruth object for simulated electron, write it to art event
        simb::MCTruth truth;
        truth.Add(simElec);
        auto outputTruth = std::make_unique< std::vector<simb::MCTruth> >();
        outputTruth->emplace_back(std::move(truth));
        e.put(std::move(outputTruth));

        return;
      }
    }

    throw cet::exception("MCElectronProducer") << "Couldn't find primary true muon. "
     << "Input should be MC CC nu_mu events if the TruthOverride option is set" << std::endl;
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
      throw cet::exception("MCElectronProducer") << "Repeated PFParticles in input!" << std::endl;
    const std::vector< art::Ptr<anab::T0> > T0_muon = pfp_muon_assn.at(pfParticle.key());
    if(!T0_muon.empty()) muonParticleKey.push_back(pfParticle->Self());
  }

  if(muonParticleKey.empty()) throw cet::exception("MCElectronProducer") <<
    "No tagged muon candidate. Run only over filtered output (NuCCfilter or similar)" << std::endl;

  if(muonParticleKey.size() > 1) throw cet::exception("MCElectronProducer") <<
    "Multiple tagged muon candidates. Run only over filtered output (NuCCfilter or similar)" << std::endl;

  const art::Ptr<recob::PFParticle>& recoMuon = pfParticleIdMap[ muonParticleKey[0] ];

  //Get reconstructed neutrino vertex
  const auto nuIterator = pfParticleIdMap.find(recoMuon->Parent());
  if(nuIterator == pfParticleIdMap.end()) throw cet::exception("MCElectronProducer")
    << "reco muon's parent is not in the input collection (should be neutrino)" << std::endl;
  const auto recoNu = nuIterator->second;
  if(!recoNu->IsPrimary() || std::abs(recoNu->PdgCode()) != pandora::NU_MU)
    throw cet::exception("MCElectronProducer") << "reco muon's parent is not a primary nu_mu" << std::endl;

  art::FindManyP<recob::Vertex> pfPartToVertexAssoc(pfParticleHandle, e, m_pfpToVtx_producer);
  const auto associatedVertices = pfPartToVertexAssoc.at(recoNu.key());
  if(associatedVertices.size() != 1) throw cet::exception("MCElectronProducer") <<
    "Number of neutrino vertices != 1. Is this supposed to happen?" << std::endl;
  recob::tracking::Point_t nuVtxPos = associatedVertices[0]->position();

  //Get reconstructed muon momentum
  art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, e, m_pfpToTrk_producer);
  const auto associatedMuTracks = pfPartToTrackAssoc.at(recoMuon.key());
  if(associatedMuTracks.size() != 1) throw cet::exception("MCElectronProducer")
    << "reconstructed muon does not have a single associated track" << std::endl;
  double recoMuMomentum = tmc.GetTrackMomentum(associatedMuTracks[0]->Length(), pandora::MU_MINUS);
  double recoMuE = std::sqrt( std::pow(recoMuMomentum,2.) + std::pow(m_muMass,2.) );
  double electronE = std::sqrt( std::pow(recoMuMomentum,2.) + std::pow(m_eMass,2.) );
  recob::tracking::Vector_t recoMuDir = associatedMuTracks[0]->StartMomentumVector();
  double recoMuPx = recoMuMomentum*recoMuDir.X();
  double recoMuPy = recoMuMomentum*recoMuDir.Y();
  double recoMuPz = recoMuMomentum*recoMuDir.Z();

  //find neutrino interaction time from optical flash data (default to center of beam gate)
  double nuInt_time = (m_beamLow + m_beamHigh)/2.0;
  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(m_trk_producer, trackHandle);
  art::FindManyP<recob::OpFlash> track_opf_assn(trackHandle, e, m_trkToOpf_producer);
  const auto muFlashes = track_opf_assn.at(associatedMuTracks[0].key());
  //Use the optical flash associated with the muon track if this association exists
  if(!muFlashes.empty()){
    nuInt_time = GetFirstFlashTime(muFlashes);
  }
  //otherwise, use earliest flash time in the collection
  else{
    art::Handle< std::vector<recob::OpFlash> > flashHandle;
    e.getByLabel(m_flash_producer, flashHandle);
    nuInt_time = GetFirstFlashTime(*flashHandle);
  }

  //Apply SCE correction to simulated particle's start position
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(nuInt_time) +
   detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
  double xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
  auto sce_offset = SCE->GetPosOffsets(geo::Point_t(nuVtxPos.X(), nuVtxPos.Y(), nuVtxPos.Z()));
  double simVtxX = nuVtxPos.X() + sce_offset.X() - xtimeoffset - 0.6;
  double simVtxY = nuVtxPos.Y() - sce_offset.Y();
  double simVtxZ = nuVtxPos.Z() - sce_offset.Z();

  //Set start position/time and 4-momentum for simulated electron to muon reco values
  double leptonE = electronE;
  if(m_runMuonTest) leptonE = recoMuE;
  //TLorentzVector vtx(nuVtxPos.X(), nuVtxPos.Y(), nuVtxPos.Z(), nuInt_time);
  TLorentzVector vtx(simVtxX, simVtxY, simVtxZ, nuInt_time);
  TLorentzVector mom(recoMuPx, recoMuPy, recoMuPz, leptonE);
  simElec.AddTrajectoryPoint(vtx, mom);

  //Create MCTruth object for simulated electron, write it to art event
  simb::MCTruth truth;
  truth.Add(simElec);
  auto outputTruth = std::make_unique< std::vector<simb::MCTruth> >();
  outputTruth->emplace_back(std::move(truth));
  e.put(std::move(outputTruth));

}


void MCElectronProducer::endJob(){
  std::cout << "No flash found & MC electron start time was set to middle of beam gate in "
   << m_noFlash_counter << " events" << std::endl;
  std::cout << "Flash time found outside of beam gate & MC electron start time was reset "
   << "to gate boundaries in" << m_adjFlash_counter << " events" << std::endl;
}

DEFINE_ART_MODULE(MCElectronProducer)
