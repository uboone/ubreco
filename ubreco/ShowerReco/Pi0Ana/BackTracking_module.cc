////////////////////////////////////////////////////////////////////////
// Class:       BackTracking
// Plugin Type: analyzer (art v2_11_03)
// File:        BackTracking_module.cc
//
// Generated at Wed Oct 17 14:05:32 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

//"larsoft" object includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

class BackTracking;


class BackTracking : public art::EDAnalyzer {
public:
  explicit BackTracking(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BackTracking(BackTracking const &) = delete;
  BackTracking(BackTracking &&) = delete;
  BackTracking & operator = (BackTracking const &) = delete;
  BackTracking & operator = (BackTracking &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  std::string fHitProducer, fBacktrackTag;

  std::map<size_t,size_t> HitMapping(const art::ValidHandle<std::vector<recob::Hit>> hitv1, const art::ValidHandle<std::vector<recob::Hit>> hitv2);

  std::map<size_t, std::vector<unsigned int> > GetMCShowerInfo(const art::ValidHandle<std::vector<simb::MCTruth> > mct_h, const art::ValidHandle<std::vector<sim::MCShower> > mcs_h);

  TTree* _tree;
  float _shr_etru;
  float _shr_edep;
  float _shr_Ehit;
  float _shr_ehit;
  float _shr_qhit;

  TTree* _hit_tree;
  float _ihit;
  int   _wire;
  float _ehit;
  float _Ehit;
  float _ehitr;
  float _Ehitr;

};


BackTracking::BackTracking(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  
  fHitProducer  = p.get<std::string>("HitProducer");
  fBacktrackTag = p.get<std::string>("BacktrackTag");

  art::ServiceHandle<art::TFileService> tfs;
  
  _tree = tfs->make<TTree>("_tree","shower completeness ttree");
  _tree->Branch("_shr_etru",&_shr_etru,"shr_etru/F"); // true energy of the electron
  _tree->Branch("_shr_edep",&_shr_edep,"shr_edep/F"); // true deposited energy by the shower
  _tree->Branch("_shr_Ehit",&_shr_Ehit,"shr_Ehit/F"); // truth energy measured from backtracker
  _tree->Branch("_shr_ehit",&_shr_ehit,"shr_ehit/F"); // truth charge measured from backtracker
  _tree->Branch("_shr_qhit",&_shr_qhit,"shr_qhit/F"); // reconstructed charge from hits -> MeV via 198. gain and 0.577 recomb.

  _hit_tree = tfs->make<TTree>("_hit_tree","hit ttree");
  _hit_tree->Branch("_ihit",&_ihit,"ihit/F"); // hit integral [ADCs]
  _hit_tree->Branch("_ehit",&_ehit,"ehit/F"); // hit electrons, from backtracker
  _hit_tree->Branch("_Ehit",&_Ehit,"Ehit/F"); // hit energy, from backtracker
  _hit_tree->Branch("_ehitr",&_ehitr,"ehitr/F"); // hit electrons, from reconstruction
  _hit_tree->Branch("_Ehitr",&_Ehitr,"Ehitr/F"); // hit energy, from reconstruction
  _hit_tree->Branch("_wire",&_wire,"wire/I");

}

void BackTracking::analyze(art::Event const & e)
{
  // Implementation of required member function here.


  // load mcshowers & mctruth
  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
  // load gaushits that have backtracking info
  auto const& gaushit_h = e.getValidHandle<std::vector<recob::Hit> > ("gaushit");
  // load hits we want to use
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitProducer);

  std::cout << "there are " << hit_h->size() << " hits" << std::endl;

  art::InputTag BacktrackTag { fBacktrackTag };

  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_h(gaushit_h,e,BacktrackTag);
  //std::cout << "ass : " << backtrack_h.size() << std::endl;
  //art::Handle<art::Assns<recob::Hit,simb::MCParticle,anab::BackTrackerHitMatchingData> > backtrack_h;
  //e.getByLabel(fBacktrackTag,backtrack_h);
  
  // create map:
  auto hitMap = HitMapping(hit_h,gaushit_h);


  auto MCShowerInfo = GetMCShowerInfo(mct_h,mcs_h);

  // for each MCshower, study charge-matching to hits of interest
  for (auto const& shrinfo : MCShowerInfo){

    auto shrIDX = shrinfo.first;
    auto trackid_v = shrinfo.second;
    
    // deposited energy
    auto edep = mcs_h->at(shrIDX).DetProfile().E();
    std::cout  << "shower has " << edep << " energy" << std::endl;

    _shr_edep = edep;
    _shr_etru = mcs_h->at(shrIDX).Start().E();

    double collectedEnergy   = 0.;
    double collectedCharge   = 0.;
    double collectedIntegral = 0.;

    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;  

    for (size_t i=0; i < hit_h->size(); i++) {
      
      // get associations for this hit
      auto const hit = hit_h->at(i);
      if (hit.WireID().Plane != 2) continue;

      _ihit = hit.Integral();
      _wire    = hit.WireID().Wire;
      _ehit = 0.;
      _Ehit = 0.;


      
      // now get the gaushit
      auto const gaushit_idx = hitMap[i];
      //std::cout << "idx gaushit : " << gaushit_idx << "\t at pt " << i << std::endl;
      //auto const gaushit = gaushit_h->at(gaushit_idx);
      
      particle_vec.clear(); match_vec.clear();
      
      backtrack_h.get(gaushit_idx, particle_vec, match_vec);
      
      // does this trackID match that of the MCShower?
      bool matchedID = false;

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
	
	auto mctrkid = particle_vec.at(i_p)->TrackId();

	_ehit += match_vec[i_p]->numElectrons;
	_Ehit += match_vec[i_p]->energy;

	// does this trackID match that of the MCShower?
	for (auto const& shrtrkid : trackid_v)  
	  if ( shrtrkid == (unsigned int)mctrkid ) { matchedID = true; break; }
	
	if (matchedID) {

	  collectedEnergy += match_vec[i_p]->energy;
	  collectedCharge += match_vec[i_p]->numElectrons;

	}// if the mcparticle belongs to the shower

      }// for all particles associated to this hit

      if (matchedID) {
	_ehitr = _ihit * 198;
	_Ehitr = _ehitr * 23.6 * (1e-6) / 0.577;
	collectedIntegral += _Ehitr;
	_hit_tree->Fill();
      }

    }// for all hits

    _shr_Ehit = collectedEnergy;
    _shr_ehit = collectedCharge;
    _shr_qhit = collectedIntegral;

    _tree->Fill();

    std::cout << " collected energy from hits : " << collectedEnergy << std::endl;
      
  }// for all MCShowers
  
  return;
}

void BackTracking::beginJob()
{
  // Implementation of optional member function here.
}

std::map<size_t,size_t> BackTracking::HitMapping(const art::ValidHandle<std::vector<recob::Hit>> hitv1, const art::ValidHandle<std::vector<recob::Hit>> hitv2) {

  //if the hits are not gaushit, create a map connecting the hit index to the gaushit hit index.
  std::map<size_t,size_t> HitHitMap;
  for (size_t h0=0; h0 < hitv1->size(); h0++){
      auto const& hitcp = hitv1->at(h0);
      bool foundmatch = false;
      for (size_t h1=0; h1 < hitv2->size(); h1++){
	auto const& hitoriginal = hitv2->at(h1);
	if ( (hitcp.PeakTime() == hitoriginal.PeakTime()) && (hitcp.WireID().Wire == hitoriginal.WireID().Wire) ){
	  HitHitMap[h0] = h1;
	  foundmatch = true;
	  break;
	}
      }
      if (foundmatch == false)
	std::cout << "\t ERROR no match found!" << std::endl;
  }

  std::cout << "there are " << hitv1->size() << " hits " << std::endl;
  std::cout << "there are " << HitHitMap.size() << " map elements" << std::endl;

  return HitHitMap;
  
}// end of hit-mapping function

std::map<size_t, std::vector<unsigned int> > BackTracking::GetMCShowerInfo(const art::ValidHandle<std::vector<simb::MCTruth> > mct_h, const art::ValidHandle<std::vector<sim::MCShower> > mcs_h) {

  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  Double_t xyz[3] = {};

  // save the trackID of the pi0
  unsigned int pi0trkId = 0;
  int npi0 = 0;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      npi0 += 1;
      pi0trkId = part.TrackId();
      xyz[0] = part.Trajectory().X(0);
      xyz[1] = part.Trajectory().Y(0);
      xyz[2] = part.Trajectory().Z(0);
      break;
    }
  }

  std::cout << "pi0 track id : " << pi0trkId << std::endl;

  // loop through MCShowers and identify those originating from the pi0
  // map connecting mcshower index to track ID vector for all e+/e- in MCShower
  std::map<size_t, std::vector<unsigned int> > event_shower_map;
  // map connecting e+/e- trackID in mcshower to mcshower index
  //std::map<unsigned int, size_t> event_mcpart_map;

  for (size_t i=0; i < mcs_h->size(); i++) {
    auto const& mcs = mcs_h->at(i);

    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (xyz[0] - x) * (xyz[0] - x) ) +
		     ( (xyz[1] - y) * (xyz[1] - y) ) +
		     ( (xyz[2] - z) * (xyz[2] - z) ) );

    if (d < 0.01) {
      std::vector<unsigned int> shrtrackIDs = mcs.DaughterTrackID();
      shrtrackIDs.push_back( mcs.TrackID() );
      // get daughter track IDs:
      auto daughterIDs = mcs.DaughterTrackID();
      for (auto const& id : daughterIDs)
	if (id != mcs.TrackID()) { shrtrackIDs.push_back(id); }
      
      event_shower_map[ i ] = shrtrackIDs;
      std::cout << "shower " << i << " has "<< shrtrackIDs.size() << " associated MCParticles" << std::endl;
    }// if mcshower matched to pi0
  }// for all mcshowers
  std::cout << "found " << event_shower_map.size() << " mcshowers associated to the pi0" << std::endl;

  return event_shower_map;
}

void BackTracking::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(BackTracking)
