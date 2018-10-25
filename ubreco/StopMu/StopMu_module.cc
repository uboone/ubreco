////////////////////////////////////////////////////////////////////////
// Class:       StopMu
// Plugin Type: analyzer (art v2_11_03)
// File:        StopMu_module.cc
//
// Generated at Sun Oct 21 21:38:30 2018 by David Caratelli using cetskelgen
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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class StopMu;


class StopMu : public art::EDAnalyzer {
public:
  explicit StopMu(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StopMu(StopMu const &) = delete;
  StopMu(StopMu &&) = delete;
  StopMu & operator = (StopMu const &) = delete;
  StopMu & operator = (StopMu &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // Additional functions
  float yzDistance(float y1, float z1, float y2, float z2);
  void clear();
  void fillCalorimetry(int pl, std::vector<double> dqdx, std::vector<double> rr);
  bool insideTPCvolume(float x, float y, float z);
  float getPitch(const TVector3 &direction, const int &pl);

private:

  std::string fTrkProducer, fCaloProducer, fMCProducer, fMCParticleLabel, fOpticalFlashFinderLabel;

  TTree* _reco_tree;
  float _trk_len;
  float _trk_start_x;
  float _trk_start_y;
  float _trk_start_z;
  float _trk_end_x;
  float _trk_end_y;
  float _trk_end_z;
  float _yz_true_reco_distance;
  float _pitch_u, _pitch_v, _pitch_y;
  std::vector<float> _dqdx_u;
  std::vector<float> _rr_u;
  std::vector<float> _dqdx_v;
  std::vector<float> _rr_v;
  std::vector<float> _dqdx_y;
  std::vector<float> _rr_y;
  float _delta_t_closest_flash;
};


StopMu::StopMu(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fTrkProducer  = p.get<std::string>("TrkProducer", "pandora");
  fCaloProducer = p.get<std::string>("CaloProducer");
  fMCProducer = p.get<std::string>("MCProducer", "generator");
  fMCParticleLabel = p.get<std::string>("MCParticleLabel", "largeant");
  fOpticalFlashFinderLabel = p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashCosmic");

  art::ServiceHandle<art::TFileService> tfs;

  // _mc_tree = tfs->make<TTree>("mc_tree", "MC Track Tree");
  // _mc_tree->Branch("_trk_len",&_mc_len,"trk_len/F");
  // _mc_tree->Branch("_trk_start_x",&_mc_start_x,"trk_start_x/F");
  // _mc_tree->Branch("_trk_start_y",&_mc_start_y,"trk_start_y/F");
  // _mc_tree->Branch("_trk_start_z",&_mc_start_z,"trk_start_z/F");
  // _mc_tree->Branch("_trk_end_x",  &_mc_end_x,  "trk_end_x/F"  );
  // _mc_tree->Branch("_trk_end_y",  &_mc_end_y,  "trk_end_y/F"  );
  // _mc_tree->Branch("_trk_end_z",  &_mc_end_z,  "trk_end_z/F"  );
  // _mc_tree->Branch("_dqdx_u","std::vector<float>",&_dqdx_u);
  // _mc_tree->Branch("_dqdx_v","std::vector<float>",&_dqdx_v);
  // _mc_tree->Branch("_dqdx_y","std::vector<float>",&_dqdx_y);
  // _mc_tree->Branch("_rr_v",  "std::vector<float>",&_rr_v  );

  _reco_tree = tfs->make<TTree>("reco_tree", "Reco Track Tree");
  _reco_tree->Branch("_trk_len",&_trk_len,"trk_len/F");
  _reco_tree->Branch("_trk_start_x",&_trk_start_x,"trk_start_x/F");
  _reco_tree->Branch("_trk_start_y",&_trk_start_y,"trk_start_y/F");
  _reco_tree->Branch("_trk_start_z",&_trk_start_z,"trk_start_z/F");
  _reco_tree->Branch("_trk_end_x",  &_trk_end_x,  "trk_end_x/F"  );
  _reco_tree->Branch("_trk_end_y",  &_trk_end_y,  "trk_end_y/F"  );
  _reco_tree->Branch("_trk_end_z",  &_trk_end_z,  "trk_end_z/F"  );
  _reco_tree->Branch("_yz_true_reco_distance",  &_yz_true_reco_distance,  "yz_true_reco_distance/F"  );
  _reco_tree->Branch("_pitch_u",  &_pitch_u,  "pitch_u/F"  );
  _reco_tree->Branch("_pitch_v",  &_pitch_v,  "pitch_v/F"  );
  _reco_tree->Branch("_pitch_y",  &_pitch_y,  "pitch_y/F"  );
  _reco_tree->Branch("_dqdx_u","std::vector<float>",&_dqdx_u);
  _reco_tree->Branch("_dqdx_v","std::vector<float>",&_dqdx_v);
  _reco_tree->Branch("_dqdx_y","std::vector<float>",&_dqdx_y);
  _reco_tree->Branch("_rr_u",  "std::vector<float>",&_rr_u  );
  _reco_tree->Branch("_rr_v",  "std::vector<float>",&_rr_v  );
  _reco_tree->Branch("_rr_y",  "std::vector<float>",&_rr_y  );
  _reco_tree->Branch("_delta_t_closest_flash",  &_delta_t_closest_flash,  "delta_t_closest_flash/F"  );

}

void StopMu::analyze(art::Event const & e)
{
  // consider only events with an interaction outside of the TPC
  auto const &generator_handle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCProducer);
  auto const &generator(*generator_handle);
  float _true_vx=-1000, _true_vy=-1000, _true_vz=-1000;
  size_t n_nu = 0;
  for (auto &gen : generator)
  {
    if (gen.Origin() == simb::kBeamNeutrino)
    {
      n_nu = 1;
      _true_vx = gen.GetNeutrino().Nu().Vx();
      _true_vy = gen.GetNeutrino().Nu().Vy();
      _true_vz = gen.GetNeutrino().Nu().Vz();
      break; // In case of events with more than one neutrino (2% of the total) we take for the moment only the first one
    }
  }

  if (n_nu==0) return;

  bool inTPCvolume = insideTPCvolume(_true_vx, _true_vy, _true_vz);
  if (inTPCvolume == true)
  {
    return;
  }

  auto const &mcparticles_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);
  auto const &mcparticles(*mcparticles_handle);

  std::vector<float> mc_muon_end_y, mc_muon_end_z;

  for (auto &mcparticle : mcparticles)
  {
    if (mcparticle.Process() == "primary" &&
        mcparticle.StatusCode() == 1 &&
        abs(mcparticle.PdgCode()) == 13)
    {
      if (mcparticle.EndE() == mcparticle.Mass())
      {
        // _mc_start_x = mcparticle.StartX();
        // _mc_start_y = mcparticle.StartY();
        // _mc_start_z = mcparticle.StartZ();
        // _mc_end_x = mcparticle.EndX();
        float _mc_end_y = mcparticle.EndY();
        float _mc_end_z = mcparticle.EndZ();

        mc_muon_end_y.push_back(_mc_end_y);
        mc_muon_end_z.push_back(_mc_end_z);
      }
    }
  }

  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // grab calorimetry objects associated to tracks
  art::FindMany<anab::Calorimetry> trk_calo_assn_v(trk_h, e, fCaloProducer);

  for (size_t t=0; t < trk_h->size(); t++)
  {
    clear();
    auto const& trk = trk_h->at(t);
    auto const& beg = trk.Vertex();
    auto const& end = trk.End();

    _trk_len = trk.Length();
    _trk_start_x = beg.X();
    _trk_start_y = beg.Y();
    _trk_start_z = beg.Z();
    _trk_end_x   = end.X();
    _trk_end_y   = end.Y();
    _trk_end_z   = end.Z();

    // look for the closest mc muon track
    std::vector<float> yz_distances;
    for (size_t k=0; k<mc_muon_end_y.size(); k++)
    {
      float aux_distance = yzDistance(_trk_end_y, _trk_end_z, mc_muon_end_y[k], mc_muon_end_z[k]);
      yz_distances.push_back(aux_distance);
    }
    _yz_true_reco_distance = *min_element(yz_distances.begin(), yz_distances.end());

    TVector3 track_direction(_trk_end_x - _trk_start_x,
                             _trk_end_y - _trk_start_y,
                             _trk_end_z - _trk_start_z);
    _pitch_u = getPitch(track_direction, 0);
    _pitch_v = getPitch(track_direction, 1);
    _pitch_y = getPitch(track_direction, 2);
    // fill calorimetry info for this track
    // grab the associated calorimetry object
    const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(t);
    for (size_t pl=0; pl < Calo_v.size(); pl++)
    {
      auto const& calo = Calo_v.at(pl);
      auto const& plane = calo->PlaneID().Plane;
      auto const& dqdx = calo->dQdx();
      auto const& rr   = calo->ResidualRange();
      fillCalorimetry(plane, dqdx, rr);
    }


    // closest flash time
    art::InputTag optical_tag_simple{fOpticalFlashFinderLabel};
    auto const &optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>(optical_tag_simple);

    for (size_t f=0; f < optical_handle->size(); f++)
    {
      std::cout << "iteration beg" << f << " delta t " << _delta_t_closest_flash << std::endl;
      auto const& flash = optical_handle->at(f);
      float ttrk = fabs(_trk_end_x / 0.1114 - flash.Time());
      if (ttrk < _delta_t_closest_flash) _delta_t_closest_flash = ttrk;
      ttrk = fabs((_trk_end_x-256.) / 0.1114 - flash.Time());
      if (ttrk < _delta_t_closest_flash) _delta_t_closest_flash = ttrk;

      std::cout << "iteration " << f << " delta t " << _delta_t_closest_flash << " ttrk " << ttrk << std::endl;
    }
    _reco_tree->Fill();
  }
  return;
}

float StopMu::yzDistance(float y1, float z1, float y2, float z2)
{
  return sqrt((y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

void StopMu::clear()
{
  _dqdx_u.clear();
  _rr_u.clear();
  _dqdx_v.clear();
  _rr_v.clear();
  _dqdx_y.clear();
  _rr_y.clear();

  _delta_t_closest_flash = 10000.;
}

void StopMu::fillCalorimetry(int pl, std::vector<double> dqdx, std::vector<double> rr)
{
  if (pl==0)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_u.push_back((float)dqdx[n]);
    _rr_u.push_back(  (float)rr[n]  );
    }
  }
  else if (pl==1)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_v.push_back((float)dqdx[n]);
    _rr_v.push_back(  (float)rr[n]  );
    }
  }
  else if (pl==2)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_y.push_back((float)dqdx[n]);
    _rr_y.push_back(  (float)rr[n]  );
    }
  }
}

bool StopMu::insideTPCvolume(float x, float y, float z)
{
  // hardcoded
  if (x>0 && x<256.35 && y>-116.5 && y<116.5 && z>0 && z<1036.8)
  {
    return true;
  }
  else
  {
    return false;
  }
}

float StopMu::getPitch(const TVector3 &direction, const int &pl)
{
  // prepare a direction vector for the plane
  TVector3 wireDir = {0., 0., 0.};
  // the direction of the plane is the vector uniting two consecutive wires
  // such that this vector is perpendicular to both wires
  // basically this is the vector perpendicular to the wire length direction,
  // and still in the wire-plane direction
  if (pl == 0)
    wireDir = {0., -sqrt(3) / 2., 1 / 2.};
  else if (pl == 1)
    wireDir = {0., sqrt(3) / 2., 1 / 2.};
  else if (pl == 2)
    wireDir = {0., 0., 1.};

  // cosine between shower direction and plane direction gives the factor
  // by which to divide 0.3, the minimum wire-spacing
  double minWireSpacing = 0.3;
  double cos = wireDir.Dot(direction);

  cos /= (wireDir.Mag() * direction.Mag());
  // if cosine is 0 the direction is perpendicular and the wire-spacing is
  // infinite
  if (cos == 0)
    return std::numeric_limits<double>::max();
  float pitch = minWireSpacing / cos;
  return pitch;
}

void StopMu::beginJob()
{
  // Implementation of optional member function here.
}

void StopMu::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(StopMu)
