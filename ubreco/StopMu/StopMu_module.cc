////////////////////////////////////////////////////////////////////////
// Class:       StopMu
// Plugin Type: analyzer (art v2_11_03)
// File:        StopMu_module.cc
//
// Generated at Sun Oct 21 21:38:30 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
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

#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"


#include "lardataobj/RecoBase/OpFlash.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


/**
   \class TruncMean
   The truncated mean class allows to compute the following quantities
   1) the truncated mean profile of an ordered vector of values, such as
   the charge profile along a particle's track.
   To create such a profile use the function CalcTruncMeanProfile()
   2) Get the truncated mean value of a distribution. This function
   iteratively hones in on the truncated mean of a distribution by
   updating the mean and cutting the tails after each iteration.
   For this functionality use CalcIterativeTruncMean()
   doxygen documentation!
*/

static const double kINVALID_DOUBLE = std::numeric_limits<double>::max();

class TruncMean;

class TruncMean{

 public:

  /// Default constructor
  TruncMean(){}

  /// Default destructor
  ~TruncMean(){}

  /**
     @brief Given residual range and dq vectors return truncated local dq.
     Input vectors are assumed to be match pair-wise (nth entry in rr_v
     corresponds to nth entry in dq_v vector).
     Input rr_v values are also assumed to be ordered: monotonically increasing
     or decreasing.
     For every dq value a truncated linear dq value is calculated as follows:
     0) all dq values within a rr range set by the class variable _rad are selected.
     1) the median and rms of these values is calculated.
     2) the subset of local dq values within the range [median-rms, median+rms] is selected.
     3) the resulting local truncated dq is the average of this truncated subset.
     @input std::vector<double> rr_v -> vector of x-axis coordinates (i.e. position for track profile)
     @input std::vector<double> dq_v -> vector of measured values for which truncated profile is requested
     (i.e. charge profile of a track)
     @input std::vector<double> dq_trunc_v -> passed by reference -> output stored here
     @input double nsigma -> optional parameter, number of sigma to keep around RMS for TM calculation
  */
  void CalcTruncMeanProfile(const std::vector<double>& rr_v, const std::vector<double>& dq_v,
			    std::vector<double>& dq_trunc_v, const double& nsigma = 1);

  /**
     @brief Iteratively calculate the truncated mean of a distribution
     @input std::vector<double> v -> vector of values for which truncated mean is asked
     @input size_t nmin -> minimum number of iterations to converge on truncated mean
     @input size_t nmax -> maximum number of iterations to converge on truncated mean
     @input size_t lmin -> minimum number of entries in vector before exiting and returning current value
     @input size_t currentiteration -> current iteration
     @input double convergencelimit -> fractional difference between successive iterations
     under which the iteration is completed, provided nmin iterations have occurred.
     @input nsigma -> number of sigma around the median value to keep when the distribution is trimmed.
   */
  double CalcIterativeTruncMean(std::vector<double> v, const size_t& nmin,
			       const size_t& nmax, const size_t& currentiteration,
			       const size_t& lmin,
			       const double& convergencelimit,
			       const double& nsigma, const double& oldmed = kINVALID_DOUBLE);

  /**
     @brief Set the smearing radius over which to take hits for truncated mean computaton.
   */
  void setRadius(const double& rad) { _rad = rad; }

 private:

  double Mean  (const std::vector<double>& v);
  double Median(const std::vector<double>& v);
  double RMS   (const std::vector<double>& v);

  /**
     Smearing radius over which charge from neighboring hits is scanned to calculate local
     truncated mean
   */
  double _rad;

};

double TruncMean::CalcIterativeTruncMean(std::vector<double> v, const size_t& nmin,
					const size_t& nmax, const size_t& currentiteration,
					const size_t& lmin,
					const double& convergencelimit,
					const double& nsigma, const double& oldmed)
{

  auto const& mean = Mean(v);
  auto const& med  = Median(v);
  auto const& rms  = RMS(v);

  // if the vector length is below the lower limit -> return
  if (v.size() < lmin)
    return mean;

  // if we have passed the maximum number of iterations -> return
  if (currentiteration >= nmax)
    return mean;

  // if we passed the minimum number of iterations and the mean is close enough to the old value
  double fracdiff = fabs(med-oldmed) / oldmed;
  if ( (currentiteration >= nmin) && (fracdiff < convergencelimit) )
    return mean;

  // if reached here it means we have to go on for another iteration

  // cutoff tails of distribution surrounding the mean
  // use erase-remove : https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
  // https://stackoverflow.com/questions/17270837/stdvector-removing-elements-which-fulfill-some-conditions
  v.erase( std::remove_if( v.begin(), v.end(),
			   [med,nsigma,rms](const double& x) { return ( (x < (med-nsigma*rms)) || (x > (med+nsigma*rms)) ); }), // lamdda condition for events to be removed
	   v.end());

  return CalcIterativeTruncMean(v, nmin, nmax, lmin, currentiteration+1, convergencelimit, nsigma, med);
}

void TruncMean::CalcTruncMeanProfile(const std::vector<double>& rr_v, const std::vector<double>& dq_v,
				     std::vector<double>& dq_trunc_v, const double& nsigma)
{

  // how many points to sample
  int Nneighbor = (int)(_rad * 3 * 2);

  dq_trunc_v.clear();
  dq_trunc_v.reserve( rr_v.size() );

  int Nmax = dq_v.size()-1;

  for (size_t n=0; n < dq_v.size(); n++) {

    // current residual range
    double rr = rr_v.at(n);

    int nmin = n - Nneighbor;
    int nmax = n + Nneighbor;

    if (nmin < 0) nmin = 0;
    if (nmax > Nmax) nmax = Nmax;

    // vector for local dq values
    std::vector<double> dq_local_v;

    for (int i=nmin; i < nmax; i++) {

      double dr = rr - rr_v[i];
      if (dr < 0) dr *= -1;

      if (dr > _rad) continue;

      dq_local_v.push_back( dq_v[i] );

    }// for all ticks we want to scan

    if (dq_local_v.size() == 0) {
      dq_trunc_v.push_back( dq_v.at(n) );
      continue;
    }

    // calculate median and rms
    double median = Median(dq_local_v);
    double rms    = RMS(dq_local_v);

    double truncated_dq = 0.;
    int npts = 0;
    for (auto const& dq : dq_local_v) {
      if ( ( dq < (median+rms * nsigma) ) && ( dq > (median-rms * nsigma) ) ){
	truncated_dq += dq;
	npts += 1;
      }
    }

    dq_trunc_v.push_back( truncated_dq / npts );
  }// for all values

  return;
}

double TruncMean::Mean(const std::vector<double>& v)
{

  double mean = 0.;
  for (auto const& n : v) mean += n;
  mean /= v.size();

  return mean;
}

double TruncMean::Median(const std::vector<double>& v)
{

  if (v.size() == 1) return v[0];

  std::vector<double> vcpy = v;

  std::sort(vcpy.begin(), vcpy.end());

  double median = vcpy[ vcpy.size() / 2 ];

  return median;
}

double TruncMean::RMS(const std::vector<double>& v)
{

  double avg = 0.;
  for (auto const& val : v) avg += val;
  avg /= v.size();
  double rms = 0.;
  for (auto const& val : v) rms += (val-avg)*(val-avg);
  rms = sqrt( rms / ( v.size() -  1 ) );

  return rms;
}


class StopMu;


class StopMu : public art::EDProducer {
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
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // Additional functions
  double yzDistance(double y1, double z1, double y2, double z2);
  double xyzDistance(double x1, double y1, double z1, double x2, double y2, double z2);
  int findEndPointHits(const double& endx, const double& endy, const double& endz,
		       const std::vector<art::Ptr<recob::Hit> >& ass_hit_v,
		       const art::ValidHandle<std::vector<recob::Hit> >& gaushit_h,
		       const unsigned int& pl);
    void clear();
  void fillCalorimetry(int pl, std::vector<float> dqdx, std::vector<float> rr);//, std::vector<TVector3> xyz);
  bool insideTPCvolume(double x, double y, double z);
  double getPitch(const TVector3 &direction, const int &pl);
  void shiftTruePosition(double true_point[3], double true_time, double true_point_shifted[3]);

private:

  std::string fTrkProducer, fCaloProducer, fMCProducer, fMCParticleLabel, fOpticalFlashFinderLabelB, fOpticalFlashFinderLabelC;
  bool fGeoCuts, fUseTruth;
  double fdT, fMinX, fMaxX, fMinZ, fMaxZ, fMinY, fMinLen;
  double fZdeadStart, fZdeadEnd;
  double fAbsXspan; // to remove anode-to-cathode crossing tracks

  double _wire2cm, _time2cm;

  TruncMean _tmean;

  TTree* _reco_tree;
  int _run, _sub, _evt;

  double _trk_len;

  double _trk_start_x;
  double _trk_start_y;
  double _trk_start_z;

  double _px, _py, _pz; // direction vector of reconstructed track measured from start to end
  
  double _trk_end_x;
  double _trk_end_y;
  double _trk_end_z;

  double _mc_trk_end_x;
  double _mc_trk_end_y;
  double _mc_trk_end_z;

  double _offsetx, _offsety, _offsetz;

  int _nhit_endpoint;

  double _xyz_true_reco_distance, _yz_true_reco_distance;
  double _true_energy;
  int _yz_trackid; // track id of closest stopping muon by YZ distance
  double _matchscore; // fraction of track hits associated to a true stopping muon MCParticle
  int   _matchtrackid; // track id of best match stopping muon from backtracker
  int   _pdg; // PDG score of true track
  double _pitch_u, _pitch_v, _pitch_y;
  std::vector<double> _dqdx_u;
  std::vector<double> _dqdx_tm_u;
  std::vector<double> _rr_u;
  std::vector<double> _x_position_u;
  std::vector<double> _y_position_u;
  std::vector<double> _z_position_u;
  std::vector<double> _dqdx_v;
  std::vector<double> _dqdx_tm_v;
  std::vector<double> _rr_v;
  std::vector<double> _x_position_v;
  std::vector<double> _y_position_v;
  std::vector<double> _z_position_v;
  std::vector<double> _dqdx_y;
  std::vector<double> _dqdx_tm_y;
  std::vector<double> _rr_y;
  std::vector<double> _x_position_y;
  std::vector<double> _y_position_y;
  std::vector<double> _z_position_y;
  double _delta_t_closest_flash, _flash_time;
};


StopMu::StopMu(fhicl::ParameterSet const & p)
  :
  EDProducer(p)  // ,
 // More initializers here.
{

  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();

  fTrkProducer  = p.get<std::string>("TrkProducer");
  fCaloProducer = p.get<std::string>("CaloProducer");
  fMCProducer = p.get<std::string>("MCProducer", "generator");
  fMCParticleLabel = p.get<std::string>("MCParticleLabel", "largeant");
  fOpticalFlashFinderLabelC = p.get<std::string>("OpticalFlashFinderLabelC", "simpleFlashCosmic");
  fOpticalFlashFinderLabelB = p.get<std::string>("OpticalFlashFinderLabelB", "simpleFlashBeam");
  fGeoCuts = p.get<bool>("GeoCuts");
  fdT = p.get<double>("dT");
  fMinX = p.get<double>("MinX");
  fMaxX = p.get<double>("MaxX");
  fMinZ = p.get<double>("MinZ");
  fMaxZ = p.get<double>("MaxZ");
  fMinY = p.get<double>("MinY");
  fZdeadStart = p.get<double>("ZdeadStart");
  fZdeadEnd = p.get<double>("ZdeadEnd");
  fAbsXspan = p.get<double>("AbsXspan");
  fMinLen = p.get<double>("MinLen");
  fUseTruth = p.get<bool>("UseTruth");
  
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,0,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  _tmean.setRadius(10);

  art::ServiceHandle<art::TFileService> tfs;

  _reco_tree = tfs->make<TTree>("reco_tree", "Reco Track Tree");

  _reco_tree->Branch("_run",&_run,"run/I");
  _reco_tree->Branch("_sub",&_sub,"sub/I");
  _reco_tree->Branch("_evt",&_evt,"evt/I");

  _reco_tree->Branch("_trk_len",&_trk_len,"trk_len/D");

  _reco_tree->Branch("_px",&_px,"px/D");
  _reco_tree->Branch("_py",&_py,"py/D");
  _reco_tree->Branch("_pz",&_pz,"pz/D");

  _reco_tree->Branch("_trk_start_x",&_trk_start_x,"trk_start_x/D");
  _reco_tree->Branch("_trk_start_y",&_trk_start_y,"trk_start_y/D");
  _reco_tree->Branch("_trk_start_z",&_trk_start_z,"trk_start_z/D");

  _reco_tree->Branch("_trk_end_x",  &_trk_end_x,  "trk_end_x/D"  );
  _reco_tree->Branch("_trk_end_y",  &_trk_end_y,  "trk_end_y/D"  );
  _reco_tree->Branch("_trk_end_z",  &_trk_end_z,  "trk_end_z/D"  );

  // truth variables
  _reco_tree->Branch("_xyz_true_reco_distance",  &_xyz_true_reco_distance,  "xyz_true_reco_distance/D"  );
  _reco_tree->Branch("_yz_true_reco_distance",  &_yz_true_reco_distance,  "yz_true_reco_distance/D"  );
  _reco_tree->Branch("_true_energy",&_true_energy,"true_energy/D");
  _reco_tree->Branch("_yz_trackid",  &_yz_trackid,  "yz_trackid/I"  );
  _reco_tree->Branch("_matchtrackid",  &_matchtrackid,  "matchtrackid/I"  );
  _reco_tree->Branch("_matchscore",  &_matchscore,  "matchscore/D"  );
  _reco_tree->Branch("_pdg",  &_pdg,  "pdg/I"  );
  _reco_tree->Branch("_mc_trk_end_x",  &_mc_trk_end_x,  "mc_trk_end_x/D"  );
  _reco_tree->Branch("_mc_trk_end_y",  &_mc_trk_end_y,  "mc_trk_end_y/D"  );
  _reco_tree->Branch("_mc_trk_end_z",  &_mc_trk_end_z,  "mc_trk_end_z/D"  );

  _reco_tree->Branch("_offsetx",  &_offsetx,  "offsetx/D"  );
  _reco_tree->Branch("_offsety",  &_offsety,  "offsety/D"  );
  _reco_tree->Branch("_offsetz",  &_offsetz,  "offsetz/D"  );

  //_reco_tree->Branch("_pitch_u",  &_pitch_u,  "pitch_u/D"  );
  //_reco_tree->Branch("_pitch_v",  &_pitch_v,  "pitch_v/D"  );
  _reco_tree->Branch("_pitch_y",  &_pitch_y,  "pitch_y/D"  );

  //_reco_tree->Branch("_dqdx_u","std::vector<double>",&_dqdx_u);
  //_reco_tree->Branch("_dqdx_v","std::vector<double>",&_dqdx_v);
  _reco_tree->Branch("_dqdx_y","std::vector<double>",&_dqdx_y);

  //_reco_tree->Branch("_dqdx_tm_u","std::vector<double>",&_dqdx_tm_u);
  //_reco_tree->Branch("_dqdx_tm_v","std::vector<double>",&_dqdx_tm_v);
  _reco_tree->Branch("_dqdx_tm_y","std::vector<double>",&_dqdx_tm_y);

  //_reco_tree->Branch("_rr_u",  "std::vector<double>",&_rr_u  );
  //_reco_tree->Branch("_rr_v",  "std::vector<double>",&_rr_v  );
  _reco_tree->Branch("_rr_y",  "std::vector<double>",&_rr_y  );

  //_reco_tree->Branch("_x_position_u",  "std::vector<double>", &_x_position_u  );
  //_reco_tree->Branch("_x_position_v",  "std::vector<double>", &_x_position_v  );
  _reco_tree->Branch("_x_position_y",  "std::vector<double>", &_x_position_y  );

  //_reco_tree->Branch("_y_position_u",  "std::vector<double>", &_y_position_u  );
  //_reco_tree->Branch("_y_position_v",  "std::vector<double>", &_y_position_v  );
  _reco_tree->Branch("_y_position_y",  "std::vector<double>", &_y_position_y  );

  //_reco_tree->Branch("_z_position_u",  "std::vector<double>", &_z_position_u  );
  //_reco_tree->Branch("_z_position_v",  "std::vector<double>", &_z_position_v  );
  _reco_tree->Branch("_z_position_y",  "std::vector<double>", &_z_position_y  );

  _reco_tree->Branch("_delta_t_closest_flash",  &_delta_t_closest_flash,  "delta_t_closest_flash/D"  );
  _reco_tree->Branch("_flash_time",  &_flash_time,  "flash_time/D"  );

  _reco_tree->Branch("_nhit_endpoint",&_nhit_endpoint,"nhit_endpoint/I");

}

void StopMu::produce(art::Event& e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // produce CosmicTag and association to tracks
  std::unique_ptr< std::vector< anab::CosmicTag > >               CT_v( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >  trk_CT_assn_v( new art::Assns<recob::Track, anab::CosmicTag> );

  // truth muon info for stopping muons
  std::vector<double> mc_muon_time;
  std::vector<double> mc_muon_end_y, mc_muon_end_z, mc_muon_end_x;
  std::vector<int>    mc_muon_pdg;
  std::vector<size_t> stop_mu_trackid_v;
  std::vector<double> mc_muon_energy_v;


  if (fUseTruth) {


    // consider only events with an interaction outside of the TPC
    auto const &generator_handle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCProducer);
    auto const &generator(*generator_handle);
    double _true_vx=-1000, _true_vy=-1000, _true_vz=-1000;
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
    
    if (n_nu==0) {
      e.put( std::move(CT_v) );
      e.put( std::move(trk_CT_assn_v) );
      return;
    }
    
    bool inTPCvolume = insideTPCvolume(_true_vx, _true_vy, _true_vz);
    if (inTPCvolume == true)
      {
	e.put( std::move(CT_v) );
	e.put( std::move(trk_CT_assn_v) );
	return;
      }
    
    auto const &mcparticles_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);
    auto const &mcparticles(*mcparticles_handle);
    
    for (auto &mcparticle : mcparticles)
      {
	if (mcparticle.Process() == "primary" &&
	    mcparticle.StatusCode() == 1 &&
	    abs(mcparticle.PdgCode()) == 13 &&
	    (mcparticle.EndE() == mcparticle.Mass()) &&
	    insideTPCvolume(mcparticle.EndX(), mcparticle.EndY(), mcparticle.EndZ()) == true)
	  {
	    // find first step in the TPC
	    unsigned int npoints = mcparticle.NumberTrajectoryPoints();
	    _true_energy = 0.;
	    for (unsigned int s=0; s < npoints; s++){
	      auto pt = mcparticle.Position(s);
	      if (insideTPCvolume( pt.X(), pt.Y(), pt.Z() ) == true) { 
		_true_energy = mcparticle.E(s);
		break;
	      }// if step in TPC
	    }// for all steps
	    double true_end[3], true_end_shifted[3];
	    true_end[0] = mcparticle.EndX();
	    true_end[1] = mcparticle.EndY();
	    true_end[2] = mcparticle.EndZ();
	    shiftTruePosition(true_end, mcparticle.EndT(), true_end_shifted);
	    mc_muon_end_x.push_back(true_end_shifted[0]);
	    mc_muon_end_y.push_back(true_end_shifted[1]);
	    mc_muon_end_z.push_back(true_end_shifted[2]);
	    stop_mu_trackid_v.push_back( mcparticle.TrackId() );
	    mc_muon_pdg.push_back( mcparticle.PdgCode() );
	  }
      }

  }// if we are using truth

  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // load hits
  auto const& gaushit_h = e.getValidHandle<std::vector<recob::Hit> > ("gaushit");
  // grab calorimetry objects associated to tracks
  art::FindMany<anab::Calorimetry> trk_calo_assn_v(trk_h, e, fCaloProducer);
  // grab hits associated to tracks
  art::FindManyP<recob::Hit> trk_hit_assn_v(trk_h, e, fTrkProducer);

  ///* backtrack
  //if (fUsTruth)
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_h(gaushit_h,e,"gaushitTruthMatch");  
  //*/
  

  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, trk_h);

  size_t t = 0;


  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();


  for (auto trk: TrkVec){

    t += 1;

    clear();

    auto const& beg = trk->Vertex();
    auto const& end = trk->End();

    _trk_len = trk->Length();

    _trk_start_x = beg.X();
    _trk_start_y = beg.Y();
    _trk_start_z = beg.Z();

    _trk_end_x   = end.X();
    _trk_end_y   = end.Y();
    _trk_end_z   = end.Z();

    if (_trk_end_y > _trk_start_y) continue;

    double trk_mag = sqrt( (_trk_end_x - _trk_start_x) * (_trk_end_x - _trk_start_x) + 
			   (_trk_end_y - _trk_start_y) * (_trk_end_y - _trk_start_y) + 
			   (_trk_end_z - _trk_start_z) * (_trk_end_z - _trk_start_z) );
    
    _px = (_trk_end_x - _trk_start_x) / trk_mag;
    _py = (_trk_end_y - _trk_start_y) / trk_mag;
    _pz = (_trk_end_z - _trk_start_z) / trk_mag;

    // look for the closest mc muon track
    _yz_true_reco_distance = 1500.;
    _yz_trackid = 0;
    for (size_t k=0; k<mc_muon_end_y.size(); k++)
    {

      auto offset = SCE->GetPosOffsets(geo::Point_t(mc_muon_end_x[k],mc_muon_end_y[k],mc_muon_end_z[k]));
      
      double this_distance = yzDistance(_trk_end_y, _trk_end_z, mc_muon_end_y[k], mc_muon_end_z[k]);
      if (this_distance < _yz_true_reco_distance)
	{
	  _yz_true_reco_distance = this_distance;
	  _xyz_true_reco_distance = xyzDistance(_trk_end_x, _trk_end_y, _trk_end_z, mc_muon_end_x[k], mc_muon_end_y[k], mc_muon_end_z[k]);
	  _yz_trackid = stop_mu_trackid_v.at(k);
	  _offsetx = offset.X();
	  _offsety = offset.Y();
	  _offsetz = offset.Z();
	  _mc_trk_end_x = mc_muon_end_x[k];
	  _mc_trk_end_y = mc_muon_end_y[k];
	  _mc_trk_end_z = mc_muon_end_z[k];
	  _pdg = mc_muon_pdg.at(k);
	}// if closest
    }// for all distances
    
    TVector3 track_direction(_trk_end_x - _trk_start_x,
                             _trk_end_y - _trk_start_y,
                             _trk_end_z - _trk_start_z);
    _pitch_u = getPitch(track_direction, 0);
    _pitch_v = getPitch(track_direction, 1);
    _pitch_y = getPitch(track_direction, 2);
    // fill calorimetry info for this track
    // grab the associated calorimetry object
    const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(t-1);
    for (size_t pl=0; pl < Calo_v.size(); pl++)
    {
      auto const& calo = Calo_v.at(pl);
      auto const& plane = calo->PlaneID().Plane;
      auto const& dqdx  = calo->dQdx();
      auto const& rr    = calo->ResidualRange();
      //auto const& xyz   = calo->XYZ();
      fillCalorimetry(plane, dqdx, rr);//, xyz);
    }

    if (fUseTruth && (backtrack_h.isValid() == true) ) {
      // get associated hits
      const std::vector<art::Ptr<recob::Hit> > ass_hit_v = trk_hit_assn_v.at(t-1);
      
      _matchscore = 0.;
      _matchtrackid = 0;
      
      for (art::Ptr<recob::Hit> hit : ass_hit_v)
	{
	  
	  auto hitidx = hit.key();
	  
	  std::vector<simb::MCParticle const*> particle_vec;
	  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
	  
	  backtrack_h.get(hitidx, particle_vec, match_vec);
	  
	  // does this trackID match that of the MCShower?
	  bool matchedID = false;
	  
	  for(size_t i_p=0; i_p<particle_vec.size(); ++i_p)
	    {
	      
	      auto mctrkid = particle_vec.at(i_p)->TrackId();
	      
	      // does this trackID match that of the MCShower?
	      for (auto const& stopmu_trkid : stop_mu_trackid_v)
		{
		  if ( stopmu_trkid == (unsigned int)mctrkid )
		    {
		      _matchtrackid = stopmu_trkid;
		      matchedID = true;
		      break;
		    }
		}
	    }
	  // this way of filling the score assumes if we find a match it will be with one of the identified true stopping muons
	  if (matchedID)
	    {
	      _matchscore += 1;
	    }// for all particles associated to this hit
	}// for all hits
      
      _matchscore /= ass_hit_v.size();
    }// if using truth and backtracking
    
    // closest flash time
    art::InputTag CFtag{fOpticalFlashFinderLabelC};
    auto const &CF = e.getValidHandle<std::vector<recob::OpFlash>>(CFtag);

    for (size_t f=0; f < CF->size(); f++) {
      //std::cout << "iteration beg" << f << " delta t " << _delta_t_closest_flash << std::endl;
      auto const& flash = CF->at(f);
      double ttrk = _trk_end_x / 0.1114 - flash.Time() - 6.0;
      if (fabs(ttrk) < fabs(_delta_t_closest_flash)) {
	_delta_t_closest_flash = ttrk;
	_flash_time = flash.Time();
      }
      ttrk = (_trk_end_x-256.35) / 0.1114 - flash.Time() - 19.3;
      if (fabs(ttrk) < fabs(_delta_t_closest_flash)) {
	_delta_t_closest_flash = ttrk;
	_flash_time = flash.Time();
      }
    }// for all flashes
    
    // closest flash time
    art::InputTag BFtag{fOpticalFlashFinderLabelB};
    auto const &BF = e.getValidHandle<std::vector<recob::OpFlash>>(BFtag);

    for (size_t f=0; f < BF->size(); f++) {
      //std::cout << "iteration beg" << f << " delta t " << _delta_t_closest_flash << std::endl;
      auto const& flash = BF->at(f);
      double ttrk = _trk_end_x / 0.1114 - flash.Time() - 6.0;
      if (fabs(ttrk) < fabs(_delta_t_closest_flash)) {
	_delta_t_closest_flash = ttrk;
	_flash_time = flash.Time();
      }
      ttrk = (_trk_end_x-256.35) / 0.1114 - flash.Time() - 19.3;
      if (fabs(ttrk) < fabs(_delta_t_closest_flash)) {
	_delta_t_closest_flash = ttrk;
	_flash_time = flash.Time();
      }
    }// for all flashes

    // apply minimum track length requirement 
    if (_trk_len < fMinLen) continue;

    // apply geometric cuts if desired
    if (fGeoCuts) {
      if (fabs(_delta_t_closest_flash) < fdT) continue;
      if ( (_trk_end_x < fMinX) or (_trk_end_x > fMaxX) ) continue;
      if (_trk_end_y < fMinY) continue;
      if ( (_trk_end_z < fMinZ) or (_trk_end_z > fMaxZ) ) continue;
      if ( (_trk_end_z > fZdeadStart) and (_trk_end_z < fZdeadEnd) ) continue;
      if ( fabs(_trk_start_x - _trk_end_x) > fAbsXspan ) continue;
    } 

    // find hits near the track end point, if any
    _nhit_endpoint = findEndPointHits(_trk_end_x, _trk_end_y, _trk_end_z, trk_hit_assn_v.at(t-1), gaushit_h, 2);

    CT_v->emplace_back(1.0);
    util::CreateAssn(*this, e, *CT_v, trk, *trk_CT_assn_v );
    
    _reco_tree->Fill();
  }// for all tracks

  e.put( std::move(CT_v) );
  e.put( std::move(trk_CT_assn_v) );

  return;
}

int StopMu::findEndPointHits(const double& endx, const double& endy, const double& endz,
			     const std::vector<art::Ptr<recob::Hit> >& ass_hit_v,
			     const art::ValidHandle<std::vector<recob::Hit> >& gaushit_h,
			     const unsigned int& pl) {

  int nearby = 0;

  std::cout << "event " << _evt << " with end-point : " << endx << ", " << endy << ", " << endz << std::endl;

  // get end point coordinates in wire/time
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto stopwirecm = geom->WireCoordinate(endy,endz,geo::PlaneID(0,0,pl)) * _wire2cm;
  auto stoptickcm = endx;

  std::cout << " track end point @ [ " << stopwirecm << ", " << stoptickcm << " ]" << std::endl;

  /*
  for (size_t h=0; h < ass_hit_v.size(); h++) {
    auto const& hitptr = ass_hit_v.at(h);
    if (hitptr->WireID().Plane != pl) continue;
    auto timecm = hitptr->PeakTime()    * _time2cm;
    auto wirecm = hitptr->WireID().Wire * _wire2cm;
    std::cout << " track hit @ [ " << wirecm << ", " << timecm << " ]" << std::endl;
  }
  */
  
  for (size_t h=0; h < gaushit_h->size(); h++) {

    // only count hits on the plane of interest
    if (gaushit_h->at(h).WireID().Plane != pl) continue;

    // grab hit index and coordinates
    auto hit = gaushit_h->at(h);

    auto timecm = hit.PeakTime()    * _time2cm - 44.56;
    // if time not within min cm -> continue
    if (fabs(timecm-stoptickcm) > 5.) continue;

    auto wirecm = hit.WireID().Wire * _wire2cm;

    //std::cout << "\t found nearby hit in time @ [ " << wirecm << ", " << timecm << " ] " << std::endl;


    auto dist2D = yzDistance(wirecm,timecm, stopwirecm, stoptickcm);

    // if within a 5 cm radius
    if (dist2D < 5.) {
      //ask if the hit is already in the track hit vector
      bool foundhit = false;
      for (auto const& ass_hit_ptr : ass_hit_v) {
	if (ass_hit_ptr.key() == h) {
	  foundhit = true;
	  break;
	}// if the hit is found in the track hit vector
      }//for all hits in the track
      if (foundhit == false) { nearby += 1; }
    }// if the hit is within 5 cm
  }// for all reconstructed hits
  
  std::cout << "nearby hits : " << nearby << std::endl;
  return nearby;
}

double StopMu::yzDistance(double y1, double z1, double y2, double z2)
{
  return sqrt((y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

double StopMu::xyzDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

void StopMu::clear()
{
  _dqdx_u.clear();
  _dqdx_tm_u.clear();
  _rr_u.clear();
  _x_position_u.clear();
  _y_position_u.clear();
  _z_position_u.clear();

  _dqdx_v.clear();
  _dqdx_tm_v.clear();
  _rr_v.clear();
  _x_position_v.clear();
  _y_position_v.clear();
  _z_position_v.clear();

  _dqdx_y.clear();
  _dqdx_tm_y.clear();
  _rr_y.clear();
  _x_position_y.clear();
  _y_position_y.clear();
  _z_position_y.clear();

  _delta_t_closest_flash = 10000.;
}

void StopMu::fillCalorimetry(int pl, std::vector<float> dqdx, std::vector<float> rr)//, std::vector<TVector3> xyz)
{
  if (pl==0)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_u.push_back((double)dqdx[n]);
    _rr_u.push_back(  (double)rr[n]  );
    _tmean.CalcTruncMeanProfile(_rr_u, _dqdx_u, _dqdx_tm_u);
    //_x_position_u.push_back((double)(xyz[n].X()));
    //_y_position_u.push_back((double)(xyz[n].Y()));
    //_z_position_u.push_back((double)(xyz[n].Z()));
    }
  }
  else if (pl==1)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_v.push_back((double)dqdx[n]);
    _rr_v.push_back(  (double)rr[n]  );
    _tmean.CalcTruncMeanProfile(_rr_v, _dqdx_v, _dqdx_tm_v);
    //_x_position_v.push_back((double)(xyz[n].X()));
    //_y_position_v.push_back((double)(xyz[n].Y()));
    //_z_position_v.push_back((double)(xyz[n].Z()));
    }
  }
  else if (pl==2)
  {
    for (size_t n=0; n < dqdx.size(); n++)
    {
    _dqdx_y.push_back((double)dqdx[n]);
    _rr_y.push_back(  (double)rr[n]  );
    _tmean.CalcTruncMeanProfile(_rr_y, _dqdx_y, _dqdx_tm_y);
    //_x_position_y.push_back((double)(xyz[n].X()));
    //_y_position_y.push_back((double)(xyz[n].Y()));
    //_z_position_y.push_back((double)(xyz[n].Z()));
    }
  }
}

bool StopMu::insideTPCvolume(double x, double y, double z)
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

double StopMu::getPitch(const TVector3 &direction, const int &pl)
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
  double pitch = minWireSpacing / cos;
  return pitch;
}

void StopMu::shiftTruePosition(double true_point[3], double true_time, double true_point_shifted[3])
{
  true_point_shifted[0] = true_point[0];
  true_point_shifted[1] = true_point[1];
  true_point_shifted[2] = true_point[2];

  ::detinfo::DetectorProperties const* _detector_properties;
  ::detinfo::DetectorClocks const* _detector_clocks;

  _detector_properties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _detector_clocks = lar::providerFrom<detinfo::DetectorClocksService>();

  double g4Ticks = _detector_clocks->TPCG4Time2Tick(true_time)
                       + _detector_properties->GetXTicksOffset(0,0,0)
                       - _detector_properties->TriggerOffset();
  double xOffset = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  true_point_shifted[0] += xOffset;

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  // auto offset = SCE->GetPosOffsets(true_point[0], true_point[1], true_point[2]);
  // if (offset.size() == 3)
  // {
  //   true_point_shifted[0] -= offset[0];
  //   true_point_shifted[1] += offset[1];
  //   true_point_shifted[2] += offset[2];
  // }
  auto offset = SCE->GetPosOffsets(geo::Point_t(true_point[0], true_point[1], true_point[2]));
  true_point_shifted[0] -= offset.X();
  true_point_shifted[1] += offset.Y();
  true_point_shifted[2] += offset.Z();
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
