///////////////////////////////////////////////////////
// BlipUtils.h
//
// Helper functions and algs. Based on RecoUtils and
// functions in AnalysisTree.
//
//////////////////////////////////////////////////////
#ifndef BLIPUTIL_H_SEEN
#define BLIPUTIL_H_SEEN

// framework
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"

// c++
#include <vector>
#include <map>

typedef std::vector<int> vi_t;
typedef std::set<int> si_t;
typedef std::map<int,float> mif_t;
typedef std::map<float,float> mff_t;
typedef std::vector<art::Ptr<sim::SimEnergyDeposit>> SEDVec_t;
typedef std::vector<simb::MCParticle const*> vmcp_t;
typedef std::vector<anab::BackTrackerHitMatchingData const*> vbt_t;

const int kNplanes  = 3;  

namespace blip {
  
  //###################################################
  //  Data structures
  //###################################################

  struct ParticleInfo {
    simb::MCParticle particle;
    int   trackId           = -9;
    int   index             = -9;
    int   isPrimary         = -9;
    float depEnergy         = -9;
    int   depElectrons      = -9;
    float numElectrons      = -9;
    float mass              = -9;
    float E                 = -9;
    float endE              = -9;
    float KE                = -9;
    float endKE             = -9;
    float P                 = -9; 
    float Px                = -9; 
    float Py                = -9; 
    float Pz                = -9; 
    float pathLength        = -9;
    float time              = -9;
    float endtime           = -9;
    TVector3 startPoint;      
    TVector3 endPoint;      
  };
  
  struct HitInfo {
    art::Ptr<recob::Hit> hit;
    int   hitid         = -9;
    int   plane         = -9;
    int   tpc           = -9;
    int   wire          = -9;
    int   trkid         = -9;
    int   shwrid        = -9;
    bool  ismatch       = false;
    bool  isneartrk     = false;
    float trkdist2D     = -9;
    int   blipid        = -9;
    int   clustid       = -9;
    int   g4id          = -9;
    float g4frac        = -99;
    float g4energy      = -999;
    float g4charge      = -999;
    float qcoll         = -999;
    float driftTicks    = -999999;
    bool  rmsCut        = false;
    bool  ratioCut      = false;
    si_t  g4ids;
    TVector2 point2D;
  };

  struct TrueBlip {
    int       ID;
    bool      isValid       = false;
    int       TPC           = -9;
    int       LeadG4ID      = -9;
    int       LeadG4Index   = -9;
    int       LeadG4PDG     = -9;
    float     LeadEnergy    = -9;
    float     Energy        = 0;
    int       DepElectrons  = 0;
    float     NumElectrons  = 0; // (post-drift)
    float     Length        = 0;
    float     Time          = -999e9;
    TVector3  Position;
    TVector3  StartPoint;
    TVector3  EndPoint;
    vi_t      G4IDs;
    vi_t      PDGs;
  };

  
  struct HitClust {
    art::Ptr<recob::Hit> LeadHit;
    bool    isValid         = false;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     LeadHitID       = -999;
    float   LeadHitCharge   = -999;
    float   LeadHitTime     = -999;
    int     G4ID            = -9;
    int     TPC             = -9;
    int     Plane           = -9;
    float   ADCs            = -999;
    float   Charge          = -999;
    float   Time            = -999;
    float   WeightedTime    = -999;
    float   TimeErr         = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     ID              = -9;
    int     BlipID          = -9;
    int     EdepID          = -9;
    si_t    HitIDs;
    si_t    G4IDs;
    si_t    Wires;
    mif_t   mapWireCharge;
    mff_t   mapTimeCharge;
  };

  struct Blip {
    int       ID              = -9;
    bool      isValid         = false;
    int       TPC             = -9;
    int       NPlanes         = -9;                     // Num. matched planes
    int       Planes[3]       = {false, false, false};  // Match status for each plane
    float     Charge[3]       = {-999, -999, -999};     // Charge (e-) on each plane
    float     Energy          = -999;                   // Energy (const dE/dx = 2 MeV/cm)
    float     EnergyESTAR     = -999;                   // Energy (ESTAR method from ArgoNeuT)
    float     DriftTime       = -999;                   // Drift time (ticks)
    float     MaxIntersectDiff= -9;                     // Max difference between wire intersection 
                                                        // points (only valid for >=3 planes)
    TVector3  Position;         // 3D position vector
    std::set<int> ClustIDs;     // Associated blip::HitClusts
    std::set<int> HitIDs;       // Associated recob::Hits
  };
}


namespace BlipUtils {
  
  // Tick period
  float   fTickPeriod;
  
  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  void      FillParticleInfo(simb::MCParticle const&, blip::ParticleInfo&, SEDVec_t&);
  void      CalcTotalDep(float&,int&,float&, SEDVec_t&);
  void      CalcPartDrift(int,float&);
  void      MakeTrueBlips(std::vector<blip::ParticleInfo> const&, std::vector<blip::TrueBlip>&);
  void      GrowTrueBlip(blip::ParticleInfo const&, blip::TrueBlip&);
  void      MergeTrueBlips(std::vector<blip::TrueBlip>&, float);
  void      GrowHitClust(blip::HitInfo const&, blip::HitClust&);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);
  bool      DoHitClustsOverlap(blip::HitClust const&, blip::HitClust const&);
  bool      DoHitClustsOverlap(blip::HitClust const&,float,float);
  bool      DoChannelsIntersect(int,int);
  bool      DoHitClustsMatch(blip::HitClust const&, blip::HitClust const&,float);
  float     ModBoxRecomb(float,float);
  float     ConvertTicksToX(float, int, int, int);
  blip::HitClust  MakeHitClust(blip::HitInfo const&);
  blip::HitClust  MergeHitClusts(blip::HitClust&, blip::HitClust&);
  blip::Blip      MakeBlip(std::vector<blip::HitClust> const&);

  //###################################################
  // General functions 
  //###################################################
  void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  bool    G4IdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  bool    DoesHitHaveSimChannel( art::Ptr<recob::Hit> const&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);
  double  DistToBoundary(const recob::Track::Point_t&);
  double  DistToLine(TVector3&, TVector3&, TVector3&);
  double  DistToLine2D(TVector2&, TVector2&, TVector2&);
  void    GetGeoBoundaries(double&,double&,double&,double&,double&,double&);
  bool    IsPointInAV(float,float,float);
  bool    IsPointInAV(TVector3&);

}

#endif
