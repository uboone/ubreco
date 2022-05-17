
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

typedef std::vector<int> vi_t;
typedef std::set<int> si_t;
typedef std::map<int,float> mif_t;
typedef std::map<float,float> mff_t;

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
    int       NPlanes         = -9;         // Num. matched planes
    bool      Planes[kNplanes]= {false};    // Match status for each plane
    float     SumADC[kNplanes]= {-999};
    float     Charge[kNplanes]= {-999};     // Charge (e-) on each plane
    float     Energy          = -999;       // Energy (const dE/dx = 2 MeV/cm)
    float     EnergyESTAR     = -999;       // Energy (ESTAR method from ArgoNeuT)
    float     DriftTime       = -999;       // Drift time (ticks)
    float     MaxIntersectDiff= -9;         // Max diff. btwn wire intersections
    float     trkdist         = -9;         // Distance to cloest track
    int       trkid           = -9;         // ID of closest track
    bool      inCylinder      = false;      // Whether this blip is within a 
                                            //   track's cone/cylinder region
    // 3D position
    TVector3  Position;                     // 3D position TVector3
    float     x               = -999;       // Reconstructed X [cm]
    float     y               = -999;       // Reconstructed Y [cm]
    float     z               = -999;       // Reconstructed Z [cm]
    
    // Plane-specific cluster info
    int       clustID[kNplanes]     = {-9};     // clust IDs per plane
    float     matchscore[kNplanes]  ={-9};
    int       nhits[kNplanes]       = {-9};     // NHits per cluster on each plane
    int       nwires[kNplanes]      = {-9};     // NWires per cluster on each plane
    float     timespan[kNplanes]    = {-9};     // cluster timespan [ticks] per plane
    std::set<int> ClustIDs;                     // IDs of associated blip::HitClusts
    std::set<int> HitIDs;                       // IDs of associated recob::Hits

    // Truth-matching information
    TrueBlip  truth;

  };
}


