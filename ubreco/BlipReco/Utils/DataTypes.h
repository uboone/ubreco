
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

typedef std::vector<int> vint_t;
typedef std::vector<bool> vbool_t;
typedef std::vector<float> vfloat_t;
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
    art::Ptr<recob::Hit> Hit;
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
    int   g4pdg         = -999;
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
    vint_t      G4IDs;
    vint_t      PDGs;
  };

  
  struct HitClust {
    art::Ptr<recob::Hit> LeadHit;
    int     ID              = -9;
    bool    isValid         = false;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     LeadHitID       = -999;
    float   LeadHitCharge   = -999;
    float   LeadHitTime     = -999;
    int     LeadHitWire     = -999;
    int     TPC             = -9;
    int     Plane           = -9;
    float   ADCs            = -999;
    float   Charge          = -999;
    float   Time            = -999;
    float   TimeErr         = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     BlipID          = -9;
    int     EdepID          = -9;
    int     G4ID            = -9;
    float   G4Energy        = -9;
    float   G4PDG           = -999;
    si_t    HitIDs;
    si_t    G4IDs;
    si_t    Wires;
    //mif_t   mapWireCharge;
    //mff_t   mapTimeCharge;
  };

  struct Blip {

    int       ID              = -9;         // Blip ID / index
    bool      isValid         = false;      // Blip passes basic checks
    int       TPC             = -9;         // TPC
    int       NPlanes         = -9;         // Num. matched planes
    float     Energy          = -999;       // Energy (const dE/dx = 2 MeV/cm)
    float     EnergyESTAR     = -999;       // Energy (ESTAR method from ArgoNeuT)
    float     DriftTime       = -999;       // Drift time (ticks)
    float     MaxIntersectDiff= -9;         // Max diff. btwn wire intersections
    float     TrkDist         = -9;         // Distance to cloest track
    int       TrkID           = -9;         // ID of closest track
    bool      inCylinder      = false;      // Is it in a cone/cylinder region? 
    TVector3  Position;                     // 3D position TVector3
    float     X               = -999;       // Reconstructed X [cm]
    float     Y               = -999;       // Reconstructed Y [cm]
    float     Z               = -999;       // Reconstructed Z [cm]
    
    // Sets of cluster and hit IDs
    std::set<int> ClustID_set;             // IDs of associated blip::HitClusts
    std::set<int> HitID_set;               // IDs of associated recob::Hits
    
    // Plane-specific information
    bool   Planes[kNplanes];  // Match status for each plane
    float   SumADC[kNplanes];   // Sum ADC per plane
    float   Charge[kNplanes];   // Charge (e-) on each plane
    int     ClustID[kNplanes];  // clust IDs per plane
    int     NHits[kNplanes];    // Number hits per plane
    int     NWires[kNplanes];   // Number wires per plane
    float   Timespan[kNplanes]; // Cluster timespan [ticks] per plane

    // Truth-matching information
    TrueBlip  truth;
    
    //----------------------------
    // Initialize the blip
    //----------------------------
    Blip() {
      for(size_t i=0; i<kNplanes; i++){
        Planes[i] = false;
        ClustID[i] = -9;
      }
    }

  };
  
}


