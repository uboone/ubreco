
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

typedef std::vector<int>        vint_t;
typedef std::vector<bool>       vbool_t;
typedef std::vector<float>      vfloat_t;
typedef std::set<int>           si_t;
typedef std::map<int,int>       mii_t;
typedef std::map<int,float>     mif_t;
typedef std::map<float,float>   mff_t;

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
    TVector3 position;
  };
  
  struct TrueBlip {
    int       ID            = -9;
    bool      isValid       = false;
    int       TPC           = -9;
    float     Energy        = 0;
    int       DepElectrons  = 0;
    float     NumElectrons  = 0; // (post-drift)
    float     Length        = 0;
    float     Time          = -999e9;
    int       LeadG4ID      = -9;
    int       LeadG4Index   = -9;
    int       LeadG4PDG     = -9;
    float     LeadEnergy    = -9;
    mif_t     G4EnergyMap;
    mif_t     G4PDGMap;
    TVector3  Position;
  };

  struct HitInfo {
    art::Ptr<recob::Hit> Hit;
    int   hitid         = -9;
    int   plane         = -9;
    int   tpc           = -9;
    int   wire          = -9;
    int   chan          = -9;
    int   trkid         = -9;
    int   shwrid        = -9;
    int   clustid       = -9;
    int   blipid        = -9;
    bool  ismatch       = false;
    int   neighbors     = -9;
    float charge        = -999;     // [e-]
    float driftTime     = -999999;  // [tick]
    int   g4id          = -9;
    int   g4pdg         = -999;
    float g4frac        = -99;      
    float g4charge      = -999;     // [e-]
    float g4energy      = -999;     // [MeV]
  };
  
  struct HitClust {
    //art::Ptr<recob::Hit>  LeadHit;
    //int     LeadHitID       = -999;
    float   LeadHitCharge   = -999;
    //float   LeadHitTime     = -999;
    int     ID              = -9;
    //int     LeadHitWire     = -999;
    //int     LeadHitChan     = -999;
    float   CentHitTime     = -999;
    int     CentHitChan     = -999;
    int     CentHitWire     = -999;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     TPC             = -9;
    int     Plane           = -9;
    int     NHits           = -9;
    int     NWires          = -9;
    float   ADCs            = -999;
    float   Amplitude       = -999;
    float   Charge          = -999;
    float   Time            = -999;
    float   TimeErr         = -999;
    float   StartHitTime    = -999;
    float   EndHitTime      = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    float   Timespan        = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     BlipID          = -9;
    int     EdepID          = -9;
    int     G4ID            = -9;
    float   G4PDG           = -999;
    //float   G4Energy        = -9;
    int     TrkID           = -9;
    si_t    HitIDs;
    si_t    Wires;
    mii_t   TrkIDs;
    
    std::map<int,TVector3> IntersectLocations;

    vfloat_t  HitTimes;
    vint_t    HitChans;
    vint_t    HitWires;
    //vint_t    HitTrkIDs;
  };

  struct Blip {
    
    int       ID              = -9;         // Blip ID / index
    bool      isValid         = false;      // Blip passes basic checks
    int       TPC             = -9;         // TPC
    int       NPlanes         = -9;         // Num. matched planes
    int       MaxWireSpan     = -9;         // Maximum span of wires on any plane cluster
    float     Energy          = -999;       // Energy (const dE/dx, fcl-configurable)
    float     EnergyESTAR     = -999;       // Energy (ESTAR method from ArgoNeuT)
    float     Time            = -999;       // Drift time [ticks]
    float     ProxTrkDist     = -9;         // Distance to cloest track
    int       ProxTrkID       = -9;         // ID of closest track
    bool      inCylinder      = false;      // Is it in a cone/cylinder region? 
    
    TVector3  Position;                     // 3D position TVector3
    float     X               = -999;       // Reconstructed X [cm]
    float     Y               = -999;       // Reconstructed Y [cm]
    float     Z               = -999;       // Reconstructed Z [cm]
    float     SigmaYZ         = -9.;        // Uncertainty in YZ intersect [cm]
    float     dX              = -9;         // Equivalent length along drift direction [cm] 
    float     dYZ             = -9;         // Approximate length scale in YZ space [cm]
    float     TrkID           = -9;         // Track ID (if hits were 3D-tracked)
    float     Length          = -9;         // Length of 3D track (if applicable)

    // Plane/cluster-specific information
    blip::HitClust clusters[kNplanes];
    float   Match_dT[kNplanes];
    float   Match_dTfrac[kNplanes];
    float   Match_overlap[kNplanes];
    float   Match_score[kNplanes];
    
    // Truth-matched energy deposition
    blip::TrueBlip truth;
    
    // Initialize the blip
    Blip() {
      for(size_t i=0; i<kNplanes; i++){
        Match_dT[i]     = -999;
        Match_dTfrac[i] = -999;
        Match_overlap[i]= -9;
        Match_score[i]  = -9;
      }
    }

  };
  
}


