
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

typedef std::vector<int>        vint_t;
typedef std::vector<bool>       vbool_t;
typedef std::vector<float>      vfloat_t;
typedef std::set<int>           si_t;
typedef std::map<int,float>     mif_t;
typedef std::map<int,TVector3>  mitv3_t;

const int kNplanes  = 3;  

namespace blipobj {
  
  //###################################################
  //  Data structures
  //###################################################

  struct ParticleInfo {
    simb::MCParticle particle;
    int   trackId           = -9;
    int   index             = -9;
    int   isPrimary         = -9;
    int   numTrajPts          = -9;
    double depEnergy         = -9;
    int   depElectrons        = -9;
    double numElectrons      = -9;
    double mass              = -9;
    double E                 = -9;
    double endE              = -9;
    double KE                = -9;
    double endKE             = -9;
    double P                 = -9; 
    double Px                = -9; 
    double Py                = -9; 
    double Pz                = -9; 
    double pathLength        = -9;
    double time              = -9;
    double endtime           = -9;
    TVector3 startPoint;      
    TVector3 endPoint;
    TVector3 position;
  };

  

  // True energy depositions
  struct TrueBlip {
    int       ID            = -9;     // unique blip ID
    int       TPC           = -9;     // TPC ID
    float     Time          = -999e9; // time [us]
    float     Energy        = -9;      // energy dep [MeV]
    int       DepElectrons  = -9;      // deposited electrons
    int       NumElectrons  = -9;      // electrons reaching wires
    float     DriftTime     = -9;     // drift time [us]
    int       LeadG4ID      = -9;     // lead G4 track ID
    int       LeadG4Index   = -9;     // lead G4 track index
    int       LeadG4PDG     = -9;     // lead G4 PDG
    float     LeadCharge    = -9;     // lead G4 charge dep
    TVector3  Position;               // XYZ position
    bool      AllChansGood  = true;    
    mif_t     G4ChargeMap;          
  };

  struct HitInfo {
    int   hitid         = -9;
    int   plane         = -9;
    int   tpc           = -9;
    int   wire          = -9;
    int   chan          = -9;
    float amp           = -9;
    float rms           = -9;
    int   trkid         = -9;
    int   shwrid        = -9;
    int   clustid       = -9;
    int   blipid        = -9;
    bool  ismatch       = false;
    bool  touchTrk      = false;
    int   touchTrkID    = -9;
    float integralADC    = -999;     // [ADCs] from integral
    float sigmaintegral = -999;
    float sumADC        = -999;     // [ADCs] from sum 
    float charge        = -999;     // [e-]
    float peakTime      = -999999;
    float driftTime     = -999999;  // [tick]
    float gof           = -9;
    int   g4trkid       = -9;
    int   g4pdg         = -999;
    int   g4charge      = -999;     // [e-]
    float g4frac        = -99;      
    float g4energy      = -999;     // [MeV]
  };

  struct HitClust {
    int     ID              = -9;
    bool    isValid         = false;
    int     CenterChan      = -999;
    int     CenterWire      = -999;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     TouchTrkID      = -9;
    int     DeadWireSep     = 99;
    int     TPC             = -9;
    int     Plane           = -9;
    int     NHits           = -9;
    int     NWires          = -9;
    int     NWiresNoisy     = -9;
    int     NWiresBad       = -9;
    float   SummedADC       = -999;
    float   Amplitude       = -999;
    float   Charge          = -999;
    float   SigmaCharge     = -999;
    float   Time            = -999;
    float   RMS             = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    float   Timespan        = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     NPulseTrainHits = -9;
    int     BlipID          = -9;
    int     EdepID          = -9;
    si_t    HitIDs;
    si_t    Wires;
    si_t    Chans;
    si_t    G4IDs;
    mitv3_t IntPts;
  };
  
  struct Blip {
    
    int       ID              = -9;         // Blip ID / index
    bool      isValid         = false;      // Blip passes basic checks
    int       TPC             = -9;         // TPC
    int       NPlanes         = -9;         // Num. matched planes
    int       MaxWireSpan     = -9;         // Maximum span of wires on any plane cluster
    float     Charge          = -9;         // Charge on calorimetry plane
    float     ChargeCorr      = -9;         // Charge on calorimetry plane (lifetime corrected)
    float     Energy          = -999;       // Energy (const dE/dx, fcl-configurable)
    float     EnergyCorr      = -999;       // Energy following SCE + lifetime correction
    float     Time            = -999;       // Drift time [ticks]
    TVector3  Position;                     // 3D position
    TVector3  PositionSCE;                  // 3D position following SCE spatial correction
    float     SigmaYZ         = -9.;        // Uncertainty in YZ intersect [cm]
    float     dX              = -9;         // Equivalent length along drift direction [cm] 
    float     dYZ             = -9;         // Approximate length scale in YZ space [cm]
                                            // (also referred to as 'dW')
    float     ProxTrkDist     = -9;         // Distance to closest track
    int       ProxTrkID       = -9;         // ID of closest track
    bool      inCylinder      = false;      // Is it in a cone/cylinder region? 
    int       TouchTrkID      = -9;         // Track ID of track that is touched

    // Plane/cluster-specific information
    blipobj::HitClust clusters[kNplanes];

    // Truth-matched energy deposition
    blipobj::TrueBlip truth;
    
    // Prototype getter functions
    double X() { return Position.X(); }
    double Y() { return Position.Y(); }
    double Z() { return Position.Z(); }
    
  };
  
}

// Prototype data product that can serve as a "lighter"
// blip object to store in artROOT events if we run into
// file size issues using the default object above.
/*
namespace blpobj {
  
  struct Blip {
    int       ID            = -9;   // Blip ID / index
    int       TPC           = -9;   // TPC
    int       NPlanes       = -9;   // Num. matched planes
    TVector3  Location;             // Reconstructed XYZ
    TVector3  LocationSCE;          // Reconstructed XYZ w/SCE offset corrections applied
    float     Charge        = -9;   // Charge on calorimetry plane
    float     ChargeCorr    = -9;   // Charge w/SCE and lifetime corrections
    float     Energy        = -999; // Energy (const dE/dx, fcl-configurable)
    float     EnergyCorr    = -999;
    float     ProxTrkDist   = -9;   // Distance to closest track
    int       ProxTrkID     = -9;   // ID of closest track
    int       TouchTrkID    = -9;   // Track ID of track that is touched
 
    blipobj::HitClust hitcluster2D[kNplanes];
    
    float     trueEnergy    = 0;    // Truth-matched energy dep [MeV]
    int       trueCharge    = 0;    // Truth-matched deposited electrons
    int       trueG4ID      = -9;   // Truth-matched lead G4 track ID
    int       truePDG       = -9;   // Truth-matched lead G4 PDG

    // To add:
    //  x nwires per plane
    //  x location w/SCE corrections
    //  x energy w/lifetime+SCE corrections
    //  x flag for dead wire adjacency --> possible thru 'cluster'
    //  x flag for "noisy" wire adjacency --> possible thru 'cluster'
  };
}
*/

