/////////////////////////////////////////////////////////////////////
// 
//  BlipRecoAlg
//
//  To include in your module, add the header file at the top:
//    #include "ubreco/BlipReco/Alg/BlipRecoAlg.h
//
//  then declare a private alg object in your class definition:
//    BlipRecoAlg fBlipAlg;
//    
//  W. Foreman, May 2022
//  wforeman @ iit.edu
//
/////////////////////////////////////////////////////////////////////
#ifndef BLIPRECOALG_H
#define BLIPRECOALG_H

// framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
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
#include "larreco/Calorimetry/CalorimetryAlg.h"

// Microboone includes
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
//#include "ubevt/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"
//#include "ubevt/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// Blip-specific utils
#include "ubreco/BlipReco/Utils/BlipUtils.h"

// ROOT stuff
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

// c++
#include <vector>
#include <iostream>
#include <memory>
#include <math.h>
#include <limits>


namespace blip {

  //--------------------------------------------
  class BlipRecoAlg{
   public:
    
    //Constructor/destructor
    BlipRecoAlg( fhicl::ParameterSet const& pset );
    BlipRecoAlg();
    ~BlipRecoAlg();
  
    void    reconfigure(fhicl::ParameterSet const& pset );
    void    RunBlipReco(const art::Event& evt);
    void    PrintConfig();
    
    // TO-DO: make these private and create getters instead
    std::vector<blip::HitInfo>      hitinfo;
    std::vector<blip::HitClust>     hitclust;
    std::vector<blip::Blip>         blips;  
    std::vector<blip::TrueBlip>     trueblips;
    std::vector<blip::ParticleInfo> pinfo;
    
    calo::CalorimetryAlg*   fCaloAlg;
    float   ModBoxRecomb(float,float);
    float   dQdx_to_dEdx(float,float);

   private:
   
    const detinfo::DetectorProperties* detProp;

    // --- FCL configs ---
    std::string         fHitProducer;
    std::string         fTrkProducer;
    std::string         fGeantProducer;
    std::string         fSimDepProducer;
    std::string         fSimChanProducer;
    float               fSimGainFactor;
    bool                fDebugMode;
    float               fTrueBlipMergeDist;
    bool                fDoHitFiltering;
    float               fMaxHitTrkLength;
    float               fMaxHitAmp;
    std::vector<float>  fMinHitAmp;
    std::vector<float>  fMinHitRMS;
    std::vector<float>  fMaxHitRMS;
    std::vector<float>  fMinHitGOF;
    std::vector<float>  fMaxHitGOF;
    std::vector<float>  fMinHitRatio;
    std::vector<float>  fMaxHitRatio;
    int                 fMaxHitMult;
    float               fHitClustWidthFact;
    int                 fHitClustWireRange;
    float               fMatchQDiffLimit;
    float               fMatchMaxQRatio;
    //std::vector<float>  fTimeOffsets;
    float               fMatchMinOverlap;
    float               fMatchSigmaFact;
    float               fMatchMaxTicks;
    int                 fMaxWiresInCluster;
    float               fMinClusterCharge;
    float               fMaxClusterSpan;
    int                 fMinMatchedPlanes;
    bool                fPickyBlips;
    bool                fApplyTrkCylinderCut;
    float               fCylinderRadius; 
 
    // --- Calorimetry configs ---
    int                 fCaloPlane;
    float               fCalodEdx;
    bool                fLifetimeCorr;
    bool                fSCECorr;
    bool                fYZUniformityCorr;
    float               fModBoxA;
    float               fModBoxB;
    float               fWion;
  
    // --- Histograms ---
    TH1D*   h_clust_nwires;
    TH1D*   h_clust_timespan;
    TH1D*   h_clust_overlap[kNplanes];
    TH1D*   h_clust_dt[kNplanes];
    TH1D*   h_clust_dtfrac[kNplanes];
    TH2D*   h_clust_q[kNplanes]; 
    TH2D*   h_clust_q_cut[kNplanes]; 
    TH1D*   h_clust_picky_overlap[kNplanes];
    TH1D*   h_clust_picky_dt[kNplanes];
    TH1D*   h_clust_picky_dtfrac[kNplanes];
    TH2D*   h_clust_picky_q[kNplanes]; 
    TH1D*   h_nmatches[kNplanes];

    TH1D*   h_recomb;
    TH1D*   h_efield;

  };

}

#endif
