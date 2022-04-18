//#####################################################################
//###  BlipAna analyzer module
//###
//###  Contains algorithms for reconstructing isolated, MeV-scale energy
//###  depositions in the TPC, called "blips." A TTree is made for offline
//###  analysis and plot-making. Algs will eventually be migrated into 
//###  dedicated alg/tool classes as appropriate.
//###
//###  Author: Will Foreman (wforeman_at_iit.edu)
//###  Date:   Sept 2021
//#####################################################################
#ifndef BLIPANA_H
#define BLIPANA_H

// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
//#include "larcorealg/Geometry/WireGeo.h"
//#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "cetlib/search_path.h"

//#include "larevt/CalibrationDBI/Interface/ElectronLifetimeService.h"
//#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"

// BlipReco-specific utility includes
#include "ubreco/BlipReco/Utils/BlipUtils.h"

// C++ includes
#include <cstring>
#include <vector>
#include <map>
#include <utility>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
//#include "TTimeStamp.h"
#include "TF1.h"


// Set global constants and max array sizes
const int kMaxHits  = 50000;
const int kMaxTrks  = 10000;
const int kMaxShwrs = 1000;
const int kMaxBlips = 10000;
const int kMaxG4    = 100000;
const int kMaxEDeps = 100000;
const int kNplanes  = 3;  
   
class BlipAna;
  
//###################################################
//  Data storage structure
//###################################################
class BlipAnaTreeDataStruct 
{
  public:

  // --- Main TTree object ---
  TTree* tree;

  // --- Configurations and switches ---
  std::string treeName = "anatree";
  bool  saveTruthInfo  = true;
  bool  saveTrkInfo    = true;
  bool  saveHitInfo    = true;
  bool  saveClustInfo  = true;

  // --- Event information ---   
  int           event;                    // event number
  int           run;                      // run number
  unsigned int  timestamp;                // unix time of event
  float         lifetime;                 // electron lifetime
  
  // --- G4 information ---
  float total_depEnergy;          // total deposited energy in AV
  float total_numElectrons;       // total electrons reaching anode wires
  //float gamma_depEnergy;          // total gamma-induced energy deposited
  //float gamma_numElectrons;       // total electrons from gamma-induced depositions
  short nparticles;               // number of G4 particles
  bool  isPrimary[kMaxG4];        // is primary particle
  short trackID[kMaxG4];          // G4 track ID
  int   pdg[kMaxG4];              // PDG
  short nDaughters[kMaxG4];       // number of daughters
  short mother[kMaxG4];           // mother particle
  float E[kMaxG4];                // initial energy (MeV)
  float KE[kMaxG4];               // initial kinetic energy (MeV)
  float endE[kMaxG4];             // final energy (MeV)
  float endKE[kMaxG4];             // final energy (MeV)
  float mass[kMaxG4];             // mass (MeV)
  float P[kMaxG4];                // momentum (MeV)
  float Px[kMaxG4];               // momentum x (MeV)
  float Py[kMaxG4];               // momentum y (MeV)
  float Pz[kMaxG4];               // momentum z (MeV)
  float startPointx[kMaxG4];      // starting x (cm)
  float startPointy[kMaxG4];      // starting y (cm)
  float startPointz[kMaxG4];      // starting y (cm)
  float endPointx[kMaxG4];        // ending x (cm)
  float endPointy[kMaxG4];        // ending y (cm)
  float endPointz[kMaxG4];        // ending y (cm)
  float startT[kMaxG4];           // starting time (us)
  float endT[kMaxG4];             // ending time (us)
  float pathlen[kMaxG4];          // path length (cm)
  float numElectrons[kMaxG4];     // electrons reaching anode wires
  float depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  std::vector<std::string> process;// process name

  // --- True energy deposit info (derived) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 track ID
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_clustid[kMaxEDeps];  // hitclust ID
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
  //float edep_energyESTAR[kMaxEDeps];   // total energy deposited (MeV)
  float edep_charge[kMaxEDeps];   // total electrons reaching anode wires
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)
  float edep_ds[kMaxEDeps];       // extent (cm)

  // --- Hit information ---
  int	  nhits;                    // number of hits
  int	  hit_tpc[kMaxHits];        // tpc number
  int	  hit_plane[kMaxHits];      // plane number
  int	  hit_wire[kMaxHits];       // wire number
  int	  hit_channel[kMaxHits];    // channel ID
  float	hit_peakT[kMaxHits];      // raw peak time (tick)
  float	hit_time[kMaxHits];       // corrected peak time (tick)
  float hit_rms[kMaxHits];        // shape RMS
  float	hit_amp[kMaxHits];         // amplitude
  float	hit_area[kMaxHits];       // charge (area) in ADC units
  float	hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
  int   hit_isgood[kMaxHits];     // does hit pass shape-based cuts?
  int   hit_isreal[kMaxHits];     // is this hit real?
  int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
  float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
  float hit_g4energy[kMaxHits];   // true energy
  float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
  int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
  int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
  float hit_gof[kMaxHits];   // goodness of fit (default -1)


  // --- Track hits ---
  float ave_trkhit_amp;            // average pulse height for hits in tracks

  // --- Track information ---
  int   ntrks;                    // number tracks
  short trk_id[kMaxTrks];         // trackID
  short trk_npts[kMaxTrks];       // number 3D trajectory points
  float trk_length[kMaxTrks];     // track length [cm]
  float trk_startx[kMaxTrks];     // starting X coordinate
  float trk_starty[kMaxTrks];     // starting Y coordinate
  float trk_startz[kMaxTrks];     // starting Z coordinate
  float trk_startd[kMaxTrks];     // starting distance to boundary
  float trk_endx[kMaxTrks];       // ending X coordinate
  float trk_endy[kMaxTrks];       // ending Y coordinate
  float trk_endz[kMaxTrks];       // ending Z coordinate
  float trk_endd[kMaxTrks];       // ending distance to boundary
  //int   trk_origin[kMaxTrks];     // 0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
  //int   trk_g4id[kMaxTrks];       // G4 TrackID of leading contributing particle
  //int   trk_g4pdg[kMaxTrks];      // G4 PDG
  //float trk_purity[kMaxTrks];     // track hit purity
  //float trk_pitch[kMaxTrks];      // track pitch
  //float trk_ke[kMaxTrks];         // track kinetic energy
  //int   trk_cosmictag[kMaxTrks];  // cosmic tagg
  //float trk_cosmicscore[kMaxTrks];// cosmic score
  //int   trk_pidpdg[kMaxTrks];     // particle PID Pdg
  //float trk_pidchi[kMaxTrks];     // chisquared for PID
  //int   trk_bestplane[kMaxTrks];  // plane w/ most hits
  
  // --- Shower information ---
  /*
  int   nshwrs;                     // number showers
  int   shwr_id[kMaxShwrs];         // shower ID
  float shwr_dirx[kMaxShwrs];       // direction X
  float shwr_diry[kMaxShwrs];       // direction Y
  float shwr_dirz[kMaxShwrs];       // direction Z
  float shwr_startx[kMaxShwrs];     // starting X coordinate
  float shwr_starty[kMaxShwrs];     // starting X coordinate
  float shwr_startz[kMaxShwrs];     // starting X coordinate
  float shwr_length[kMaxShwrs];     // shower length
  float shwr_openangle[kMaxShwrs];  // opening angle
  */

  // --- Hit cluster information ---
  int   nclusts;                      // total clusters made
  short clust_tpc[kMaxHits];          // cluster TPC ID
  short clust_plane[kMaxHits];        // cluster plane 
  short clust_wire[kMaxHits];         // cluster wire (lead hit wire)
  short clust_startwire[kMaxHits];
  short clust_endwire[kMaxHits];
  short clust_chan[kMaxHits];         // cluster channel (lead hit wire)
  int   clust_id[kMaxHits];           // cluster ID (index)
  short clust_nwires[kMaxHits];       // number of wires in this cluster
  short clust_nhits[kMaxHits];        // number of hits
  int   clust_lhit_id[kMaxHits];      // lead hit ID (index for hit_X[i] branches)
  float clust_lhit_amp[kMaxHits];   // lead hit peak amplitude [ADC]
  float clust_lhit_rms[kMaxHits];   // lead hit RMS [ADC]
  float clust_lhit_time[kMaxHits];  // lead hit time-tick, corrected for offsets
  float clust_lhit_peakT[kMaxHits];   // lead hit time-tick, uncorrected (hit->PeakT)
  //float clust_lhit_gof[kMaxHits];     // lead hit goodness-of-fit
  bool  clust_lhit_isfit[kMaxHits];   // is there a valid goodness of fit for lead hit?
  bool  clust_ismatch[kMaxHits];      // was this cluster plane-matched?
  float clust_charge[kMaxHits];       // cluster charge at anode [e-]
  float clust_time[kMaxHits];         // cluster time-tick
  float clust_time_w[kMaxHits];       // cluster time-tick (charge-weighted)
  float clust_time_err[kMaxHits];     // cluster time uncertainty
  float clust_startTime[kMaxHits];    // cluster start tick
  float clust_endTime[kMaxHits];      // cluster end tick
  float clust_timespan[kMaxHits];     // cluster timespan
  float clust_g4energy[kMaxHits];     // true cluster energy from G4
  float clust_g4charge[kMaxHits];     // true cluster charge at anode
  int   clust_g4id[kMaxHits];         // true MCParticle ID (index for particle branches)
  int   clust_blipid[kMaxHits];       // blip ID for this clusteer (if it was made into one)
  int   clust_edepid[kMaxHits];       // true energy dep ID

  // --- 3D Blip information ---
  float total_blip_energy;            // total summed blip energy in event [MeV]
  int   nblips;                       // number of blips in event
  short blip_id[kMaxBlips];           // blip ID / index
  short blip_tpc[kMaxBlips];          // blip TPC
  short blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_maxdiff[kMaxBlips];      // difference in wire intersection points
  float blip_charge[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip energy [MeV]
  //float blip_energyESTAR[kMaxBlips];
  int   blip_edepid[kMaxBlips];       // true energy dep ID
  int   blip_clustid[kNplanes][kMaxBlips];  // cluster ID per plane
  //float blip_dTmatch[kNplanes][kMaxBlips];
  float blip_trkdist[kMaxBlips];      // distance to nearest track
  
  // --- lifetime calibration info ---
  float calib_lifetime;
  float calib_lifetime_err;
  
  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    lifetime              = -999;
    timestamp             = -999;
    total_depEnergy       = -999;
    total_numElectrons    = -999;
    //gamma_depEnergy       = -999;
    //gamma_numElectrons    = -999;
    nparticles            = 0;    // --- G4 particles ---
    FillWith(isPrimary,   false);
    FillWith(trackID,     -999);
    FillWith(pdg,         -99999);
    FillWith(nDaughters,  -999);
    FillWith(mother,      -999);
    FillWith(E,           -999.);
    FillWith(endE,        -999.);
    FillWith(KE,           -999.);
    FillWith(endKE,        -999.);
    FillWith(mass,        -999.);
    FillWith(P,           -999.);
    FillWith(Px,          -999.);
    FillWith(Py,          -999.);
    FillWith(Pz,          -999.);
    FillWith(startPointx, -99999.);
    FillWith(startPointy, -99999.);
    FillWith(startPointz, -99999.);
    FillWith(endPointx,   -99999.);
    FillWith(endPointy,   -99999.);
    FillWith(endPointz,   -99999.);
    FillWith(startT,      -99999.);
    FillWith(endT,        -99999.);
    FillWith(pathlen,     -999.);
    FillWith(numElectrons,-999);
    FillWith(depEnergy,   -999.);
    FillWith(process,     "");
    nedeps                = 0;    // --- EDeps ---
    FillWith(edep_tpc,    -9);
    FillWith(edep_energy, -999);
    //FillWith(edep_energyESTAR, -999);
    FillWith(edep_charge, -999);
    FillWith(edep_x,      -99999.);
    FillWith(edep_y,      -99999.);
    FillWith(edep_z,      -99999.);
    FillWith(edep_ds,     -999);
    FillWith(edep_g4id,   -9);
    FillWith(edep_pdg,   -999);
    FillWith(edep_clustid,-9);
    FillWith(edep_blipid, -9);
    nhits                 = 0;    // --- TPC hits ---
    FillWith(hit_tpc,     -9);
    FillWith(hit_plane,   -9);
    FillWith(hit_wire,    -999);
    FillWith(hit_channel, -999);
    FillWith(hit_peakT,   -999);
    FillWith(hit_time,    -999);
    FillWith(hit_rms,     -999);
    FillWith(hit_amp,      -999);
    FillWith(hit_area,    -999);
    FillWith(hit_sumadc,  -999);
    FillWith(hit_mult,    -999);
    FillWith(hit_charge,  -999);
    FillWith(hit_isreal,  -9);
    FillWith(hit_ismatch, -9);
    FillWith(hit_isgood, -9);
    FillWith(hit_trkid,   -9);
    FillWith(hit_g4id,    -999);
    FillWith(hit_g4frac,  -9);
    FillWith(hit_g4energy,-999);
    FillWith(hit_g4charge,-999);
    FillWith(hit_clustid, -9);
    FillWith(hit_blipid,  -9);
    FillWith(hit_gof,  -9);
    ave_trkhit_amp         = -9;
    ntrks                 = 0;    // --- Tracks --- 
    FillWith(trk_id,      -999); 
    FillWith(trk_npts,    -999); 
    FillWith(trk_length,  -999);    
    FillWith(trk_startx,  -999);    
    FillWith(trk_starty,  -999);    
    FillWith(trk_startz,  -999);    
    FillWith(trk_startd,  -999);   
    FillWith(trk_endx,    -999);      
    FillWith(trk_endy,    -999);      
    FillWith(trk_endz,    -999);      
    FillWith(trk_endd,    -999);      
    //FillWith(trk_origin,  -999);
    //FillWith(trk_g4id,    -999);      
    //FillWith(trk_g4pdg,   -999);     
    //FillWith(trk_purity,  -999);    
    //FillWith(trk_pitch,   -999);     
    //FillWith(trk_ke,      -999);        
    //FillWith(trk_cosmictag,-999); 
    //FillWith(trk_cosmicscore,-999);
    //FillWith(trk_pidpdg,  -999);     
    //FillWith(trk_pidchi,  -999);   
    //FillWith(trk_bestplane,-999);  
    //nshwrs                = 0;    // --- Showers --- 
    //FillWith(shwr_id,     -999);
    //FillWith(shwr_dirx,   -999);
    //FillWith(shwr_diry,   -999);
    //FillWith(shwr_dirz,   -999);
    //FillWith(shwr_startx, -999);
    //FillWith(shwr_starty, -999);
    //FillWith(shwr_startz, -999);
    //FillWith(shwr_length, -999);
    //FillWith(shwr_openangle, -999);
    nclusts               = 0;    // --- Hit Clusters ---
    FillWith(clust_id,  -9);
    FillWith(clust_tpc,  -9);
    FillWith(clust_plane, -9);
    FillWith(clust_wire, -9);
    FillWith(clust_startwire, -9);
    FillWith(clust_endwire, -9);
    FillWith(clust_chan, -9);
    FillWith(clust_nwires, -9);
    FillWith(clust_nhits, -9);
    FillWith(clust_lhit_id, -9);
    FillWith(clust_lhit_amp, -9);
    FillWith(clust_lhit_rms, -9);
    FillWith(clust_lhit_time, -9);
    FillWith(clust_lhit_peakT, -9);
    //FillWith(clust_lhit_gof, -9);
    FillWith(clust_lhit_isfit, false);
    FillWith(clust_charge, -999);
    FillWith(clust_time, -999);
    FillWith(clust_time_w, -999);
    FillWith(clust_time_err, -999);
    FillWith(clust_startTime, -999);
    FillWith(clust_endTime, -999);
    FillWith(clust_timespan, -9);
    FillWith(clust_g4id, -9);
    FillWith(clust_g4charge, -999);
    FillWith(clust_g4energy, -999);
    FillWith(clust_ismatch, false);
    FillWith(clust_edepid, -9);
    FillWith(clust_blipid, -9);
    total_blip_energy     = -9;  // --- Blips ---
    nblips                = 0;
    FillWith(blip_id,   -9);
    FillWith(blip_tpc,  -9);
    FillWith(blip_nplanes,  -9);
    FillWith(blip_x, -99999);
    FillWith(blip_y, -99999);
    FillWith(blip_z, -99999);
    FillWith(blip_maxdiff, -9);
    FillWith(blip_charge, -999);
    FillWith(blip_energy, -999);
    FillWith(blip_trkdist, -9);
    //FillWith(blip_energyESTAR, -999);
    FillWith(blip_edepid, -9);
    for(int i=0; i<kNplanes; i++){
      FillWith(blip_clustid[i], -9);
      //FillWith(blip_dTmatch[i], -999);
      //FillWith(blip_charge[i], -999);
    }
   
    calib_lifetime = -9999;
    calib_lifetime_err = -9999;

  }

  // === Function for resizing vectors (if necessary) ===
  // To be called after numbers of hits/tracks/particles
  // in the event has been determined
  void Resize() {
    if(nparticles) process.assign(nparticles,"");
  }
      
  // === Function for initializing tree branches ===
  void MakeTree(){
    art::ServiceHandle<art::TFileService> tfs;
    tree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    tree->Branch("event",&event,"event/I");
    tree->Branch("run",&run,"run/I");
    tree->Branch("timestamp",&timestamp,"timestamp/i");
    tree->Branch("lifetime",&lifetime,"lifetime/F");
    
    if( saveTruthInfo ) {
    tree->Branch("total_depEnergy",&total_depEnergy,"total_depEnergy/F");
    tree->Branch("total_numElectrons",&total_numElectrons,"total_numElectrons/F");
    tree->Branch("nparticles",&nparticles,"nparticles/S");
    tree->Branch("isPrimary",isPrimary,"isPrimary[nparticles]/O");
    tree->Branch("trackID",trackID,"trackID[nparticles]/S");
    tree->Branch("pdg",pdg,"pdg[nparticles]/I");
    tree->Branch("nDaughters",nDaughters,"nDaughters[nparticles]/S");
    tree->Branch("mother",mother,"mother[nparticles]/S");
    tree->Branch("KE",KE,"KE[nparticles]/F");
    tree->Branch("endKE",endKE,"endKE[nparticles]/F");
    //tree->Branch("mass",mass,"mass[nparticles]/F");
    tree->Branch("P",P,"P[nparticles]/F");
    //tree->Branch("Px",Px,"Px[nparticles]/F");
    //tree->Branch("Py",Py,"Py[nparticles]/F");
    //tree->Branch("Pz",Pz,"Pz[nparticles]/F");
    tree->Branch("startPointx",startPointx,"startPointx[nparticles]/F");
    tree->Branch("startPointy",startPointy,"startPointy[nparticles]/F");
    tree->Branch("startPointz",startPointz,"startPointz[nparticles]/F");
    tree->Branch("endPointx",endPointx,"endPointx[nparticles]/F");
    tree->Branch("endPointy",endPointy,"endPointy[nparticles]/F");
    tree->Branch("endPointz",endPointz,"endPointz[nparticles]/F");
    tree->Branch("startT",startT,"startT[nparticles]/F");
    tree->Branch("endT",endT,"endT[nparticles]/F");
    tree->Branch("pathlen",pathlen,"pathlen[nparticles]/F");
    tree->Branch("numElectrons",numElectrons,"numElectrons[nparticles]/F");
    tree->Branch("depEnergy",depEnergy,"depEnergy[nparticles]/F");
    tree->Branch("process",&process);
    
    tree->Branch("nedeps",&nedeps,"nedeps/I");
    tree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
    tree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
    tree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
    tree->Branch("edep_clustid",edep_clustid,"edep_clustid[nedeps]/I"); 
    tree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
    tree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/F"); 
    tree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
    tree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
    tree->Branch("edep_z",edep_z,"edep_z[nedeps]/F"); 
    tree->Branch("edep_ds",edep_ds,"edep_ds[nedeps]/F"); 
    }
  
    tree->Branch("nhits",&nhits,"nhits/I");
    if( saveHitInfo ) {
    tree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I"); 
    tree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I"); 
    tree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I"); 
    tree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F"); 
    tree->Branch("hit_time",hit_time,"hit_time[nhits]/F"); 
    tree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F"); 
    tree->Branch("hit_amp",hit_amp,"hit_amp[nhits]/F"); 
    tree->Branch("hit_area",hit_area,"hit_area[nhits]/F"); 
    tree->Branch("hit_sumadc",hit_sumadc,"hit_sumadc[nhits]/F"); 
    tree->Branch("hit_mult",hit_mult,"hit_mult[nhits]/I"); 
    tree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
    tree->Branch("hit_ismatch",hit_ismatch,"hit_ismatch[nhits]/I");
    tree->Branch("hit_isgood",hit_isgood,"hit_isgood[nhits]/I");
    tree->Branch("hit_isreal",hit_isreal,"hit_isreal[nhits]/I"); 
    tree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
    if( saveTruthInfo ) {
    tree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
    tree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F"); 
    tree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F"); 
    tree->Branch("hit_g4charge",hit_g4charge,"hit_g4charge[nhits]/F"); 
    }
    tree->Branch("hit_clustid",hit_clustid,"hit_clustid[nhits]/I"); 
    tree->Branch("hit_blipid",hit_blipid,"hit_blipid[nhits]/I");
    tree->Branch("hit_gof",hit_gof,"hit_gof[nhits]/F");
    tree->Branch("ave_trkhit_amp",&ave_trkhit_amp,"ave_trkhit_amp/F");
    }
   
    if( saveTrkInfo ) {
    tree->Branch("ntrks",&ntrks,"ntrks/I");
    tree->Branch("trk_id",trk_id,"trk_id[ntrks]/S");       
    //tree->Branch("trk_npts",trk_npts,"trk_npts[ntrks]/S");
    tree->Branch("trk_length",trk_length,"trk_length[ntrks]/F");
    tree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
    tree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
    tree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
    //tree->Branch("trk_startd",trk_startd,"trk_startd[ntrks]/F");
    tree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
    tree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
    tree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
    //tree->Branch("trk_endd",trk_endd,"trk_endd[ntrks]/F");
    //tree->Branch("trk_origin",trk_origin,"trk_origin[ntrks]/I");
    //tree->Branch("trk_g4id",trk_g4id,"trk_g4id[ntrks]/I");
    //tree->Branch("trk_g4pdg",trk_g4pdg,"trk_g4pdg[ntrks]/I");
    //tree->Branch("trk_purity",trk_purity,"trk_purity[ntrks]/F");
    //tree->Branch("trk_pitch",trk_pitch,"trk_pitch[ntrks]/F");
    //tree->Branch("trk_ke",trk_ke,"trk_ke[ntrks]/F");
    //tree->Branch("trk_cosmictag",trk_cosmictag,"trk_cosmictag[ntrks]/I"); 
    //tree->Branch("trk_cosmicscore",trk_cosmicscore,"trk_cosmicscore[ntrks]/F");
    //tree->Branch("trk_pidpdg",trk_pidpdg,"trk_pidpdg[ntrks]/I");
    //tree->Branch("trk_pidchi",trk_pidchi,"trk_pidchi[ntrks]/F");
    //tree->Branch("trk_bestplane",trk_bestplane,"trk_bestplane[ntrks]/I");
    }

    //tree->Branch("nshwrs",&nshwrs,"nshwrs/I");
    //tree->Branch("shwr_id",shwr_id,"shwr_id[nshwrs]/I");
    //tree->Branch("shwr_dirx",shwr_dirx,"shwr_dirx[nshwrs]/F");
    //tree->Branch("shwr_diry",shwr_diry,"shwr_diry[nshwrs]/F");
    //tree->Branch("shwr_dirz",shwr_dirz,"shwr_dirz[nshwrs]/F");
    //tree->Branch("shwr_startx",shwr_startx,"shwr_startx[nshwrs]/F");
    //tree->Branch("shwr_starty",shwr_starty,"shwr_starty[nshwrs]/F");
    //tree->Branch("shwr_startz",shwr_startz,"shwr_startz[nshwrs]/F");
    //tree->Branch("shwr_length",shwr_length,"shwr_length[nshwrs]/F");
    //tree->Branch("shwr_openangle",shwr_openangle,"shwr_openangle[nshwrs]/F");

    if( saveClustInfo ) {
    tree->Branch("nclusts",&nclusts,"nclusts/I");
    tree->Branch("clust_id",&clust_id,"clust_id[nclusts]/I");
    //tree->Branch("clust_tpc",clust_tpc,"clust_tpc[nclusts]/S");
    tree->Branch("clust_plane",clust_plane,"clust_plane[nclusts]/S");
    tree->Branch("clust_wire",clust_wire,"clust_wire[nclusts]/S");
    tree->Branch("clust_startwire",clust_startwire,"clust_startwire[nclusts]/S");
    tree->Branch("clust_endwire",clust_endwire,"clust_endwire[nclusts]/S");
    //tree->Branch("clust_chan",clust_chan,"clust_chan[nclusts]/S");
    tree->Branch("clust_nwires",clust_nwires,"clust_nwires[nclusts]/S");
    tree->Branch("clust_nhits",clust_nhits,"clust_nhits[nclusts]/S");
    //tree->Branch("clust_lhit_id",clust_lhit_id,"clust_lhit_id[nclusts]/I");
    tree->Branch("clust_lhit_amp",clust_lhit_amp,"clust_lhit_amp[nclusts]/F");
    tree->Branch("clust_lhit_rms",clust_lhit_rms,"clust_lhit_rms[nclusts]/F");
    //tree->Branch("clust_lhit_peakT",clust_lhit_peakT,"clust_lhit_peakT[nclusts]/F");
    //tree->Branch("clust_lhit_gof",clust_lhit_gof,"clust_lhit_gof[nclusts]/F");
    tree->Branch("clust_lhit_isfit",clust_lhit_isfit,"clust_lhit_isfit[nclusts]/O");
    tree->Branch("clust_charge",clust_charge,"clust_charge[nclusts]/F");
    tree->Branch("clust_time",clust_time,"clust_time[nclusts]/F");
    //tree->Branch("clust_time_w",clust_time_w,"clust_time_w[nclusts]/F");
    //tree->Branch("clust_time_err",clust_time_err,"clust_time_err[nclusts]/F");
    //tree->Branch("clust_startTime",clust_startTime,"clust_startTime[nclusts]/F");
    //tree->Branch("clust_endTime",clust_endTime,"clust_endTime[nclusts]/F");
    tree->Branch("clust_timespan",clust_timespan,"clust_timespan[nclusts]/F");
    if( saveTruthInfo ) {
      tree->Branch("clust_g4charge",clust_g4charge,"clust_g4charge[nclusts]/F");
      tree->Branch("clust_g4energy",clust_g4energy,"clust_g4energy[nclusts]/F");
      tree->Branch("clust_edepid",clust_edepid,"clust_edepid[nclusts]/I");
    }
    tree->Branch("clust_ismatch",clust_ismatch,"clust_ismatch[nclusts]/O");
    tree->Branch("clust_blipid",clust_blipid,"clust_blipid[nclusts]/I");
    }

    
    tree->Branch("nblips",&nblips,"nblips/I");
    //tree->Branch("blip_tpc",blip_tpc,"blip_tpc[nblips]/S");
    tree->Branch("blip_id",blip_id,"blip_id[nblips]/S");
    tree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/S");
    tree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
    tree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
    tree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
    tree->Branch("blip_maxdiff",blip_maxdiff,"blip_maxdiff[nblips]/F");
    tree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
    tree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
    tree->Branch("blip_trkdist",blip_trkdist,"blip_trkdist[nblips]/F");
    //tree->Branch("blip_energyESTAR",blip_energyESTAR,"blip_energyESTAR[nblips]/F");
    if( saveTruthInfo ) tree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
    for(int i=0;i<kNplanes;i++) tree->Branch(Form("blip_clustid_pl%i",i),blip_clustid[i],Form("blip_clustid_pl%i[nblips]/I",i));
    //for(int i=0;i<kNplanes;i++) tree->Branch(Form("blip_dTmatch_pl%i",i),blip_dTmatch[i],Form("blip_dTmatch_pl%i[nblips]/F",i));
    tree->Branch("total_blip_energy",&total_blip_energy,"total_blip_energy/F");
   
  }
  
  void AddCalibBranches(){
    tree->Branch("calib_lifetime",&calib_lifetime,"calib_lifetime/F");
    tree->Branch("calib_lifetime_err",&calib_lifetime_err,"calib_lifetime_err/F");
  }
    
    
  // === Function for filling tree ===
  void FillTree(){ tree->Fill(); }

};//BlipAnaTreeDataStruct class



//###################################################
//  BlipAna class definition
//###################################################
class BlipAna : public art::EDAnalyzer 
{ 
  public:
  explicit BlipAna(fhicl::ParameterSet const& pset);
  virtual ~BlipAna();
  
  void beginJob();                      // called once, at start of job
  void endJob();                        // called once, at end of job
  void analyze(const art::Event& evt);  // called per event

  private:
  void PrintParticleInfo(size_t);
  void PrintG4Info(const simb::MCParticle&);
  void PrintClusterInfo(const BlipUtils::HitClust&);
  void PrintHitInfo(const BlipUtils::HitInfo&);

  // --- Data and calo objects ---
  BlipAnaTreeDataStruct*  fData;
  calo::CalorimetryAlg    fCaloAlg;
  //TGraph2D*               ESTAR;

  // --- FCL configs ---
  //std::string         fAnaTreeName;
  std::string         fHitProducer;
  std::string         fTrkProducer;
  //std::string         fShowerProducer;
  //std::string         fCaloProducer;
  //bool                fSaveTruthInfo;
  //bool                fSaveHitInfo;
  bool                fDebugMode;
  bool                fMakeDiagHists;
  float               fTrueBlipMergeDist;
  bool                fDoTrkCalibrations;
  bool                fDoHitFiltering;
  float               fMaxHitTrkLength;
  float               fMaxHitAmp;
  std::vector<float>  fMinHitRMS;
  std::vector<float>  fMaxHitRMS;
  std::vector<float>  fMinHitGOF;
  std::vector<float>  fMaxHitGOF;
  float               fHitClustWidthFact;
  int                 fHitClustWireRange;
  float               fHitMatchWidthFact;
  float               fHitMatchMaxTicks;
  float               fClustMatchMinScore;
  int                 fMaxWiresInCluster;
  float               fMaxClusterSpan;
  int                 fCaloPlane;         // use this plane for calorimetry
  bool                fPickyBlips;
  
  // --- Detector properties
  float fTickPeriod;
  float fElectronLifetime;

  // --- Counters and such ---
  bool  fIsRealData         = false;
  bool  fIsMC               = false;
  int   fNumEvents          = 0;
  int   fNumHits[3]         = {};
  int   fNumHitsTrue[3]     = {};
  int   fNum3DBlips         = 0;
  int   fNum3DBlipsTrue     = 0;

  // --- Histograms ---
  TH1D*   h_nhits[kNplanes];
  TH1D*   h_nhits_m[kNplanes];
  TH1D*   h_nhits_tm[kNplanes];
  TH1D*   h_hitamp[kNplanes];
  TH1D*   h_hitrms[kNplanes];
  TH1D*   h_hitratio[kNplanes];
  TH1D*   h_hitint[kNplanes];
  TH2D*   h_nelec_TrueVsReco[kNplanes];
  TH1D*   h_nelec_Resolution[kNplanes];
  TH1D*   h_chargecomp[kNplanes];
  TH1D*   h_hitpur[kNplanes];
  TH1D*   h_hit_dt[kNplanes];
  TH1D*   h_hit_dt_coh[kNplanes];
  TH1D*   h_hit_dtfrac[kNplanes];
  TH1D*   h_clust_matchScore[kNplanes];
  TH1D*   h_clust_matchScore_best[kNplanes];
  TH1D*   h_nmatches[3];
  TH1D*   h_hitgof_2D[3];
  TH1D*   h_hitgof_3D[3];
  TH1D*   h_clust_nwires;
  TH1D*   h_clust_timespan;

  TH1D*   h_nblips;
  TH1D*   h_nblips_tm;
  TH1D*   h_blip_qcomp;
  TH1D*   h_blip_pur;

  TH2D*   h_blip_zy;
  TH1D*   h_blip_charge;
  
  TH1D*   h_trkhit_amp;  // To be reset each event

  TH1D*   h_calib_nctrks;
  TH1D*   h_calib_trkdy;
  TH1D*   h_calib_trkdx;
  TH1D*   h_calib_trknpts;
  TH1D*   h_calib_lifetime;
  TH1D*   h_calib_lifetime_inverse;
  TH1D*   h_calib_lifetime_fracerr;
  //TH1D*   h_calib_lts;
  
  // Some truth metrics for debugging
  TH1D*   h_qres_electrons;
  TH1D*   h_qres_alphas;
  TH1D*   h_adc_factor;
 
  TH1D*   h_alpha_qdep;

  // Initialize histograms
  void InitializeHistograms(){
    art::ServiceHandle<art::TFileService> tfs;
  
    h_blip_zy   = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy   ->SetOption("COLZ");
    h_blip_charge = tfs->make<TH1D>("blip_charge","3D blip charge;Charge [e-]",500,0,100e3);
    
    h_qres_electrons= tfs->make<TH1D>("qres_electrons","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_qres_alphas   = tfs->make<TH1D>("qres_alphas","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_adc_factor    = tfs->make<TH1D>("adc_per_e","Collection plane;ADC per electron",200,0,0.01);
    h_alpha_qdep    = tfs->make<TH1D>("alpha_qdep","True charge deposited by alpha",500,0,10000);
    
    if( fMakeDiagHists ) {

      art::TFileDirectory diagDir = tfs->mkdir("Diagnostics");
      float hitMax  = 3000;    int hitBins  = 3000;
      float ampMax = 200;     int ampBins   = 1000;
      float rmsMax = 30;      int rmsBins   = 600;
      float areaMax = 1000;   int areaBins  = 1000;
      float ratioMax = 5.0;   int ratioBins = 250;
      float blipMax = 300;    int blipBins  = 300;

      if( fDoTrkCalibrations ) {
      h_calib_nctrks    = diagDir.make<TH1D>("calib_nctrks","Usable calibration tracks;Number of trks per evt",10,0,10);
      h_calib_trkdy     = diagDir.make<TH1D>("calib_trkdy","dY [cm]",300,0,300);
      h_calib_trkdx     = diagDir.make<TH1D>("calib_trkdx","dX [cm]",300,0,300);
      h_calib_trknpts   = diagDir.make<TH1D>("calib_trknpts","Number hits in calibration tracks",200,0,1000);
      h_calib_lifetime  = diagDir.make<TH1D>("calib_lifetime","Quick and dirty lifetime;Fitted lifetime per track [#mus]",1000,0,1000e3);
      h_calib_lifetime_inverse  = diagDir.make<TH1D>("calib_lifetime_inverse","Quick and dirty lifetime;Fitted 1/#tau per track [#mus]",1000,0,0.001);
      h_calib_lifetime_fracerr  = diagDir.make<TH1D>("calib_lifetime_fracerr","Quick and dirty lifetime;Fitted lifetime fractional error",1000,0,1.0);
      }
    
      h_clust_nwires = diagDir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);
      h_clust_timespan = diagDir.make<TH1D>("clust_timespan","Clusters (pre-cut);Time span",100,0,100);
      
      h_nblips        = diagDir.make<TH1D>("nblips","Reconstructed 3D blips per event",blipBins,0,blipMax);
      h_nblips_tm     = diagDir.make<TH1D>("nblips_tm","Truth-matched 3D blips per event",blipBins,0,blipMax);
      h_blip_qcomp   = diagDir.make<TH1D>("blip_qcomp","Fraction of true charge (at anode) reconstructed into 3D blips",202,0,1.01);
      h_blip_pur      = diagDir.make<TH1D>("blip_pur","Fraction of truth-matched blips",202,0,1.01);
      
      for(int i=0; i<kNplanes; i++) {
    
        if( i != fCaloPlane ) {
          h_hit_dt[i]     = diagDir.make<TH1D>(Form("pl%i_hit_dt",i),    Form("Plane %i hits;dT [ticks]",i),200,-20,20);
          
          h_hit_dtfrac[i] = diagDir.make<TH1D>(Form("pl%i_hit_dtfrac",i),Form("Plane %i hits;dT/RMS",i),200,-10,10);
          //h_hit_dt_coh[i]     = diagDir.make<TH1D>(Form("pl%i_hit_dt_coh",i),    Form("Plane %i hits;PeakTime difference [ticks]",i),200,-20,20);
          h_nmatches[i] = diagDir.make<TH1D>(Form("pl%i_nmatches",i),Form("number of plane%i clusters matched to coll plane",i),10,0,10);
          h_clust_matchScore[i] = diagDir.make<TH1D>(Form("pl%i_clust_matchScore",i),Form("match score for clusters on plane %i",i),110,0,1.1);
          h_clust_matchScore_best[i] = diagDir.make<TH1D>(Form("pl%i_clust_matchScore_best",i),Form("best match score for clusters on plane %i",i),110,0,1.1);
        }

        h_nhits[i]    = diagDir.make<TH1D>(Form("pl%i_nhits",i),  Form("Plane %i;number of hits",i),hitBins,0,hitMax);
        h_nhits_m[i]   = diagDir.make<TH1D>(Form("pl%i_nhits_planematched",i), Form("Plane %i;number of plane-matched hits",i),hitBins,0,hitMax);
        h_nhits_tm[i]   = diagDir.make<TH1D>(Form("pl%i_nhits_truthmatched",i), Form("Plane %i;number of truth-matched hits",i),hitBins,0,hitMax);
        h_chargecomp[i] = diagDir.make<TH1D>(Form("pl%i_charge_completeness",i),Form("charge completness, plane %i",i),101,0,1.01);
        h_hitpur[i]     = diagDir.make<TH1D>(Form("pl%i_hit_purity",i),Form("hit purity, plane %i",i),101,0,1.01);
        h_hitamp[i]   = diagDir.make<TH1D>(Form("pl%i_hit_amp",i), Form("Plane %i hits;hit amplitude [ADC]",i),ampBins,0,ampMax);
        h_hitrms[i]   = diagDir.make<TH1D>(Form("pl%i_hit_rms",i), Form("Plane %i hits;hit RMS [ADC]",i),rmsBins,0,rmsMax);
        h_hitratio[i] = diagDir.make<TH1D>(Form("pl%i_hit_ratio",i), Form("Plane %i hits;hit RMS/Amplitude ratio",i),ratioBins,0,ratioMax);
        h_hitint[i]   = diagDir.make<TH1D>(Form("pl%i_hit_integral",i), Form("Plane %i hits;hit integral [ADC]",i),areaBins,0,areaMax);
        h_hitgof_2D[i]= diagDir.make<TH1D>(Form("pl%i_hit_gof_2d",i), Form("Plane %i hits (NOT plane-matched);goodness of fit [ADC]",i),600,-1,10);
        h_hitgof_3D[i]= diagDir.make<TH1D>(Form("pl%i_hit_gof_3d",i), Form("Plane %i hits (plane-matched);goodness of fit [ADC]",i),600,-1,10);
       
        h_nelec_TrueVsReco[i] = diagDir.make<TH2D>( Form("pl%i_nelec_TrueVsReco",i),
          Form("Plane %i;True hit charge [ #times 10^{3} electrons ];reconstructed hit charge [ #times 10^{3} electrons ]",i),60,0,30, 60,0,30);
        h_nelec_TrueVsReco[i] ->SetOption("colz");
        h_nelec_Resolution[i] = diagDir.make<TH1D>( Form("pl%i_nelec_res",i),Form("Plane %i;hit charge resolution: (reco-true)/true",i),300,-3,3);

      }//endloop over planes
        
      
    
    }
 
 }

};//class BlipAna


//###################################################
//  BlipAna constructor and destructor
//###################################################
BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData  (nullptr)
  ,fCaloAlg (pset.get<fhicl::ParameterSet>("CaloAlg"))
{
  
  // read in fcl parameters
  fDebugMode            = pset.get<bool>        ("DebugMode",false);
  fMakeDiagHists        = pset.get<bool>        ("MakeDiagHists");
  fTrueBlipMergeDist    = pset.get<float>       ("TrueBlipMergeDist");
  fDoHitFiltering       = pset.get<bool>        ("DoHitFiltering");
  fDoTrkCalibrations    = pset.get<bool>        ("DoTrkCalibrations",false);
  fHitClustWidthFact    = pset.get<float>       ("HitClustWidthFact");
  fHitClustWireRange    = pset.get<int>         ("HitClustWireRange");
  fHitMatchWidthFact    = pset.get<float>       ("HitMatchWidthFact");
  fHitMatchMaxTicks     = pset.get<float>       ("HitMatchMaxTicks");
  fClustMatchMinScore   = pset.get<float>       ("ClustMatchMinScore",0);
  fCaloPlane            = pset.get<int>         ("CaloPlane");
  fHitProducer        = pset.get<std::string>   ("HitProducer");
  fTrkProducer        = pset.get<std::string>   ("TrkProducer");
  fMaxHitTrkLength    = pset.get<float>         ("MaxHitTrkLength");
  fMaxWiresInCluster  = pset.get<int>           ("MaxWiresInCluster");
  fMaxClusterSpan     = pset.get<float>         ("MaxClusterSpan");
  fPickyBlips         = pset.get<bool>          ("PickyBlips");
  fMaxHitAmp          = pset.get<float>         ("MaxHitAmp");
  fMinHitRMS          = pset.get<std::vector<float>>  ("MinHitRMS",  {-9999,-9999,-9999});
  fMaxHitRMS          = pset.get<std::vector<float>>  ("MaxHitRMS",  { 9999, 9999, 9999});
  fMinHitGOF          = pset.get<std::vector<float>>  ("MinHitGOF",  {-9999,-9999,-9999});
  fMaxHitGOF          = pset.get<std::vector<float>>  ("MaxHitGOF",  { 9999, 9999, 9999});

  // data tree object
  fData = new BlipAnaTreeDataStruct();
  fData ->saveTruthInfo   = pset.get<bool>        ("SaveTruthInfo",true);
  fData ->saveTrkInfo     = pset.get<bool>        ("SaveTrkInfo",true);
  fData ->saveHitInfo     = pset.get<bool>        ("SaveHitInfo",true);
  fData ->saveClustInfo   = pset.get<bool>        ("SaveClustInfo",true);
  fData ->treeName        = pset.get<std::string> ("AnaTreeName","anatree");
  fData ->Clear();
  fData ->MakeTree();

  if( fDoTrkCalibrations ) fData->AddCalibBranches();

  // initialize histograms
  InitializeHistograms();
  
  // determine tick period
  auto const* detClock = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const clockData = detClock->TPCClock();
  fTickPeriod = clockData.TickPeriod();

  
}
BlipAna::~BlipAna(){}



//###################################################
//  beginJob: retrieve relevant detector and clock
//  data here, and save them into class variables
//###################################################
void BlipAna::beginJob() {
  
  BlipUtils::InitializeUtils();

  art::ServiceHandle<geo::Geometry> geom;
  std::cout<<"Checking geometrty\n";
  std::cout<<"  DetHalfWidth "<<geom->DetHalfWidth()<<"\n";
  std::cout<<"  DetHalfHeight "<<geom->DetHalfHeight()<<"\n";
  std::cout<<"  DetLength "<<geom->DetLength()<<"\n";
  
  h_trkhit_amp = new TH1D("h_trkhit_amp","h_trkhit_amp",1000,0,1000);
  //h_calib_lts = new TH1D("h_calib_lts","h_calib_lts",1000,0,20000); 

  // Load the ESTAR lookup table
  /*
  std::string fname;
  //std::string path="ShowerEnergyReco/ESTAREnergyLookupCurve.root";
  std::string path="ESTAREnergyLookupCurve.root";
  std::string gname="ESTAR_energy_lookup_curve";
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(path,fname);
  if( fname.empty() ) {
    throw cet::exception("BlipAna") << "Could not find ESTAR lookup curve.\n";
  } else {
    TFile fin(fname.c_str(),"READ");
    if(!fin.IsOpen()) {
      throw cet::exception("BlipAna") << "Could not open ESTAR file.\n";
    } else {
      ESTAR = dynamic_cast<TGraph2D*>(fin.Get(gname.c_str()));
      if(!ESTAR) throw cet::exception("BlipAna") << "Could not read the ESTAR TGraph";
    }
  }
  */

}


//###################################################
//  Main event-by-event analysis
//###################################################
void BlipAna::analyze(const art::Event& evt)
{
  
  //*****************************************
  // New event
  //*****************************************
  fData            ->Clear();
  fData->event      = evt.id().event();
  fData->run        = evt.id().run();
  fIsRealData       = evt.isRealData();
  
  // Get timestamp
  unsigned long long int tsval = evt.time().value();
  const unsigned long int mask32 = 0xFFFFFFFFUL;
  fData->timestamp = ( tsval >> 32 ) & mask32;
  
  // Get det properties object
  auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
  //std::cout<<"TicksOffset plane "<<0<<" = "<<detProp->GetXTicksOffset(0,0,0)<<"\n";
  //std::cout<<"TicksOffset plane "<<1<<" = "<<detProp->GetXTicksOffset(1,0,0)<<"\n";
  //std::cout<<"TicksOffset plane "<<2<<" = "<<detProp->GetXTicksOffset(2,0,0)<<"\n";
  
  /*
  float tickOffsets[3];
  tickOffsets[0] = detProp->TimeOffsetU();
  tickOffsets[1] = detProp->TimeOffsetV();
  tickOffsets[2] = detProp->TimeOffsetZ();
  std::cout<<"TimeOffsetU = "<<tickOffsets[0]<<"\n";
  std::cout<<"TimeOffsetV = "<<tickOffsets[1]<<"\n";
  std::cout<<"TimeOffsetZ = "<<tickOffsets[2]<<"\n";
  */

  // Retrieve lifetime
  //fData->lifetime = detProp->ElectronLifetime(); 
  const lariov::UBElectronLifetimeProvider& elifetime_provider = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
  fElectronLifetime = elifetime_provider.Lifetime() * /*convert ms->mus*/ 1e3;
  fData->lifetime   = fElectronLifetime;

  // Tell us what's going on!
  std::cout<<"\n"
  <<"=========== BlipAna =========================\n"
  <<"Processing event "<<evt.id().event()<<" in run "<<evt.id().run()<<"; total: "<<fNumEvents<<"\n";
  std::cout<<"Lifetime is "<<fElectronLifetime<<"\n";
  fNumEvents++;
  

  //=======================================
  // Get data products for this event
  //========================================
  
  // -- geometry
  art::ServiceHandle<geo::Geometry> geom;
  
  // -- G4 particles
  art::Handle< std::vector<simb::MCParticle> > plistHandle;
  std::vector<art::Ptr<simb::MCParticle> > plist;
  if (evt.getByLabel("largeant",plistHandle))
    art::fill_ptr_vector(plist, plistHandle);
  
  // -- SimEnergyDeposits
  /*
  art::Handle< std::vector<sim::SimEnergyDeposit> > simdeplistHandle;
  std::vector<art::Ptr<simb::MCParticle> > simdeplist;
  if (evt.getByLabel("ionization",simdeplistHandle))
    art::fill_ptr_vector(simdeplist, simdeplistHandle);
  */
 
  //std::vector<sim::SimEnergyDeposit> simdepvec;
  //auto const& sim_edep_handle = evt.getValidHandle<std::vector<sim::SimEnergyDeposit>>("ionization");
  //for (auto const& edep : *sim_edep_handle ) simdepvec.push_back(edep);
  art::Handle<std::vector<sim::SimEnergyDeposit> > sedlistHandle;
  std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
  if (evt.getByLabel("ionization",sedlistHandle)) art::fill_ptr_vector(sedlist, sedlistHandle);


  // -- hits (from input module)
  art::Handle< std::vector<recob::Hit> > hitlistHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitProducer,hitlistHandle))
    art::fill_ptr_vector(hitlist, hitlistHandle);
  
  // -- hits from gaushit
  art::Handle< std::vector<recob::Hit> > hitlistHandleGH;
  std::vector<art::Ptr<recob::Hit> > hitlistGH;
  if (evt.getByLabel("gaushit",hitlistHandleGH))
    art::fill_ptr_vector(hitlistGH, hitlistHandleGH);
  
  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrkProducer,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  
  // -- hit<->track associations
  art::FindManyP<recob::Track> fmtrk(hitlistHandle,evt,fTrkProducer);
  art::FindManyP<recob::Track> fmtrkGH(hitlistHandleGH,evt,fTrkProducer);

  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = (int)plist.size();
  fData->ntrks      = (int)tracklist.size();
  fData->Resize();
  
  // flag this data as MC
  if( fData->nparticles && !fIsMC ) fIsMC = true;

  std::cout
  <<"Found "<<fData->nparticles<<" G4 particles, "
            <<sedlist.size()<<" sim energy deposits, "
            <<fData->nhits<<" hits from "<<fHitProducer
  <<"\n";
  

  //====================================
  // Save G4 particle information
  //====================================

  // Find total visible energy and number electrons drifted to wires
  BlipUtils::CalcTotalDep(fData->total_depEnergy,fData->total_numElectrons);
  std::cout<<"Total energy deposited: "<<fData->total_depEnergy<<" MeV \n";

  // Create empty vector to save all the "true" blips in the event
  std::vector<BlipUtils::TrueBlip> trueblips;

  // Create empty vector to save all particle information
  size_t nParticles = plist.size();
  std::vector<BlipUtils::ParticleInfo> pinfo(nParticles);

  // Loop through the MCParticles
  //sim::ParticleList::const_iterator itPart = plist.begin();
  for(size_t i = 0; i<nParticles; i++){
    //const simb::MCParticle* pPart = (itPart++)->second;
    auto pPart = plist[i];

    FillParticleInfo( *pPart, pinfo[i] );
    
    // Save to TTree object
    if(i<kMaxG4){
      fData->trackID[i]         = pPart->TrackId();
      fData->pdg[i]             = pPart->PdgCode();
      fData->nDaughters[i]      = pPart->NumberDaughters();
      fData->mother[i]          = pPart->Mother();
      fData->E[i]               = pinfo[i].E;
      fData->endE[i]            = pinfo[i].endE;
      fData->mass[i]            = pinfo[i].mass;
      fData->KE[i]              = pinfo[i].KE;
      fData->endKE[i]           = pinfo[i].endKE;
      fData->P[i]               = pinfo[i].P;
      fData->Px[i]              = pinfo[i].Px;
      fData->Py[i]              = pinfo[i].Py;
      fData->Pz[i]              = pinfo[i].Pz;
      fData->startPointx[i]     = pPart->Vx();
      fData->startPointy[i]     = pPart->Vy();
      fData->startPointz[i]     = pPart->Vz();
      fData->endPointx[i]       = pPart->EndPosition()[0];
      fData->endPointy[i]       = pPart->EndPosition()[1];
      fData->endPointz[i]       = pPart->EndPosition()[2];
      fData->startT[i]          = pPart->T();
      fData->endT[i]            = pPart->EndT();
      fData->pathlen[i]         = pinfo[i].pathLength;
      fData->process[i]         = pPart->Process();
      fData->depEnergy[i]       = pinfo[i].energyDep;
      fData->numElectrons[i]    = pinfo[i].numElectrons;
      //if( pinfo.pathLength ) PrintParticleInfo(i);
      fData->isPrimary[i]       = pinfo[i].isPrimary;
      if( fDebugMode ) PrintParticleInfo(i);
    }

    
    if( pPart->PdgCode() == 1000020040 ) {
      int qdep = 0;
      for(size_t ii=0; ii<sedlist.size(); ii++){
        if( sedlist[ii]->TrackID() == pPart->TrackId() ) {
          qdep += sedlist[ii]->NumElectrons();
        }
      }
      std::cout<<"alpha deposited "<<qdep<<" electrons\n";
      if( qdep > 0 ) h_alpha_qdep->Fill(qdep);
    }



    //if( pPart->PdgCode() == 11 ) h_energyDep_electrons->Fill(pinfo[i].energyDep);
    //if( pPart->PdgCode() == 22 ) h_energyDep_gammas   ->Fill(pinfo[i].energyDep);

  } // endloop over G4 particles

  // Calculate the true blips
  BlipUtils::MakeTrueBlips(pinfo, trueblips);

  // Merge and save true blip information
  MergeTrueBlips(trueblips, fTrueBlipMergeDist); 
  fData->nedeps = (int)trueblips.size();
  std::cout<<"Found "<<trueblips.size()<<" true energy deposits\n";
    for(size_t i=0; i<trueblips.size(); i++ ) {
      trueblips[i].ID       = i;
      fData->edep_tpc[i]    = trueblips.at(i).TPC;
      fData->edep_energy[i] = trueblips.at(i).Energy;
      fData->edep_charge[i] = trueblips.at(i).NumElectrons;
      fData->edep_ds[i]     = trueblips.at(i).Length;
      fData->edep_x[i]      = trueblips.at(i).Position.X();
      fData->edep_y[i]      = trueblips.at(i).Position.Y();
      fData->edep_z[i]      = trueblips.at(i).Position.Z();
      fData->edep_g4id[i]   = trueblips.at(i).LeadG4ID;
      fData->edep_pdg[i]    = trueblips.at(i).LeadG4PDG;
      if( fDebugMode ) {
        std::cout
        <<"   ~ "<<i<<"  "<<trueblips.at(i).Energy<<" MeV, "
        <<" ds= "<<trueblips.at(i).Length<<" cm, "
        <<" trkID= "<<trueblips.at(i).LeadG4ID<<", pdg "<<trueblips.at(i).LeadG4PDG<<"\n";
      }
  }
  
  
  //====================================
  // Save track information
  //====================================
  std::cout<<"Looping over tracks...\n";
  std::map<size_t,size_t> trkindexmap;
  std::map<size_t,std::vector<size_t>> trkhitMap;
  //art::FindManyP<anab::CosmicTag> fmct(tracklistHandle,evt,"pandoratag");
  for(size_t i=0; i<tracklist.size(); i++){
    const auto& startPt = tracklist[i]->Vertex();
    const auto& endPt   = tracklist[i]->End();
    fData->trk_id[i]    = tracklist[i]->ID();
    fData->trk_npts[i]  = tracklist[i]->NumberTrajectoryPoints();
    fData->trk_length[i]= tracklist[i]->Length();
    fData->trk_startx[i]= startPt.X();
    fData->trk_starty[i]= startPt.Y();
    fData->trk_startz[i]= startPt.Z();
    fData->trk_endx[i]  = endPt.X();
    fData->trk_endy[i]  = endPt.Y();
    fData->trk_endz[i]  = endPt.Z();
    fData->trk_startd[i]= BlipUtils::DistToBoundary(startPt);
    fData->trk_endd[i]  = BlipUtils::DistToBoundary(endPt);
    trkindexmap[tracklist.at(i)->ID()] = i;
    //if( fmct.isValid() ) {
    //  if( fmct.at(i).size() == 1 ) {
    //    fData->trk_cosmictag[i] = fmct.at(i).at(0)->CosmicScore();
    //  }
    //}
    //if (fmtrkGH.isValid()){
    //  if( fmtrkGH.at(i).size() ) trkhitMap[tracklist[i]->ID()].push_back(fmtrkGH.at(i)[0]->ID());
    //}
  }





  //====================================
  // Fast-and-dirty lifetime calculation 
  // using thru-going cosmics
  //====================================
  if( fDoTrkCalibrations ) {
    
    size_t caloPlane = 2;
    int   minHits   = 150;
    float dY_min    = 225.; // cm
    float dY_max    = 235.; // cm
    float dZ_min    = 50.; // cm
    float dX_min    = 100.; // cm
    float fitLim    = 300e4; // us
    int   nctrks    = 0;

    // First identify the usable tracks, if any
    std::vector<size_t> usableTrks;
    for(size_t itrk=0; itrk<tracklist.size(); itrk++){
      auto vec = (tracklist[itrk]->End()-tracklist[itrk]->Vertex());
      float dy = fabs(vec.Y());
      float dx = fabs(vec.X());
      float dz = fabs(vec.Z());
      h_calib_trkdy->Fill(dy);
      if( dy > dY_min && dy < dY_max ) {
        h_calib_trkdx->Fill(dx);
        if( dx > dX_min && dz > dZ_min ) usableTrks.push_back(itrk);
      }
    }

    if( usableTrks.size() && fmtrkGH.isValid() ) {

      float sum_lt = 0;
      float sum_lterr = 0;
      float sumSq_lterr = 0;
      int   num_lt = 0;


      // make map of track IDs <--> GausHit IDs
      for(size_t ihit=0; ihit<hitlistGH.size(); ihit++){
        if( hitlistGH[ihit]->WireID().Plane != caloPlane ) continue;
        if (fmtrkGH.at(ihit).size()) trkhitMap[fmtrkGH.at(ihit)[0]->ID()].push_back(ihit);
      }
      std::cout<<trkhitMap.size()<<"\n";

      // loop usable tracks
      for(auto& itrk : usableTrks ) {
        
        std::vector<float> ctrk_q;
        std::vector<float> ctrk_t;
        
        // loop through hits
        float t0    = 999e9;
        float tmax  = -999e9;
        for( auto& hitID : trkhitMap[tracklist[itrk]->ID()] ) {
          if( hitlistGH[hitID]->GoodnessOfFit() < 0 ) continue;
          float q = hitlistGH[hitID]->Integral();
          //float q = hitlistGH[hitID]->PeakAmplitude();
          float t = hitlistGH[hitID]->PeakTime()*fTickPeriod;
          ctrk_q.push_back(q);
          ctrk_t.push_back(t);
          if( t < t0 )    { t0    = t; }
          if( t > tmax )  { tmax  = t; }
        }

        // skip if there weren't enough hits
        int N = ctrk_q.size();
        if( N < minHits ) continue;
        h_calib_trknpts->Fill( N );
        nctrks++;
        std::cout<<"Track has "<<N<<" pts\n";

        // make temporary histogram, applying 
        // time offset so that T0 = 0
        int nbins = int(N/50.);
        TH1D h("h","h",nbins,0,tmax-t0);
        TH1D hn("hn","hn",nbins,0,tmax-t0);
        for(size_t i=0; i<ctrk_q.size(); i++){
          h.Fill(ctrk_t.at(i)-t0,ctrk_q.at(i));
          hn.Fill(ctrk_t.at(i)-t0);
        }
       
        for(int j = 1; j <= h.GetXaxis()->GetNbins(); j++){
          int n = hn.GetBinContent(j);
          float bc = h.GetBinContent(j);
          float bc_err = h.GetBinError(j);
          if( n ) {
            float bc2 = bc/n;
            h.SetBinContent(j, bc2);
            h.SetBinError(j, bc2*sqrt( pow(bc_err/bc,2) + 1./n ) );
          }
        }

        std::cout<<"The bins are: \n";
        for(int j = 1; j <= h.GetXaxis()->GetNbins(); j++){
          float xx = h.GetBinCenter(j);
          float yy = h.GetBinContent(j);
          std::cout<<"   t= "<<xx<<",  q= "<<yy<<"\n";
        }
        
        // fit the histogram
        TF1 fit("expfit","[0]*exp(-x/[1])",0,tmax-t0);
        fit.SetParameter(0,h.GetBinContent(1));
        fit.SetParameter(1,10e3);
        //fit.SetParLimits(1,0,fitLim);
        h.Fit("expfit","WR");
        float a = fit.GetParameter(1);
        float da = fit.GetParError(1);
        //std::cout<<"Chi2/NDF = "<<fit.GetChisquare()/fit.GetNDF()<<"\n";
        if( a > 0 ) {
          if( a < fitLim-1e-1 ) std::cout<<"FITGOOD\n";
            std::cout<<"Fit lifetime: "<<a<<" +/- "<<da<<"\n";
            std::cout<<"Frac error  : "<<da/a<<"\n";
            h_calib_lifetime->Fill(a);
            h_calib_lifetime_inverse->Fill(1./a);
            h_calib_lifetime_fracerr->Fill( da/a );
            sum_lt    += a;
            sum_lterr += da;
            sumSq_lterr += pow(da,2);
            num_lt++; 
            //std::cout<<"FITGOOD\n";
          //} else {
            //std::cout<<"FITBAD\n";
          //}
        }

      } // endloop over usable tracks in this event

    
      if( num_lt ) {
        fData->calib_lifetime     = sum_lt / num_lt;
        //fData->calib_lifetime_err = sum_lterr / num_lt;
        fData->calib_lifetime_err = sqrt(sumSq_lterr) / num_lt;
      }
    
    }
    h_calib_nctrks->Fill(nctrks);
  
  }





  //====================================
  // Save hit information
  //====================================
  std::vector<BlipUtils::HitInfo> hitinfo(hitlist.size());
  std::map<int,std::vector<int>> chanhitsMap;
  std::map<int,std::vector<int>> chanhitsMap_untracked;
  std::map<int,std::vector<int>> planehitsMap;
  std::cout<<"Looping over the hits...\n";
  for(size_t i=0; i<hitlist.size(); i++){
    int   chan  = hitlist[i]->Channel();
    int   wire  = hitlist[i]->WireID().Wire;
    int   plane = hitlist[i]->WireID().Plane;
    int   tpc   = hitlist[i]->WireID().TPC;
    float peakT = hitlist[i]->PeakTime();
    hitinfo[i].hitid      = i;
    hitinfo[i].hit        = hitlist[i];
    hitinfo[i].wire       = wire;
    hitinfo[i].tpc        = tpc;
    hitinfo[i].plane      = plane;
    hitinfo[i].driftTicks = peakT - detProp->GetXTicksOffset(plane,0,0);
    hitinfo[i].qcoll      = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),plane);
    
    //hitinfo[i].point2D.SetX( wire * 0.3 );            // 0.3cm wire spacing
    //hitinfo[i].point2D.SetY( peakT * 0.5 * 0.1041 ); // ~0.1041 cm/us drift
    
    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(hitlist[i]) ){
      BlipUtils::HitTruth( hitlist[i], hitinfo[i].g4id, hitinfo[i].g4frac, hitinfo[i].g4energy, hitinfo[i].g4charge);
      hitinfo[i].g4ids = BlipUtils::HitTruthIds(hitlist[i]);
      hitinfo[i].isreal = (hitinfo[i].g4id > 0);
    }
   
    // Find associated track
    if (fmtrk.isValid()){ 
      if (fmtrk.at(i).size())  hitinfo[i].trkid = trkindexmap[fmtrk.at(i)[0]->ID()];
    }
    
    // Add to the channel hit map
    chanhitsMap[chan].push_back(i);
    planehitsMap[plane].push_back(i);
    if( hitinfo[i].trkid <= 0 ) chanhitsMap_untracked[chan].push_back(i);

  }
  
  // Analyze hit statistics per plane
  int num_hits[kNplanes]={0};
  int num_hits_true[kNplanes]={0};
  float total_chargeInHits[kNplanes]={0};
  for(size_t i=0; i<hitlist.size(); i++){
    int   pl    = hitinfo[i].plane;
    int   real  = hitinfo[i].isreal;
    float q     = hitinfo[i].qcoll;
    float qtrue = hitinfo[i].g4charge;
    fNumHits[pl]++;
    num_hits[pl]++;
    h_hitamp[pl]    ->Fill(hitlist[i]->PeakAmplitude());
    h_hitrms[pl]    ->Fill(hitlist[i]->RMS());
    h_hitratio[pl]  ->Fill(hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
    h_hitint[pl]    ->Fill(hitlist[i]->Integral());
    if( real && q > 0 && qtrue > 0 ){
      fNumHitsTrue[pl]++;
      num_hits_true[pl]++;
      total_chargeInHits[pl] += qtrue;
      h_nelec_TrueVsReco[pl]->Fill(qtrue/1e3,q/1e3);
      h_nelec_Resolution[pl]->Fill((q-qtrue)/qtrue);
    }
  }
  for(size_t pl=0; pl<kNplanes; pl++){
    float qcomp = -9, pur = -9;
    if( num_hits_true[pl] ) {
      if( fData->total_numElectrons ) qcomp = total_chargeInHits[pl]/fData->total_numElectrons;
      if( num_hits[pl] )              pur   = num_hits_true[pl] / float(num_hits[pl]);
      h_chargecomp[pl] ->Fill( qcomp );
      h_hitpur[pl]     ->Fill( pur );
      std::cout<<"* plane "<<pl<<": "<<num_hits[pl]<<" hits ("<<num_hits_true[pl]<<" truth-matched) -- completenes: "<<qcomp<<", purity "<<pur<<"\n";
    }
  }
 
  //====================================
  // Save shower information
  //====================================
  /*
  for(size_t i=0; i<showerlist.size(); i++){
    fData->shwr_id[i]     = showerlist[i]->ID();
    fData->shwr_startx[i] = showerlist[i]->ShowerStart()[0]; 
    fData->shwr_starty[i] = showerlist[i]->ShowerStart()[1]; 
    fData->shwr_startz[i] = showerlist[i]->ShowerStart()[2]; 
    fData->shwr_dirx[i]   = showerlist[i]->Direction()[0];
    fData->shwr_diry[i]   = showerlist[i]->Direction()[1];
    fData->shwr_dirz[i]   = showerlist[i]->Direction()[2];
    fData->shwr_length[i] = showerlist[i]->Length();
    fData->shwr_openangle[i] = showerlist[i]->OpenAngle();
  }
  */


  //=================================================================
  // Blip Reconstruction
  //================================================================
  //  
  //  Will eventually move these into separate alg class.
  //
  //  Procedure
  //  [x] Look for hits that were not included in a track 
  //  [x] Filter hits based on hit width, etc
  //  [x] Merge together closely-spaced hits on same wires, save average peakT +/- spread
  //  [x] Merge together clusters on adjacent wires (if they match up in time)
  //  [x] Plane-to-plane time matching
  //  [x] Wire intersection check to get XYZ
  //  [x] Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_g4energy)
  
  // Create a series of masks that we'll update as we go along
  std::vector<bool> hitIsGood(hitlist.size(), true);
  std::vector<bool> hitIsClustered(hitlist.size(),false);
  
  // Track length cut
  for(size_t i=0; i<hitlist.size(); i++){
    int trkid = hitinfo[i].trkid;
    if( trkid >= 0 ) {
      if( tracklist[trkid]->Length() > fMaxHitTrkLength ) hitIsGood[i] = false;
    }
  }

  // Filter based on hit properties
  if( fDoHitFiltering ) {
    for(size_t i=0; i<hitlist.size(); i++){
      hitIsGood[i] = false;
      int plane = hitlist[i]->WireID().Plane;
      // goodness of fit
      if( hitlist[i]->GoodnessOfFit() < fMinHitGOF[plane] ) continue;
      if( hitlist[i]->GoodnessOfFit() > fMaxHitGOF[plane] ) continue;
      // hit width and amplitude
      if( hitlist[i]->RMS() < fMinHitRMS[plane] ) continue;
      if( hitlist[i]->RMS() > fMaxHitRMS[plane] ) continue;
      if( hitlist[i]->PeakAmplitude() > fMaxHitAmp ) continue;
      // we survived the gauntlet of cuts -- hit is good!
      hitIsGood[i] = true;
    }
  }


  // ---------------------------------------------------
  // Create collection of hit clusters on same wires
  // ---------------------------------------------------
  std::vector<BlipUtils::HitClust> hitclust;
  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
  
  //std::map<int,std::vector<int>> chan_clustsMap;
  std::cout<<"Doing first-pass cluster reco\n";
  for(auto const& chanhits : chanhitsMap){
    for(auto const& hi : chanhits.second ){
      //if( !hitinfo[hi].isgood || hitinfo[hi].isclustered ) continue;
      if( !hitIsGood[hi] || hitIsClustered[hi] ) continue; 
      // cluster this hit
      BlipUtils::HitClust hc = BlipUtils::MakeHitClust(hitinfo[hi]);
      if( !hc.isValid ) continue;
      //hitinfo[hi].isclustered=true;
      hitIsClustered[hi] = true;
      // see if we can add other hits to it; continue until 
      // no new hits can be lumped in with this clust
      int hitsAdded;
      do{
        hitsAdded = 0;  
        for(auto const& hj : chanhits.second ) {
          //if( !hitinfo[hj].isgood || hitinfo[hj].isclustered ) continue;
          if( !hitIsGood[hj] || hitIsClustered[hj] ) continue; 
          if( hi == hj ) continue;
          float range = fHitClustWidthFact*hitlist[hj]->RMS();
          float t1 = hitinfo[hj].driftTicks - range;
          float t2 = hitinfo[hj].driftTicks + range;
          if( (t1 > hc.StartTime && t1 < hc.EndTime )
            ||(t2 > hc.StartTime && t2 < hc.EndTime ) ){
            BlipUtils::GrowHitClust(hitinfo[hj],hc);
            //hitinfo[hj].isclustered = true;
            hitIsClustered[hj] = true;
            hitsAdded++;
          }
        }
      } while ( hitsAdded!=0 );

      hitclust.push_back(hc);
      tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(hitclust.size()-1);
    }
  }
  //std::cout<<"Reconstructed "<<hitclust.size()<<" hit clusters\n";

  // TODO: combine this next step into the previous for a more generalized merging.
  
  // -------------------------------------------------------------------------------
  // Look for clusters on adjacent wires (but same plane) and merge them together 
  // to account for scenarios where charge from a single blip is shared between wires
  // -------------------------------------------------------------------------------
  std::vector<BlipUtils::HitClust> hitclust_merged;
  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap_merged;
  for(auto const& tpcMap : tpc_planeclustsMap ) {
    for(auto const& planeclusts : tpcMap.second ) {
      for( auto const& ci : planeclusts.second ){
        if( hitclust[ci].isMerged ) continue;
        BlipUtils::HitClust hc = hitclust[ci];
        hitclust[ci].isMerged = true;
        int clustsAdded;
        do{
          clustsAdded = 0;
          for(auto const& cj : planeclusts.second ){
            if( hitclust[cj].isMerged ) continue;
            // check if clusters are on adjacent wires 
            // (within allowable wire skip range)
            int dw1 = (hitclust[cj].StartWire-hc.EndWire);
            int dw2 = (hc.StartWire-hitclust[cj].EndWire);
            if( (dw1 >= 0 && dw1 <= fHitClustWireRange) || 
                (dw2 >= 0 && dw2 <= fHitClustWireRange) ){
              // check for overlap
              if( BlipUtils::DoHitClustsOverlap(hc,hitclust[cj]) ) {
                hc = BlipUtils::MergeHitClusts(hc,hitclust[cj]);
                if( hitclust[cj].isMerged ) { clustsAdded++;}
              }
            }
          }
        } while ( clustsAdded!=0);
        
        float span = hc.EndTime - hc.StartTime;
        h_clust_nwires->Fill(hc.Wires.size());
        h_clust_timespan->Fill(span);
        
        if( (int)hc.Wires.size() <= fMaxWiresInCluster && span <= fMaxClusterSpan ) {
          hitclust_merged.push_back(hc);
          tpc_planeclustsMap_merged[hc.TPC][hc.Plane].push_back(hitclust_merged.size()-1);
        }
      }
    }
  }
  hitclust = hitclust_merged;
  tpc_planeclustsMap = tpc_planeclustsMap_merged;
  
  std::cout<<"Reconstructed "<<hitclust.size()<<" hit clusts\n";

  //--------------------------------------------
  // Save hit cluster info
  //--------------------------------------------
  fData->nclusts = (int)hitclust.size();
  for(size_t i=0; i<hitclust.size(); i++){
    hitclust[i].ID            = i;
    fData->clust_id[i]        = i;
    fData->clust_tpc[i]       = hitclust[i].TPC;
    fData->clust_plane[i]     = hitclust[i].Plane;
    fData->clust_wire[i]      = hitclust[i].LeadHit->WireID().Wire;
    fData->clust_startwire[i] = hitclust[i].StartWire;
    fData->clust_endwire[i] = hitclust[i].EndWire;
    fData->clust_chan[i]      = hitclust[i].LeadHit->Channel();
    fData->clust_nwires[i]    = (int)hitclust[i].Wires.size();
    fData->clust_nhits[i]     = (int)hitclust[i].HitIDs.size();
    fData->clust_lhit_id[i]   = hitclust[i].LeadHitID;
    fData->clust_lhit_amp[i]  = (float)hitclust[i].LeadHit->PeakAmplitude();
    fData->clust_lhit_rms[i]  = hitclust[i].LeadHit->RMS();
    fData->clust_lhit_peakT[i]= hitclust[i].LeadHit->PeakTime();
    //fData->clust_lhit_gof[i]  = hitclust[i].LeadHit->GoodnessOfFit();
    fData->clust_lhit_isfit[i]= (hitclust[i].LeadHit->GoodnessOfFit()>=0);
    fData->clust_lhit_time[i] = hitclust[i].LeadHitTime;
    fData->clust_charge[i]    = hitclust[i].Charge;
    fData->clust_time[i]      = hitclust[i].Time;
    fData->clust_time_w[i]    = hitclust[i].WeightedTime;
    fData->clust_time_err[i]  = hitclust[i].TimeErr;
    fData->clust_startTime[i] = hitclust[i].StartTime;
    fData->clust_endTime[i]   = hitclust[i].EndTime;
    fData->clust_timespan[i]  = (hitclust[i].EndTime-hitclust[i].StartTime);
    //fData->clust_isneartrk[i] = hitinfo[hitclust[i].LeadHitID].isneartrk;
    //fData->clust_trkdist2D[i] = hitinfo[hitclust[i].LeadHitID].trkdist2D;
    // tag associated hits and get true G4 energy/charge
    for(auto const& hitID : hitclust[i].HitIDs) hitinfo[hitID].clustid = i;
    for(size_t j=0; j< trueblips.size(); j++){
      int tbG4 = trueblips[j].LeadG4ID;
      if( tbG4 >= 0 && tbG4 == hitclust[i].G4ID ) {
        hitclust[i].EdepID = trueblips[j].ID;
        fData->edep_clustid[j] = hitclust[i].ID;
        fData->clust_edepid[i] = trueblips[j].ID;
        fData->clust_g4energy[i] = trueblips[j].Energy; 
        fData->clust_g4charge[i] = trueblips[j].NumElectrons;
        
        float q = hitclust[i].Charge;
        if( q > 0 && hitclust[i].Plane == 2 ) {
          int pdg = trueblips[j].LeadG4PDG;
          float qt = trueblips[j].NumElectrons;
          if( qt > 0 && trueblips[j].Energy > 2. ) {
            if( fabs(pdg) == 11 )         h_qres_electrons->Fill( (q-qt)/qt );
            if( fabs(pdg) == 1000020040 ) h_qres_alphas   ->Fill( (q-qt)/qt );
            h_adc_factor->Fill( hitclust[i].ADCs / qt );
          }
        }

        break;
      }
    }
    if( fDebugMode && hitclust[i].Plane == 2 ) PrintClusterInfo(hitclust[i]);
  }
  


  // =============================================================================
  // Plane matching and 3D blip formation
  // =============================================================================
  std::vector<BlipUtils::Blip> blips;

  // --------------------------------------
  // Method 1A: Require match between calo plane (typically collection) and
  //            1 or 2 induction planes. For every hitclust on the calo plane,
  //            do the following:
  //              1. Loop over hitclusts in one of the other planes (same TPC)
  //              2. Check for wire intersections
  //              3. Find closest-matched clust and add it to the histclust group
  //              4. Repeat for remaining plane(s)
  
  for(auto const& tpcMap : tpc_planeclustsMap ) { // loop on TPCs
   
    std::cout
    <<"Performing cluster matching in TPC "<<tpcMap.first<<", which has "<<tpcMap.second.size()<<" planes\n";
    auto planeMap = tpcMap.second;
    if( planeMap.find(fCaloPlane) != planeMap.end() ){
      int   planeA            = fCaloPlane;
      auto  hitclusts_planeA  = planeMap[planeA];
      std::cout<<"using plane "<<fCaloPlane<<" as reference/calo plane ("<<planeMap[planeA].size()<<" clusts)\n";
      for(auto const& i : hitclusts_planeA ) {
       
        //std::cout<<"Looking for matches to calo plane cluster: "; PrintClusterInfo(hitclust[i]);
        
        // initiate hit-cluster group
        std::vector<BlipUtils::HitClust> hcGroup;
        hcGroup.push_back(hitclust[i]);

        // store the best-matched hitclusts
        //std::vector<float> best_dT(kNplanes, 999);
        //std::vector<float> best_dTfrac(kNplanes, 999);
        //std::map<int,std::set<int>> cands;
        //std::map<int,int> bestMatchedClusts;
        
        // for each of the other planes, make a map of potential
        // cluster matches and for each, create a "match score"
        // that incorporates these two metrics:
        //   - # matched hits between them
        //   - average match dT-diff

        std::map<int, std::set<int>> cands;
        std::map<int, float> clust_score;
        
        //std::map<int, std::set<int>> cands_coherent;

        // loop over other planes
        for(auto  hitclusts_planeB : planeMap ) {
          int planeB = hitclusts_planeB.first;
          if( planeB == planeA ) continue;
          
          //std::cout
          //<<"checking plane "<<planeB<<": "<<hitclusts_planeB.second.size()<<" clusts\n";

          // Loop over all non-matched clusts on this plane
          for(auto const& j : hitclusts_planeB.second ) {
            if( hitclust[j].isMatched ) continue;
            if( hitclust[j].Charge > 2e6 ) continue; // guard against absurdly large clusters
             

            // Check if the two channels actually intersect
            double y,z;
            const int chA = hitclust[i].LeadHit->Channel();
            const int chB = hitclust[j].LeadHit->Channel();
            bool doTheyCross = art::ServiceHandle<geo::Geometry>()
                                ->ChannelsIntersect(chA,chB,y,z);
            
            if( !doTheyCross ) continue;
            
            // best match-score per hit for clust A
            //std::map<int, float> clustA_hit_score;
            float score_clustA = 0;

            //Now check the time separation for each hit; if any are within 
            //the maximum dT threshold, the clust, is a match candidate.
            std::set<int> hitsA = hitclust[i].HitIDs;
            std::set<int> hitsB = hitclust[j].HitIDs; 
            for(auto hitA : hitsA){
              
              float score_hitA = 0;
              //clust_i_hit_score[ hitA ] = 0;

              for(auto hitB : hitsB){
                float tA = hitinfo[hitA].driftTicks;
                float tB = hitinfo[hitB].driftTicks;
                float dT = tB-tA;
                float maxWidth = std::max(hitlist[hitA]->RMS(),hitlist[hitB]->RMS());
                float matchTol = std::min(fHitMatchMaxTicks, fHitMatchWidthFact*maxWidth);
                
                //float ctA = hitlist[hitA]->PeakTime() +  tickOffsets[planeA];
                //float ctB = hitlist[hitB]->PeakTime() +  tickOffsets[planeB];
                //float cdT = ctB-ctA;
                //std::cout<<"dT "<<dT<<"    cdT "<<cdT<<"\n";
                
                // ...........................................
                // fill dT histograms only if there is only 
                // 1 hit on each wire to be matched
                if( chanhitsMap_untracked[chA].size() == 1 && chanhitsMap_untracked[chB].size() == 1 ) {
                  h_hit_dt[planeB]->Fill(dT);
                  //h_hit_dt_coh[planeB]->Fill(cdT);
                  h_hit_dtfrac[planeB]->Fill(dT/matchTol);
                }

                //............................................
                if( fabs(dT) < matchTol ) {
                  float score = std::max(1. - fabs(dT)/matchTol,0.);
                  if( score > score_hitA ) score_hitA = score;
                  //if( score >= clust_i_hit_score[hitA] ) clust_i_hit_score[hitA] = score;
                  //std::cout<<"Found match between clust "<<i<<" and clust "<<j<<", dT "<<dT<<", maxWidth "<<maxWidth<<", matchTol "<<matchTol<<", score "<<score<<"\n";
                  cands[planeB].insert(j);
                  //clust_score[j] += score;
                  //std::cout<<"Total score so far for clust "<<j<<" = "<<clust_score[j]<<"\n";
                  //if( fabs(dT) < best_dT[planeB] ) {
                    //std::cout<<"    --> matching hit found! "<<dT<<"\n";
                    //update vectors
                    //best_dT[planeB] = dT;
                    //best_dTfrac[planeB] = dT / maxWidth;
                    //bestMatchedClusts[planeB] = j;
                  //}
                }//endif dT < hitMatchTol
                
                //if( fabs(cdT) < matchTol ) {
                //  cands_coherent[planeB].insert(j);
                //}

              }//endloop over hits in cluster j/B

              score_clustA += score_hitA;
            
            }//endloop over hits in cluster i/A (calo plane)
            
            // save final match score for clust i <--> clust j
            clust_score[j] = score_clustA / hitclust[i].HitIDs.size();
            
          }//endloop over B clusters
        }//endloop over other planes
        
        //std::cout
        //<<"    Found "<<cands.size()<<" planes with candidate matches\n";
        
        // were 2 coherent cands found?
        //if( cands_coherent.size() == 2 ) continue; 

        // loop over the candidates found on each plane
        // and select the one with the largest score
        if( cands.size() ) {
          for(auto c : cands ) {
            // cand.first = plane
            // cand.second = # cand clusters
            h_nmatches[c.first]->Fill(c.second.size());
            float bestScore = 0;
            int   bestID = -9;
            for(auto cid : c.second) {
              h_clust_matchScore[c.first]->Fill( clust_score[cid] );
              if( clust_score[cid] > bestScore ) {
                bestScore = clust_score[cid];
                bestID = cid;
              }
            }
            h_clust_matchScore_best[c.first]->Fill(bestScore);
            if( bestID >= 0 && bestScore >= fClustMatchMinScore ) {
              hcGroup.push_back(hitclust[bestID]);
              hitclust[bestID].isMatched = true;
              hitclust[i].isMatched = true;
              for(auto hit : hitclust[bestID].HitIDs) hitinfo[hit].ismatch = true;
              for(auto hit : hitclust[i].HitIDs) hitinfo[hit].ismatch = true;
            }
          }
          
          // make our new blip
          BlipUtils::Blip newBlip = BlipUtils::MakeBlip(hcGroup);

          // if it isn't valid, forget it and move on
          if( !newBlip.isValid ) continue;

          // if we are being picky...
          if( fPickyBlips ) {
            if( newBlip.NPlanes < 3           ) continue;
            if( newBlip.MaxIntersectDiff > 5  ) continue;
          }
          
          // if we made it this far, the blip is good!
          blips.push_back(newBlip);
          for(auto hc : hcGroup ) hitclust[hc.ID].BlipID = blips.size()-1;
        
        }


      /*
        if( cands.size() ) {
          
          // add matching candidates to the hitcluster group hcGroup
          hitclust[i].isMatched = true;
          for(auto hit : hitclust[i].HitIDs) hitinfo[hit].ismatch = true;
          for(auto c : cands ) {
            h_nmatches[c.first]->Fill(c.second.size());
            int bestID = bestMatchedClusts[c.first];
            // For now: pick only the match closest in dT
            if( bestID >= 0 ){
              hcGroup.push_back(hitclust[bestID]);
              hitclust[bestID].isMatched = true;
              for(auto hit : hitclust[bestID].HitIDs) hitinfo[hit].ismatch = true;
            }
          }
        
          // make our new blip
          BlipUtils::Blip newBlip = BlipUtils::MakeBlip(hcGroup);

          // if it isn't valid, forget it and move on
          if( !newBlip.isValid ) continue;

          // if we are being picky...
          if( fPickyBlips ) {
            if( newBlip.NPlanes < 3           ) continue;
            if( newBlip.MaxIntersectDiff > 1  ) continue;
          }
          
          // if we made it this far, the blip is good!
          blips.push_back(newBlip);
          for(auto hc : hcGroup ) hitclust[hc.ID].BlipID = blips.size()-1;

        }//endif check for existence of candidates on other planes
      */



      }//endloop over caloplane ("Plane A") clusters
    }//endif calo plane has clusters
  }//endloop over TPCs
  

  // -----------------------------------------------------
  // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
  // like associating blip with some nearby track/shower and using its tagged T0)
  //    Method 1: Assume a dE/dx = 2 MeV/cm for electrons, use that + local E-field
  //              calculate recombination.
  //    Method 2: ESTAR lookup table method ala ArgoNeuT
  for(size_t i=0; i<blips.size(); i++){

    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    float qColl = blips[i].Charge[fCaloPlane];
    float td    = blips[i].DriftTime;
    float depEl = qColl * exp( td / fElectronLifetime ); 
    auto const blipPos = blips[i].Position;
    float Efield = detProp->Efield(0);
    if( SCE->EnableSimEfieldSCE() ) {
      geo::Point_t point = {double(blipPos.X()), double(blipPos.Y()), double(blipPos.Z())};
      auto const EfieldOffsets = SCE->GetEfieldOffsets(point);
      Efield *= std::hypot(1+EfieldOffsets.X(), EfieldOffsets.Y(), EfieldOffsets.Z());
    }
    
    // METHOD 1
    float dEdx = 2.0; // MeV/cm
    float recomb = BlipUtils::ModBoxRecomb(dEdx,Efield);
    blips[i].Energy = depEl * (1./recomb) * 23.6e-6;

    // METHOD 2
    //std::cout<<"Calculating ESTAR energy dep...  "<<depEl<<", "<<Efield<<"\n";
    //blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 

  }

  

  // Loop over 3D tracks and determine the impact parameter for each found blip
  for(size_t i=0; i<tracklist.size(); i++){
    if( tracklist[i]->Length() > 10. ) {
      auto& a = tracklist[i]->Vertex();
      auto& b = tracklist[i]->End();
      TVector3 p1(a.X(), a.Y(), a.Z() );
      TVector3 p2(b.X(), b.Y(), b.Z() );
      for(size_t j=0; j<blips.size(); j++){
        TVector3& bp = blips[j].Position;
        float d = BlipUtils::DistToLine(p1,p2,bp);
        if(   fData->blip_trkdist[j] < 0 
          ||  d < fData->blip_trkdist[j] ) fData->blip_trkdist[j] = d;
      }
    }
  }



  // Save blip info to tree
  int nblips_matched = 0;
  int nblips_total = 0;
  float true_blip_charge = 0;
  fData->total_blip_energy = 0;
  fData->nblips = blips.size();
  for(size_t i=0; i<blips.size(); i++){
    fNum3DBlips++;
    nblips_total++;
    auto const& b = blips[i];
    fData->blip_id[i]         = i;
    fData->blip_tpc[i]        = b.TPC;
    fData->blip_nplanes[i]    = b.NPlanes;
    fData->blip_charge[i]       = b.Charge[fCaloPlane];
    fData->blip_energy[i]     = b.Energy;
    //fData->blip_energyESTAR[i]= b.EnergyESTAR;
    fData->blip_maxdiff[i]    = b.MaxIntersectDiff;
    fData->blip_x[i]          = (float)b.Position.X();
    fData->blip_y[i]          = (float)b.Position.Y();
    fData->blip_z[i]          = (float)b.Position.Z();
    
    h_blip_zy->Fill(b.Position.Z(),b.Position.Y());
    h_blip_charge->Fill(b.Charge[fCaloPlane]);
    fData->total_blip_energy += b.Energy;

    // std::cout<<"Finding true energy dep for blip "<<i<<" (associated clusters: "<<blips[i].ClustIDs.size()<<")\n";
    // find associated true edep
    float max = 0;
    for(auto hitID : b.HitIDs ) {
      if( hitID >= 0 ) hitinfo[hitID].blipid = i;
    }
    for(auto clustID : b.ClustIDs ) {
      int clustPlane  = fData->clust_plane[clustID];
      int edepid      = fData->clust_edepid[clustID];
      fData->blip_clustid[clustPlane][i] = clustID;
      if( edepid >= 0 && edepid < fData->nedeps ) {
        float E = trueblips[edepid].Energy;
        if( E > max ) {
          fData->blip_edepid[i] = edepid;
          fData->edep_blipid[edepid] = i;
          max = E;
        }
      }
    }
    int b_eid = fData->blip_edepid[i];
    if( b_eid >= 0 ) {
      nblips_matched++;
      true_blip_charge+= trueblips[b_eid].NumElectrons;
    }
  }
  h_nblips->Fill(nblips_total);
  h_nblips_tm->Fill(nblips_matched);
  if( fData->total_numElectrons ) {
    h_blip_qcomp->Fill(true_blip_charge / fData->total_numElectrons);
    if( nblips_total ) h_blip_pur    ->Fill( float(nblips_matched) / float(nblips_total) );
  }

  std::cout<<"Reconstructed "<<blips.size()<<" 3D blips";
  if( nblips_matched ) std::cout<<" ("<<nblips_matched<<" were truth-matched)";
  std::cout<<"\n";
  fNum3DBlipsTrue += nblips_matched;

  if( fDebugMode ) {
  for(size_t i=0; i<blips.size(); i++){
    auto const& b = blips[i];
    std::cout
    <<"   -- "<<i<<", TPC: "<<b.TPC
    <<"; charge: "<<b.Charge[2]
    <<"; recoEnergy: "<<b.Energy<<" MeV"
    <<"; Position: "<<b.Position.X()<<", "<<b.Position.Y()<<", "<<b.Position.Z()
    <<"; MaxIntersectDiff: "<<b.MaxIntersectDiff
    <<"; EdepID: "<<fData->blip_edepid[i]
    <<"\n";
  }
  }

  // Update clust data in Tree, so we can tie reconstructed hit clusters
  // to the blips they were grouped into
  for(size_t i=0; i<hitclust.size(); i++){
    fData->clust_ismatch[i] = hitclust[i].isMatched;
    fData->clust_blipid[i] = hitclust[i].BlipID;
  }

  // For every 3D track found in the event, check which 3D blips
  // are within a certain distance of it.
  // (TODO)
 
  // ========================================
  // Evaluate hit cuts
  // ========================================
  int num_hits_pmatch[kNplanes]={0};
  for(size_t i=0; i<hitlist.size(); i++){
    if( hitinfo[i].ismatch ) num_hits_pmatch[hitinfo[i].plane]++;
    // Save TTree data for hits (doing this here so we can
    // incorporate blip and cluster associations that were found
    // through the blip reconstruction process)
    if( i < kMaxHits ){
        fData->hit_plane[i]     = hitinfo[i].plane;
        fData->hit_wire[i]      = hitinfo[i].wire;
        fData->hit_tpc[i]       = hitinfo[i].tpc;
        fData->hit_trkid[i]     = hitinfo[i].trkid;
        fData->hit_channel[i]   = hitlist[i]->Channel();
        fData->hit_peakT[i]     = hitlist[i]->PeakTime();
        fData->hit_gof[i]       = hitlist[i]->GoodnessOfFit();
        fData->hit_rms[i]       = hitlist[i]->RMS();
        fData->hit_amp[i]	      = hitlist[i]->PeakAmplitude();
        fData->hit_area[i]      = hitlist[i]->Integral();
        fData->hit_sumadc[i]    = hitlist[i]->SummedADC();
        fData->hit_mult[i]      = hitlist[i]->Multiplicity();
        fData->hit_time[i]      = hitinfo[i].driftTicks;
        fData->hit_charge[i]    = hitinfo[i].qcoll;
        fData->hit_ismatch[i]   = hitinfo[i].ismatch;
        fData->hit_isreal[i]    = hitinfo[i].isreal;
        fData->hit_g4id[i]      = hitinfo[i].g4id;
        fData->hit_g4frac[i]    = hitinfo[i].g4frac;
        fData->hit_g4energy[i]  = hitinfo[i].g4energy;
        fData->hit_g4charge[i]  = hitinfo[i].g4charge;
        fData->hit_blipid[i]    = hitinfo[i].blipid;
        fData->hit_clustid[i]   = hitinfo[i].clustid;
        fData->hit_isgood[i]    = hitIsGood[i];
        if( hitinfo[i].blipid >= 0 )  h_hitgof_3D[hitinfo[i].plane]->Fill(fData->hit_gof[i]);
        else                          h_hitgof_2D[hitinfo[i].plane]->Fill(fData->hit_gof[i]);
      }
  }
  for(size_t pl=0; pl<kNplanes; pl++){
    h_nhits[pl]->Fill(num_hits[pl]);
    h_nhits_m[pl]->Fill(num_hits_pmatch[pl]);
    h_nhits_tm[pl]->Fill(num_hits_true[pl]);
  }
  
  

  fData->ave_trkhit_amp = h_trkhit_amp->GetMean();
  
  h_trkhit_amp->Reset();
  //h_calib_lts->Reset();

  //====================================
  // Fill TTree
  //====================================
  fData->FillTree();
}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipAna::endJob(){

  /*
  // Divide hitcut histos
  for(size_t i=0; i<kNplanes; i++) {
    h_hitCuts[i][0]->Scale(1./fNumEvents);
    h_hitCuts[i][1]->Scale(1./fNumEvents);
    h_hitCuts[i][2]->Scale(1./fNumEvents);
  }
  */

  printf("\n=============================================\n");
  printf("BlipAna Summary\n\n");
  printf("  Input hit collection      : %s\n",          fHitProducer.c_str());
  printf("  Input trk collection      : %s\n",          fTrkProducer.c_str());
  printf("  Max match window          : %3.1f ticks\n", fHitMatchMaxTicks);
  printf("  Hit match RMS fact        : x%3.1f\n",      fHitMatchWidthFact);
  printf("\n");
  printf("  Total events              : %i\n\n",        fNumEvents);
  printf("  Blips per evt             : %6.2f\n",       fNum3DBlips/float(fNumEvents));
  if(fIsMC){
  printf("  MC-matched blips per evt  : %6.2f\n",       fNum3DBlipsTrue/float(fNumEvents));
  printf("  Blip charge completeness  : %.4f\n",        h_blip_qcomp->GetMean());
  }
  printf("  Mean blip charge [e-]     : %8f\n",       h_blip_charge->GetMean());
  printf("\n");
  for(size_t i=0; i<kNplanes; i++){
  printf("  Plane %lu -------------------------\n",i);
  printf("   * total hits/evt        : %.2f\n",fNumHits[i]/(float)fNumEvents);
  if(fIsMC) {
  printf("   * true-matched hits/evt : %.2f\n",fNumHitsTrue[i]/(float)fNumEvents);
  printf("   * charge completeness   : %.4f\n",h_chargecomp[i]->GetMean());
  printf("   * hit purity            : %.4f\n",h_hitpur[i]->GetMean());
  }
  } 
  
  printf("\n=============================================\n");
  std::cout<<fDoTrkCalibrations<<"\n"; 
  if( fDoTrkCalibrations ){
    std::cout<<"h_calib_lifetime mean: "<<h_calib_lifetime->GetMean()<<"\n";
    std::cout<<"h_calib_lifetime_inv mean: "<<h_calib_lifetime_inverse->GetMean()<<"\n";
  }

}




//###################################################
//  Printouts for debugging
//###################################################

void BlipAna::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%8.2f, Edep=%8.3f, T=%10.2f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->trackID[i],
   fData->pdg[i],
   fData->startPointx[i],
   fData->startPointy[i],
   fData->startPointz[i],
   fData->pathlen[i], 
   fData->E[i]-fData->mass[i],
   fData->depEnergy[i],
   fData->startT[i]/1e3,
   fData->mother[i],
   fData->process[i].c_str(),
   fData->nDaughters[i]
  ); 
}

void BlipAna::PrintG4Info(const simb::MCParticle& part){
  printf("  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, endXYZ: %7.1f %7.1f %7.1f, KE0=%8.1f, moth=%5i, %12s, ND=%i\n",
   (int)part.TrackId(),
   (int)part.PdgCode(),
   (float)part.Vx(),
   (float)part.Vy(),
   (float)part.Vz(),
   (float)(part.EndPosition()[0]),
   (float)(part.EndPosition()[1]),
   (float)(part.EndPosition()[2]),
   float(1e3*(part.E()-part.Mass())),
   (int)part.Mother(),
   part.Process().c_str(),
   (int)part.NumberDaughters()
  ); 
}

void BlipAna::PrintHitInfo(const BlipUtils::HitInfo& hi){
  printf("  hitID: %4i, TPC: %i, plane: %i, driftTicks: %7.2f, leadWire: %3i, isreal: %2i, G4ID: %4i, recoTrack: %4i\n",
    hi.hitid,
    hi.tpc,
    hi.plane,
    hi.driftTicks,
    hi.wire,
    (int)hi.isreal,
    hi.g4id,
    hi.trkid
  );
}

void BlipAna::PrintClusterInfo(const BlipUtils::HitClust& hc){
  printf("  ID: %4i, TPC: %i, plane: %i, time: %7.2f, leadWire: %3i, mult: %3i, leadG4: %4i, edepid: %i, isMatched: %i\n",
    hc.ID,
    hc.TPC,
    hc.Plane,
    hc.Time,
    hc.LeadHit->WireID().Wire,
    (int)hc.HitIDs.size(),
    hc.G4ID,
    hc.EdepID,
    hc.isMatched
  );
}



DEFINE_ART_MODULE(BlipAna)

#endif
