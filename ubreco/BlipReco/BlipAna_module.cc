//#####################################################################
//###  BlipAna analyzer module
//###
//###  Contains algorithms for reconstructing isolated, MeV-scale energy
//###  depositions in the TPC, called "blips." A TTree is made for offline
//###  analysis and plot-making. Algs will eventually be migrated into 
//###  dedicated alg/tool classes as appropriate.
//###
//###  Author: Will Foreman
//###  Date:   Sept 2021
//#####################################################################
#ifndef BLIPANA_H
#define BLIPANA_H

// Framework includes
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "cetlib/search_path.h"

// MicroBooNE-specific includes
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include "ubreco/BlipReco/Utils/NuSelectionToolBase.h"
#include "ubreco/BlipReco/Utils/NuSelectionSCECorrections.h"
#include "ubreco/BlipReco/Utils/NuSelectionTrackShowerScoreFuncs.h"

// C++ includes
#include <cstring>
#include <utility>
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

// Helper templates for initializing arrays
namespace{  
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, ITER to, TYPE value) 
    { std::fill(from, to, value); }
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }
}


// Set global constants and max array sizes
const int kMaxHits    =  30000;
const int kMaxClusts  =  10000; 
const int kMaxTrks    =  1000;
const int kMaxBlips   = 10000;
const int kMaxG4      = 100000;
const int kMaxEDeps   = 10000;
//const int kMaxTrkPts  =   2000;  

class BlipAna;
  
//###################################################
//  Data storage structure
//###################################################
class BlipAnaTreeDataStruct 
{
  public:

  // --- TTrees
  TTree* evtTree;

  // --- Configurations and switches ---
  std::string treeName      = "anatree";
  bool  saveTruthInfo       = true;
  bool  saveTrueEDeps       = true;
  bool  savePrimaries       = false;
  bool  saveTrueParticles   = false;
  bool  saveTrkInfo         = true;
  bool  saveHitInfo         = true;
  bool  saveClustInfo       = true;
  bool  saveNuInfo          = true;

  // --- Event information ---   
  int           event;                // event number
  int           run;                  // run number
  int           subrun;               // subrun number
  unsigned int  timestamp;            // unix time of event
  float         lifetime;             // electron lifetime
  int           badchans;             // #bad chans according to wirecell
  int           longtrks;             // tracks > 5 cm
  
  // --- MCTruth neutrino info ---
  int   mctruth_nu_pdg;         // Neutrino PDG (if present, otherwise = 0)
  int   mctruth_nu_ccnc;        // CC (0) or NC (1)
  int   mctruth_nu_mode;        // interaction mode from Genie
  float mctruth_nu_vtx_x;       // vertex X 
  float mctruth_nu_vtx_y;       // vertex Y
  float mctruth_nu_vtx_z;       // vertex Z
  float mctruth_nu_KE;          // kinetic energy 

  // --- Primary particles ---
  // these are grabbed from the G4 MCParticles
  // list with more details saved (XYZ, P, etc)
  int   nprimaries;
  int   primary_g4id[kMaxG4];
  int   primary_pdg[kMaxG4];
  float primary_x0[kMaxG4];
  float primary_y0[kMaxG4];
  float primary_z0[kMaxG4];
  float primary_Px[kMaxG4];
  float primary_Py[kMaxG4];
  float primary_Pz[kMaxG4];
  float primary_T0[kMaxG4];
  float primary_xAV[kMaxG4];
  float primary_yAV[kMaxG4];
  float primary_zAV[kMaxG4];

  // --- G4 information ---
  int   nparticles;               // number of G4 particles
  bool  part_isPrimary[kMaxG4];        // is primary particle
  bool  part_isContained[kMaxG4]; // particle is contained in active volume 
  int   part_g4id[kMaxG4];          // G4 track ID
  int   part_pdg[kMaxG4];              // PDG
  int   part_nDaughters[kMaxG4];       // number of daughters
  int   part_mother[kMaxG4];           // mother particle
  float part_E[kMaxG4];                // initial energy (MeV)
  float part_KE[kMaxG4];               // initial kinetic energy (MeV)
  float part_endE[kMaxG4];             // final energy (MeV)
  float part_endKE[kMaxG4];             // final energy (MeV)
  float part_mass[kMaxG4];             // mass (MeV)
  float part_P[kMaxG4];                // momentum (MeV)
  float part_Px[kMaxG4];               // momentum x (MeV)
  float part_Py[kMaxG4];               // momentum y (MeV)
  float part_Pz[kMaxG4];               // momentum z (MeV)
  float part_startPointx[kMaxG4];      // starting x (cm)
  float part_startPointy[kMaxG4];      // starting y (cm)
  float part_startPointz[kMaxG4];      // starting y (cm)
  float part_endPointx[kMaxG4];        // ending x (cm)
  float part_endPointy[kMaxG4];        // ending y (cm)
  float part_endPointz[kMaxG4];        // ending y (cm)
  float part_startT[kMaxG4];           // starting time (us)
  float part_endT[kMaxG4];             // ending time (us)
  float part_pathlen[kMaxG4];          // path length (cm)
  int   part_numTrajPts[kMaxG4];       // number traj points
  float part_depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  int   part_depElectrons[kMaxG4];     // electrons deposited
  std::vector<std::string> part_process;// process name

  // --- True energy deposit info (derived from SimChannels and SimEnergyDeposits) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 index ("part_variable[g4id]")
  bool  edep_allchansgood[kMaxEDeps]; // charge hits all good channels
  float edep_g4qfrac[kMaxEDeps];  // fraction of total charge from lead particle
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited [MeV]
  int   edep_electrons[kMaxEDeps];// total ionization electrons deposited (e-)
  //int   edep_charge[kMaxEDeps];   // total electrons reaching anode wires (e-)
  int   edep_tdrift[kMaxEDeps];   // drift time for this energy dep (us)
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)
  float edep_dx[kMaxEDeps];       // dx (cm)
  float edep_dz[kMaxEDeps];       // dz (cm)
  int   edep_proc[kMaxEDeps];     // encodes particle process
                                  //  0 = primary
                                  //  1 = compton scatter ("compt")
                                  //  2 = photoelectric effect ("phot")
                                  //  3 = e+e- pair production ("conv")
                                  //  4 = other
  
  // --- keep track of particles that made hits/clusters on collection plane
  // note: find better way to do this ...
  bool  part_madeClustCol[kMaxG4];
  bool  edep_madeClustCol[kMaxEDeps];  // did this deposition end up in a 2D cluster? (post track-mask)

  // --- Hit information ---
  int	  nhits;                    // number of hits
  int	  hit_tpc[kMaxHits];        // tpc number
  int	  hit_plane[kMaxHits];      // plane number
  int	  hit_wire[kMaxHits];       // wire number
  int	  hit_channel[kMaxHits];    // channel ID
  float	hit_peakT[kMaxHits];      // raw peak time (tick)
  float	hit_time[kMaxHits];       // corrected peak time (tick)
  float hit_rms[kMaxHits];        // shape RMS
  float	hit_amp[kMaxHits];        // amplitude
  float	hit_area[kMaxHits];       // charge (area) in ADC units
  float hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
  int	  hit_g4trkid[kMaxHits];    // G4 TrackID of leading particle
  float hit_g4frac[kMaxHits];     // fraction of hit charge from leading MCParticle
  float hit_g4energy[kMaxHits];   // true energy
  float hit_g4charge[kMaxHits];   // true number of electrons at wire
  int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
  int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
  float hit_gof[kMaxHits];        // goodness of fit (default -1)

  // --- Track information ---
  int   ntrks;                    // number tracks
  int   trk_id[kMaxTrks];         // trackID
  int   trk_npts[kMaxTrks];       // number 3D trajectory points
  bool  trk_isMC[kMaxTrks];       // 10% of track hits matched to MC
  int   trk_g4id[kMaxTrks];       // G4 track ID of this reconstructed track 
  float trk_length[kMaxTrks];     // track length [cm]
  float trk_startx[kMaxTrks];     // starting X coordinate
  float trk_starty[kMaxTrks];     // starting Y coordinate
  float trk_startz[kMaxTrks];     // starting Z coordinate
  float trk_startd[kMaxTrks];     // starting distance to boundary
  float trk_endx[kMaxTrks];       // ending X coordinate
  float trk_endy[kMaxTrks];       // ending Y coordinate
  float trk_endz[kMaxTrks];       // ending Z coordinate
  float trk_endd[kMaxTrks];       // ending distance to boundary

  // --- Hit cluster information ---
  int   nclusts;                      // total clusters made
  int   clust_id[kMaxClusts];           // cluster ID (index)
  int   clust_tpc[kMaxClusts];          // cluster TPC ID
  int   clust_plane[kMaxClusts];        // cluster plane
  int   clust_wire[kMaxClusts];         // central-most wire of cluster
  int   clust_startwire[kMaxClusts];    // starting wire
  int   clust_endwire[kMaxClusts];      // ending wire
  int   clust_nwires[kMaxClusts];       // number of wires in this cluster
  int   clust_nwiresNoisy[kMaxClusts];
  int   clust_nwiresBad[kMaxClusts];
  //int   clust_deadwiresep[kMaxClusts];  // separation from nearest dead region (0=adjacent)
  bool  clust_bydeadwire[kMaxClusts];     // is cluster adjacent to a dead wire
  //int   clust_nticks[kMaxClusts];       // timespan in ticks
  int   clust_nhits[kMaxClusts];        // number of hits
  int   clust_charge[kMaxClusts];       // cluster charge at anode [e-]
  int   clust_chargeErr[kMaxClusts];    // cluster charge uncertainty 
  float clust_amp[kMaxClusts];          // maximum hit amplitude [ADC]
  float clust_time[kMaxClusts];         // charge-weighted time
  float clust_timespan[kMaxClusts];         // timespan
  //float clust_rms[kMaxClusts];          // charge-weighted RMS
  float clust_starttime[kMaxClusts];    // cluster start tick
  float clust_endtime[kMaxClusts];      // cluster end tick
  //int   clust_nnfhits[kMaxClusts];      // number of non-fitted hits (ie, pulse trains)
  bool  clust_pulsetrain[kMaxClusts];   // does this cluster include pulse-trains?
  //float clust_gof[kMaxClusts];          // mean goodness of fit for hits
  int   clust_blipid[kMaxClusts];       // blip ID for this nlusteer (if it was made into one)
  int   clust_edepid[kMaxClusts];       // true energy dep ID
  bool  clust_ismatch[kMaxClusts];      // was this cluster plane-matched?
  bool  clust_touchtrk[kMaxClusts];     // does this cluster connect to a trk (like a delta ray)?
  int   clust_touchtrkid[kMaxClusts];   // connected track ID

  // --- 3D Blip information ---
  int   nblips;                       // number of blips in event
  int   blip_id[kMaxBlips];           // blip ID / index
  int   blip_tpc[kMaxBlips];          // blip TPC
  int   blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_sigmayz[kMaxBlips];      // difference in wire intersection points
  float blip_dx[kMaxBlips];           // length along drift direction [cm]
  float blip_dw[kMaxBlips];           // length projected onto axis perpendicular to wire orientation
  float blip_size[kMaxBlips];         // rough size estimation based on time-tick extent and wire span
  int   blip_charge[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip reco energy [MeV]
  float blip_energyTrue[kMaxBlips];   // blip truth energy [MeV]
  float blip_yzcorr[kMaxBlips];       // YZ uniformity correction factor (already applied)
  int   blip_edepid[kMaxBlips];       // truth-matched energy dep index ("edep_variable[id]")
  int   blip_g4id[kMaxBlips];         // truth-matched MC particle G4 ID
  float blip_proxtrkdist[kMaxBlips];  // distance to nearest track
  int   blip_proxtrkid[kMaxBlips];    // index of nearest trk
  bool  blip_touchtrk[kMaxBlips];     // is blip touching track?
  int   blip_touchtrkid[kMaxBlips];   // track ID of touched track
  bool  blip_incylinder[kMaxBlips];   // is blip within a cylinder near a track
  int   blip_clustid[kNplanes][kMaxBlips];     // cluster ID per plane

  // --- Reconstructed neutrino slice information (Pandora) --
  bool      nu_isNeutrino;            // neutrino slice identified by Pandora
  vfloat_t  nu_nuscore;               // neutrino score
  int       nu_pfp_pdg;               // PDG particle best matching with reco slice
  float     nu_reco_vtx_x;            // reconstructed vertex X [cm]
  float     nu_reco_vtx_y;            // reconstructed vertex Y [cm]
  float     nu_reco_vtx_z;            // reconstructed vertex Z [cm]
  vint_t    nu_trk_id;        // trackIDs for tracks in this PFP
  vfloat_t  nu_trk_score;     // track scores for tracks in this PFP
  vfloat_t  nu_shwr_score;    // shower scores in this PFP
 
  TTree*  calibTree;
  int     acptrk_npts;
  float   acptrk_theta_xz;
  float   acptrk_theta_yz;
  float   acptrk_qratio_median;
  float   acptrk_qratio_mean;
  //float   acptrk_dEdx[kMaxTrkPts];
  //float   acptrk_tdrift[kMaxTrkPts];

  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    subrun                = -999; 
    lifetime              = -999;
    badchans              = -99;
    longtrks              = -99;
    timestamp             = -999;

    mctruth_nu_pdg    = 0;
    mctruth_nu_ccnc   = -9;
    mctruth_nu_mode   = -9;
    mctruth_nu_vtx_x  = -999;
    mctruth_nu_vtx_y  = -999;
    mctruth_nu_vtx_z  = -999;
    mctruth_nu_KE     = -999;
   
    /*
    nprimaries            = 0;    // --- G4 primaries ---
    FillWith(primary_g4id, -999);
    FillWith(primary_pdg, -99999);
    FillWith(primary_x0,  -999.);
    FillWith(primary_y0,  -999.);
    FillWith(primary_z0,  -999.);
    FillWith(primary_Px,  -999.);
    FillWith(primary_Py,  -999.);
    FillWith(primary_Pz,  -999.);
    FillWith(primary_T0,  -99999.);
    FillWith(primary_xAV,  -999.);
    FillWith(primary_yAV,  -999.);
    FillWith(primary_zAV,  -999.);
    */

    nparticles            = 0;    // --- G4 particles ---
    FillWith(part_isPrimary,   false);
    FillWith(part_isContained, false);
    //FillWith(part_madeHitCol,        false);
    FillWith(part_madeClustCol,      false);
    FillWith(part_g4id,     -999);
    FillWith(part_pdg,         -99999);
    FillWith(part_nDaughters,  -999);
    FillWith(part_mother,      -999);
    FillWith(part_E,           -999.);
    FillWith(part_endE,        -999.);
    FillWith(part_KE,           -999.);
    FillWith(part_endKE,        -999.);
    FillWith(part_mass,        -999.);
    FillWith(part_P,           -999.);
    FillWith(part_Px,          -999.);
    FillWith(part_Py,          -999.);
    FillWith(part_Pz,          -999.);
    FillWith(part_startPointx, -99999.);
    FillWith(part_startPointy, -99999.);
    FillWith(part_startPointz, -99999.);
    FillWith(part_endPointx,   -99999.);
    FillWith(part_endPointy,   -99999.);
    FillWith(part_endPointz,   -99999.);
    FillWith(part_startT,      -99999.);
    FillWith(part_endT,        -99999.);
    FillWith(part_pathlen,     -999.);
    FillWith(part_numTrajPts,   -9);
    FillWith(part_depElectrons,-999);
    FillWith(part_depEnergy,   -999.);
    FillWith(part_process,     "");
    nedeps                = 0;    // --- EDeps ---
    FillWith(edep_tpc,    -9);
    FillWith(edep_energy, -999);
    FillWith(edep_electrons,  -999);
    //FillWith(edep_charge, -999);
    FillWith(edep_tdrift, -999);
    FillWith(edep_x,      -99999.);
    FillWith(edep_y,      -99999.);
    FillWith(edep_z,      -99999.);
    FillWith(edep_dx,      -99999.);
    FillWith(edep_dz,      -99999.);
    FillWith(edep_allchansgood,  true);
    FillWith(edep_g4id,     -9);
    FillWith(edep_g4qfrac,  -9);
    FillWith(edep_pdg,   -999);
    FillWith(edep_proc,   -9);
    //FillWith(edep_madeHitCol,  false);
    FillWith(edep_madeClustCol,  false);
    FillWith(edep_blipid, -9);
    nhits                 = 0;    // --- TPC hits ---
    if( saveHitInfo ) {
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
      FillWith(hit_ismatch, -9);
      FillWith(hit_trkid,   -9);
      FillWith(hit_g4trkid, -999);
      FillWith(hit_g4frac,  -9);
      FillWith(hit_g4energy,-999);
      FillWith(hit_g4charge,-999);
      FillWith(hit_clustid, -9);
      FillWith(hit_blipid,  -9);
      FillWith(hit_gof,     -9);
    }
    ntrks                 = 0;    // --- Tracks --- 
    if( saveTrkInfo ) {
      FillWith(trk_id,      -999); 
      FillWith(trk_npts,    -999); 
      FillWith(trk_length,  -999);    
      FillWith(trk_isMC,    false);      
      FillWith(trk_g4id,    -999);      
      FillWith(trk_startx,  -999);    
      FillWith(trk_starty,  -999);    
      FillWith(trk_startz,  -999);    
      FillWith(trk_startd,  -999);   
      FillWith(trk_endx,    -999);      
      FillWith(trk_endy,    -999);      
      FillWith(trk_endz,    -999);      
      FillWith(trk_endd,    -999);      
    }
    nclusts                   = 0;    // --- Hit Clusters ---
    FillWith(clust_id,        -9);
    FillWith(clust_tpc,       -9);
    FillWith(clust_plane,     -9);
    FillWith(clust_nwires,    -9);
    FillWith(clust_nwiresNoisy,    -9);
    FillWith(clust_nwiresBad,    -9);
    FillWith(clust_bydeadwire, false);
    FillWith(clust_nhits,     -9);
    FillWith(clust_wire,      -9);
    FillWith(clust_startwire, -9);
    FillWith(clust_endwire,   -9);
    FillWith(clust_charge,    -999);
    FillWith(clust_chargeErr,    -999);
    FillWith(clust_time,      -999);
    FillWith(clust_timespan,  -9);
    FillWith(clust_starttime, -999);
    FillWith(clust_endtime,   -999);
    FillWith(clust_amp,       -9);
    FillWith(clust_pulsetrain,  false);
    FillWith(clust_edepid,    -9);
    FillWith(clust_blipid,    -9);
    FillWith(clust_ismatch,   false);
    FillWith(clust_touchtrk,   false);
    FillWith(clust_touchtrkid,  -9);
    nblips                    = 0;
    FillWith(blip_id,         -9);
    FillWith(blip_tpc,        -9);
    FillWith(blip_nplanes,    -9);
    FillWith(blip_x,          -9999);
    FillWith(blip_y,          -9999);
    FillWith(blip_z,          -9999);
    FillWith(blip_sigmayz,    -9);
    FillWith(blip_dx,         -9);
    FillWith(blip_dw,        -9);
    FillWith(blip_size,       -9);
    FillWith(blip_charge,     -999);
    FillWith(blip_energy,     -999);
    FillWith(blip_energyTrue, -999);
    FillWith(blip_yzcorr,     -9);
    FillWith(blip_proxtrkdist,-99);
    FillWith(blip_proxtrkid,  -9);
    FillWith(blip_touchtrk,   false);
    FillWith(blip_touchtrkid,  -9);
    FillWith(blip_incylinder, false);
    FillWith(blip_edepid,     -9);
    FillWith(blip_g4id,     -9);
    for(int i=0; i<kNplanes; i++){ 
      FillWith(blip_clustid[i],-9);
    }

    nu_isNeutrino           = false;
    nu_nuscore              .clear();
    nu_pfp_pdg              = -999;
    nu_reco_vtx_x           = -999;
    nu_reco_vtx_y           = -999;
    nu_reco_vtx_z           = -999;
    nu_trk_score            .clear();
    nu_trk_id               .clear();
    nu_shwr_score           .clear();  



}

  // === Function for resizing vectors (if necessary) ===
  // To be called after numbers of hits/tracks/particles
  // in the event has been determined
  void Resize() {
    if(nparticles) part_process.assign(nparticles,"");
  }
      
  // === Function for initializing tree branches ===
  void MakeTree(){
    
    art::ServiceHandle<art::TFileService> tfs;
   
    evtTree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    evtTree->Branch("event",&event,"event/I");
    evtTree->Branch("run",&run,"run/I");
    evtTree->Branch("subrun",&subrun,"subrun/I");
    evtTree->Branch("timestamp",&timestamp,"timestamp/i");
    evtTree->Branch("lifetime",&lifetime,"lifetime/F");
    evtTree->Branch("badchans",&badchans,"badchans/I");
    evtTree->Branch("longtrks",&longtrks,"longtrks/I");
      
    if( saveHitInfo ) {
      evtTree->Branch("nhits",&nhits,"nhits/I");
      //evtTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I"); 
      evtTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I"); 
      evtTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I"); 
      evtTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F"); 
      evtTree->Branch("hit_time",hit_time,"hit_time[nhits]/F"); 
      evtTree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F"); 
      evtTree->Branch("hit_amp",hit_amp,"hit_amp[nhits]/F"); 
      evtTree->Branch("hit_area",hit_area,"hit_area[nhits]/F"); 
      evtTree->Branch("hit_sumadc",hit_sumadc,"hit_sumadc[nhits]/F"); 
      evtTree->Branch("hit_mult",hit_mult,"hit_mult[nhits]/I"); 
      evtTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
      evtTree->Branch("hit_ismatch",hit_ismatch,"hit_ismatch[nhits]/I");
      evtTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
      if( saveTruthInfo ) {
      evtTree->Branch("hit_g4trkid",hit_g4trkid,"hit_g4trkid[nhits]/I");
      evtTree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F"); 
      evtTree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F"); 
      evtTree->Branch("hit_g4charge",hit_g4charge,"hit_g4charge[nhits]/F"); 
      }
      evtTree->Branch("hit_clustid",hit_clustid,"hit_clustid[nhits]/I"); 
      evtTree->Branch("hit_blipid",hit_blipid,"hit_blipid[nhits]/I");
      evtTree->Branch("hit_gof",hit_gof,"hit_gof[nhits]/F");
    }
 
    if( saveTrkInfo ) {
      evtTree->Branch("ntrks",&ntrks,"ntrks/I");
      evtTree->Branch("trk_id",trk_id,"trk_id[ntrks]/I");       
      evtTree->Branch("trk_length",trk_length,"trk_length[ntrks]/F");
      evtTree->Branch("trk_isMC",trk_isMC,"trk_isMC[ntrks]/O");
      evtTree->Branch("trk_g4id",trk_g4id,"trk_g4id[ntrks]/I");
      evtTree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
      evtTree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
      evtTree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
      evtTree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
      evtTree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
      evtTree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
    }

    if( saveClustInfo ) {
      evtTree->Branch("nclusts",        &nclusts,       "nclusts/I");
      //evtTree->Branch("clust_id",       clust_id,       "clust_id[nclusts]/I");
      evtTree->Branch("clust_plane",    clust_plane,    "clust_plane[nclusts]/I");
      //evtTree->Branch("clust_wire",     clust_wire,     "clust_wire[nclusts]/I");
      evtTree->Branch("clust_nhits",    clust_nhits,    "clust_nhits[nclusts]/I");
      evtTree->Branch("clust_nwires",   clust_nwires,   "clust_nwires[nclusts]/I");
      evtTree->Branch("clust_nwiresNoisy",   clust_nwiresNoisy,   "clust_nwiresNoisy[nclusts]/I");
      //evtTree->Branch("clust_nticks",   clust_nticks,   "clust_nticks[nclusts]/I");
      evtTree->Branch("clust_startwire",clust_startwire,"clust_startwire[nclusts]/I");
      evtTree->Branch("clust_endwire",  clust_endwire,  "clust_endwire[nclusts]/I");
      evtTree->Branch("clust_bydeadwire",   clust_bydeadwire,   "clust_bydeadwire[nclusts]/O");
      evtTree->Branch("clust_time",     clust_time,     "clust_time[nclusts]/F");
      evtTree->Branch("clust_timespan", clust_timespan, "clust_timespan[nclusts]/F");
      //evtTree->Branch("clust_deadwiresep",   clust_deadwiresep,   "clust_deadwiresep[nclusts]/I");
      //evtTree->Branch("clust_nnfhits",  clust_nnfhits,  "clust_nnfhits[nclusts]/I");
      //evtTree->Branch("clust_pulsetrain",  clust_pulsetrain,  "clust_pulsetrain[nclusts]/O");
      //evtTree->Branch("clust_starttime",clust_starttime,"clust_starttime[nclusts]/F");
      //evtTree->Branch("clust_endtime",  clust_endtime,"clust_endtime[nclusts]/F");
      evtTree->Branch("clust_charge",   clust_charge,   "clust_charge[nclusts]/I");
      evtTree->Branch("clust_chargeErr",   clust_chargeErr,   "clust_chargeErr[nclusts]/I");
      evtTree->Branch("clust_amp",      clust_amp,      "clust_amp[nclusts]/F");
      //evtTree->Branch("clust_gof",      clust_gof,      "clust_gof[nclusts]/F");
      //evtTree->Branch("clust_ratio",      clust_ratio,  "clust_ratio[nclusts]/F");
      evtTree->Branch("clust_ismatch",  clust_ismatch,  "clust_ismatch[nclusts]/O");
      evtTree->Branch("clust_touchtrk",  clust_touchtrk,  "clust_touchtrk[nclusts]/O");
      if( saveTrkInfo ) evtTree->Branch("clust_touchtrkid",  clust_touchtrkid,  "clust_touchtrkid[nclusts]/I");
      evtTree->Branch("clust_blipid",   clust_blipid,   "clust_blipid[nclusts]/I");
      if( saveTrueEDeps ) evtTree->Branch("clust_edepid",   clust_edepid,   "clust_edepid[nclusts]/I");
    }

    evtTree->Branch("nblips",&nblips,"nblips/I");
    evtTree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/I");
    evtTree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
    evtTree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
    evtTree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
    //evtTree->Branch("blip_sigmayz",blip_sigmayz,"blip_sigmayz[nblips]/F");
    evtTree->Branch("blip_dx",blip_dx,"blip_dx[nblips]/F");
    evtTree->Branch("blip_dw",blip_dw,"blip_dw[nblips]/F");
    evtTree->Branch("blip_size",blip_size,"blip_size[nblips]/F");
    evtTree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/I");
    evtTree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
    evtTree->Branch("blip_yzcorr",blip_yzcorr,"blip_yzcorr[nblips]/F");
    //evtTree->Branch("blip_energyTrue",blip_energyTrue,"blip_energyTrue[nblips]/F");
    evtTree->Branch("blip_incylinder",blip_incylinder,"blip_incylinder[nblips]/O");
    evtTree->Branch("blip_proxtrkdist",blip_proxtrkdist,"blip_proxtrkdist[nblips]/F");
    evtTree->Branch("blip_touchtrk",blip_touchtrk,"blip_touchtrk[nblips]/O");
    if( saveTrkInfo ) {
      evtTree->Branch("blip_proxtrkid",blip_proxtrkid,"blip_proxtrkid[nblips]/I");
      evtTree->Branch("blip_touchtrkid",blip_touchtrkid,"blip_touchtrkid[nblips]/I");
    }
    evtTree->Branch("blip_g4id",blip_g4id,"blip_g4id[nblips]/I");
    if( saveTrueEDeps ) evtTree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
    for(int i=0;i<kNplanes;i++) evtTree->Branch(Form("blip_pl%i_clustid",i),blip_clustid[i],Form("blip_pl%i_clustid[nblips]/I",i));
   
    if( saveNuInfo ) {
    auto vf = "std::vector<float>";
    auto vi = "std::vector<int>";
    evtTree->Branch("nu_isNeutrino",&nu_isNeutrino,"nu_isNeutrino/O");
    //evtTree->Branch("nu_nuscore",&nu_nuscore,"nu_nuscore/F");
    evtTree->Branch("nu_nuscore", vf, &nu_nuscore);
    evtTree->Branch("nu_pfp_pdg",&nu_pfp_pdg,"nu_pfp_pdg/I");
    evtTree->Branch("nu_reco_vtx_x",&nu_reco_vtx_x,"nu_reco_vtx_x/F");
    evtTree->Branch("nu_reco_vtx_y",&nu_reco_vtx_y,"nu_reco_vtx_y/F");
    evtTree->Branch("nu_reco_vtx_z",&nu_reco_vtx_z,"nu_reco_vtx_z/F");
    evtTree->Branch("nu_trk_id", vi, &nu_trk_id);
    evtTree->Branch("nu_trk_score", vf, &nu_trk_score);
    evtTree->Branch("nu_shwr_score", vf, &nu_shwr_score);
    evtTree->Branch("mctruth_nu_pdg",&mctruth_nu_pdg,"mctruth_nu_pdg/I");
    evtTree->Branch("mctruth_nu_ccnc",&mctruth_nu_ccnc,"mctruth_nu_ccnc/I");
    evtTree->Branch("mctruth_nu_mode",&mctruth_nu_mode,"mctruth_nu_mode/I");
    evtTree->Branch("mctruth_nu_vtx_x",&mctruth_nu_vtx_x,"mctruth_nu_vtx_x/F");
    evtTree->Branch("mctruth_nu_vtx_y",&mctruth_nu_vtx_y,"mctruth_nu_vtx_y/F");
    evtTree->Branch("mctruth_nu_vtx_z",&mctruth_nu_vtx_z,"mctruth_nu_vtx_z/F");
    evtTree->Branch("mctruth_nu_KE",&mctruth_nu_KE,"mctruth_nu_KE/F");
    }
 
    
    if( saveTruthInfo ) {
      
      if( savePrimaries ) {
      evtTree->Branch("nprimaries",&nprimaries,"nprimaries/I");
      evtTree->Branch("primary_g4id",primary_g4id,"primary_g4id[nprimaries]/I");
      evtTree->Branch("primary_pdg",primary_pdg,"primary_pdg[nprimaries]/I");
      evtTree->Branch("primary_Px",primary_Px,"primary_Px[nprimaries]/F");
      evtTree->Branch("primary_Py",primary_Py,"primary_Py[nprimaries]/F");
      evtTree->Branch("primary_Pz",primary_Pz,"primary_Pz[nprimaries]/F");
      evtTree->Branch("primary_x0",primary_x0,"primary_x0[nprimaries]/F");
      evtTree->Branch("primary_y0",primary_y0,"primary_y0[nprimaries]/F");
      evtTree->Branch("primary_z0",primary_z0,"primary_z0[nprimaries]/F");
      evtTree->Branch("primary_xAV",primary_xAV,"primary_xAV[nprimaries]/F");
      evtTree->Branch("primary_yAV",primary_yAV,"primary_yAV[nprimaries]/F");
      evtTree->Branch("primary_zAV",primary_zAV,"primary_zAV[nprimaries]/F");
      evtTree->Branch("primary_T0",primary_T0,"primary_T0[nprimaries]/F");
      }
      
      if( saveTrueParticles ) {
      evtTree->Branch("nparticles",&nparticles,"nparticles/I");
      evtTree->Branch("part_isPrimary",part_isPrimary,"part_isPrimary[nparticles]/O");
      evtTree->Branch("part_isContained",part_isContained,"part_isContained[nparticles]/O");
      //evtTree->Branch("part_madeHitCol",part_madeHitCol,"part_madeHitCol[nparticles]/O");
      //evtTree->Branch("part_madeClustCol",part_madeClustCol,"part_madeClustCol[nparticles]/O");
      evtTree->Branch("part_g4id",part_g4id,"part_g4id[nparticles]/I");
      evtTree->Branch("part_pdg",part_pdg,"part_pdg[nparticles]/I");
      evtTree->Branch("part_nDaughters",part_nDaughters,"part_nDaughters[nparticles]/I");
      evtTree->Branch("part_mother",part_mother,"part_mother[nparticles]/I");
      evtTree->Branch("part_KE",part_KE,"part_KE[nparticles]/F");
      //evtTree->Branch("part_endKE",part_endKE,"part_endKE[nparticles]/F");
      //evtTree->Branch("part_mass",part_mass,"part_mass[nparticles]/F");
      //evtTree->Branch("part_P",part_P,"part_P[nparticles]/F");
      //evtTree->Branch("part_Px",part_Px,"part_Px[nparticles]/F");
      //evtTree->Branch("part_Py",part_Py,"part_Py[nparticles]/F");
      //evtTree->Branch("part_Pz",part_Pz,"part_Pz[nparticles]/F");
      //evtTree->Branch("part_startPointx",part_startPointx,"part_startPointx[nparticles]/F");
      //evtTree->Branch("part_startPointy",part_startPointy,"part_startPointy[nparticles]/F");
      //evtTree->Branch("part_startPointz",part_startPointz,"part_startPointz[nparticles]/F");
      //evtTree->Branch("part_endPointx",part_endPointx,"part_endPointx[nparticles]/F");
      //evtTree->Branch("part_endPointy",part_endPointy,"part_endPointy[nparticles]/F");
      //evtTree->Branch("part_endPointz",part_endPointz,"part_endPointz[nparticles]/F");
      evtTree->Branch("part_startT",part_startT,"part_startT[nparticles]/F");
      //evtTree->Branch("part_endT",part_endT,"part_endT[nparticles]/F");
      evtTree->Branch("part_pathlen",part_pathlen,"part_pathlen[nparticles]/F");
      evtTree->Branch("part_depEnergy",part_depEnergy,"part_depEnergy[nparticles]/F");
      //evtTree->Branch("part_depElectrons",part_depElectrons,"part_depElectrons[nparticles]/I");
      //evtTree->Branch("part_numElectrons",part_numElectrons,"part_numElectrons[nparticles]/F");
      evtTree->Branch("part_process",&part_process);
      }
      
      if( saveTrueEDeps ) {
      evtTree->Branch("nedeps",&nedeps,"nedeps/I");
      evtTree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
      evtTree->Branch("edep_g4qfrac",edep_g4qfrac,"edep_g4qfrac[nedeps]/F"); 
      evtTree->Branch("edep_allchansgood",edep_allchansgood,"edep_allchansgood[nedeps]/O"); 
      //evtTree->Branch("edep_madeHitCol",edep_madeHitCol,"edep_madeHitCol[nedeps]/O"); 
      //evtTree->Branch("edep_madeClustCol",edep_madeClustCol,"edep_madeClustCol[nedeps]/O"); 
      evtTree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
      evtTree->Branch("edep_proc",edep_proc,"edep_proc[nedeps]/I"); 
      evtTree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
      //evtTree->Branch("edep_clustid",edep_clustid,"edep_clustid[nedeps]/I"); 
      evtTree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
      evtTree->Branch("edep_electrons",edep_electrons,"edep_electrons[nedeps]/I"); 
      //evtTree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/I"); 
      evtTree->Branch("edep_tdrift",edep_tdrift,"edep_tdrift[nedeps]/I"); 
      evtTree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
      evtTree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
      evtTree->Branch("edep_z",edep_z,"edep_z[nedeps]/F"); 
      evtTree->Branch("edep_dx",edep_dx,"edep_dx[nedeps]/F"); 
      evtTree->Branch("edep_dz",edep_dz,"edep_dz[nedeps]/F"); 
      }
    }
  }

  void MakeCalibTree(){
      art::ServiceHandle<art::TFileService> tfs;
      calibTree = tfs->make<TTree>("calibtree","ACPT calibration tree");
      calibTree->Branch("run",&run,"run/I");
      calibTree->Branch("subrun",&subrun,"subrun/I");
      calibTree->Branch("timestamp",&timestamp,"timestamp/i");
      calibTree->Branch("acptrk_theta_xz",&acptrk_theta_xz,"acptrk_theta_xz/F");
      calibTree->Branch("acptrk_theta_yz",&acptrk_theta_yz,"acptrk_theta_yz/F");
      calibTree->Branch("acptrk_qratio_median",&acptrk_qratio_median,"acptrk_qratio_median/F");
      calibTree->Branch("acptrk_qratio_mean",&acptrk_qratio_mean,"acptrk_qratio_mean/F");
      calibTree->Branch("acptrk_npts",&acptrk_npts,"acptrk_npts/I");
      //calibTree->Branch("acptrk_dEdx",&acptrk_dEdx,"acptrk_dEdx[acptrk_npts]/F");
      //calibTree->Branch("acptrk_tdrift",&acptrk_tdrift,"acptrk_tdrift[acptrk_npts]/F");
  }
  

};//BlipAnaTreeDataStruct class



//###################################################
//  BlipAna class definition
//###################################################
class BlipAna : public art::EDAnalyzer 
{ 
  public:
  explicit BlipAna(fhicl::ParameterSet const& pset);
  virtual ~BlipAna();
  
  //void beginJob();                      // called once, at start of job
  void endJob();                        // called once, at end of job
  void analyze(const art::Event& evt);  // called per event

  private:
  void    PrintParticleInfo(size_t);
  void    PrintTrueBlipInfo(const blipobj::TrueBlip&);
  void    PrintClusterInfo(const blipobj::HitClust&);
  void    PrintHitInfo(const blipobj::HitInfo&);
  void    PrintBlipInfo(const blipobj::Blip&);
  float   Truncate(float, double = 0.1);

  // --- Data and calo objects ---
  BlipAnaTreeDataStruct*  fData;
  blip::BlipRecoAlg       fBlipAlg;

  // --- FCL configs ---
  bool                fDebugMode;
  art::InputTag       fHitProducer;
  std::string         fTrkProducer;
  std::string         fGeantProducer;
  std::string         fSimDepProducer;
  int                 fCaloPlane;
  std::vector<bool>   fSavePlaneInfo;
  bool                fDoACPTrkCalib;
  art::InputTag       fACPTCaliProducer;

  // --- Counters and such ---
  bool  fIsRealData         = false;
  bool  fIsMC               = false;
  int   fNumEvents          = 0;
  int   fNumHits[3]         = {};
  int   fNumHitsUntracked[3]= {};
  int   fNumHitsMatched[3]  = {};
  int   fNumHitsTrue[3]     = {};
  int   fNumHitsMatchedTrue[3]  = {};
  int   fNum3DBlips         = 0;
  int   fNum3DBlips3Plane   = 0;
  int   fNum3DBlipsPicky    = 0;
  int   fNum3DBlipsTrue     = 0;

  // --- Neutrino selection tools
  using ProxyPfpColl_t = selection::ProxyPfpColl_t;
  using ProxyPfpElem_t = selection::ProxyPfpElem_t;
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fHITproducer; // hit associated to cluster
  art::InputTag fSHRproducer; // shower associated to PFP
  art::InputTag fVTXproducer; // vertex associated to PFP
  art::InputTag fPCAproducer; // PCAxis associated to PFP
  art::InputTag fMCTproducer;
  art::InputTag fTRKproducer;

  void    BuildPFPMap(const ProxyPfpColl_t&);
  void    AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                    const ProxyPfpColl_t &pfp_pxy_col,
                    std::vector<ProxyPfpElem_t> &slice_v);

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // selection tool
  //std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

  // --- Histograms ---
  TH1D*   h_part_process;
 
  TH1D*   h_nhits[kNplanes];
  TH1D*   h_nclusts[kNplanes];
  TH1D*   h_nclusts_pm[kNplanes];

  TH1D*   h_hitamp[kNplanes];
  TH1D*   h_hitamp_mip[kNplanes]; 
  TH1D*   h_hitamp_alpha[kNplanes];
  
  TH1D*   h_hitrms[kNplanes];
  TH1D*   h_hitrms_mip[kNplanes];            // hits in MIP-like cosmic tracks (data overlay only)
  TH1D*   h_hitrms_iso[kNplanes];            // isolated hits, unmatched, no truth match
  TH1D*   h_hitrms_isomatch[kNplanes];       // isolated, but plane-matched, no truth match
  TH1D*   h_hitrms_iso_true[kNplanes];       // isolated hits, unmatched, truth-matched
  TH1D*   h_hitrms_isomatch_true[kNplanes];  // isolated hits, matched, truth-matched
  TH1D*   h_hitrms_alpha[kNplanes];
  TH1D*   h_hitrms_electron[kNplanes];
  
  TH1D*   h_hitgof_mip[kNplanes];             // hits in trks with dY > 220 cm
  TH1D*   h_hitgof_iso[kNplanes];             // hits not in tracks
  TH1D*   h_hitgof_isomatch_true[kNplanes];             // hits not in tracks
  TH1D*   h_hitgof_isomatch[kNplanes];        // hit not in tracks, but plane-matched
  TH1D*   h_hitgof_alpha[kNplanes];
  TH1D*   h_hit_sigmaint[kNplanes];
  TH1D*   h_hit_adcdiff[kNplanes];

  TH2D*   h_hit_charge_vs_rms[kNplanes];
  TH2D*   h_hit_charge_vs_gof[kNplanes];
 

  TH1D*   h_hitqres[kNplanes];
  TH2D*   h_hitqres_scatter[kNplanes];
  TH2D*   h_hitqres_vs_q[kNplanes];
  //TH1D*   h_hitqres_alpha[kNplanes];
  
  TH1D*   h_hitmult[kNplanes];
  TH1D*   h_hitmult_mip[kNplanes];            // hits in MIP-like cosmic tracks (data overlay only)
  TH1D*   h_hitmult_iso[kNplanes];            // isolated hits, unmatched, no truth match
  TH1D*   h_hitmult_isomatch[kNplanes];       // isolated, but plane-matched, no truth match
  TH1D*   h_hitmult_iso_true[kNplanes];       // isolated hits, unmatched, truth-matched
  TH1D*   h_hitmult_isomatch_true[kNplanes];  // isolated hits, matched, truth-matched
  TH1D*   h_hitmult_alpha[kNplanes];   // hit not in tracks, but plane-matched
  
  //TH1D*   h_hitfit[kNplanes];
  //TH1D*   h_hitfit_mip[kNplanes];            // hits in MIP-like cosmic tracks (data overlay only)
  //TH1D*   h_hitfit_iso[kNplanes];            // isolated hits, unmatched, no truth match
  //TH1D*   h_hitfit_isomatch[kNplanes];       // isolated, but plane-matched, no truth match
  //TH1D*   h_hitfit_iso_true[kNplanes];       // isolated hits, unmatched, truth-matched
  //TH1D*   h_hitfit_isomatch_true[kNplanes];  // isolated hits, matched, truth-matched
  
 
  //TH1D*   h_hitratio_alpha[kNplanes];
  //TH1D*   h_hitfit_alpha[kNplanes];   // hit not in tracks, but plane-matched

  
  //TH1D*   h_hitqres_mip[kNplanes];            // hits in MIP-like cosmic tracks (data overlay only)
  //TH1D*   h_hitqres_iso_true[kNplanes];       // isolated hits, unmatched, truth-matched
  //TH1D*   h_hitqres_isomatch_true[kNplanes];  // isolated hits, matched, truth-matched
  
  TH2D*     h_hit_tdrift_vs_RMS;
  TH2D*     h_hit_tdrift_vs_RMSratio;

  TH1D*   h_chargecomp[kNplanes];
  
  TH1D*   h_trk_length;
  TH1D*   h_trk_xspan;
  TH1D*   h_trk_yspan;
  TH1D*   h_trks_100cm[2];
  TH1D*   h_nblips;
  TH1D*   h_nblips_picky;
  TH2D*   h_blip_zy;
  TH2D*   h_blip_zy_picky;
  TH1D*   h_nblips_tm;
  TH1D*   h_blip_nplanes;
  TH1D*   h_blip_qcomp;
  TH2D*   h_blip_reszy;
  TH1D*   h_blip_resx;
  TH2D*   h_blip_resE;
  TH1D*   h_blip_charge;
  TH1D*   h_blip_charge_picky;
  TH2D*   h_blip_charge_YU;
  TH2D*   h_blip_charge_YV;
  TH2D*   h_blip_charge_UV;
  TH2D*   h_blip_charge_YU_picky;
  TH2D*   h_blip_charge_YV_picky;
  TH2D*   h_blip_charge_UV_picky;
  TH1D*   h_clust_qres_anode;
  TH1D*   h_clust_qres_dep;
  TH2D*   h_clust_qres_vs_q;
  TH2D*   h_qratio_vs_time_sim;
  //TH2D*   h_efield_distortion_yz;
  //TH2D*   h_efield_distortion_xz;
  TH1D*   h_ACPtrk_theta_xz;
  TH1D*   h_ACPtrk_theta_yz;
  TH1D*   h_ACPtrk_dEdx;
  TH1D*   h_ACPtrk_dEdx_near;
  TH1D*   h_ACPtrk_dEdx_far;
  TH1D*   h_ACPtrk_qratio;
  TH2D*   h_ACPtrk_yz;


  // Initialize histograms
  void InitializeHistograms(){
    
    art::ServiceHandle<art::TFileService> tfs;
    
    float blipMax   = 500;    
    int blipBins    = 500;
    h_nblips        = tfs->make<TH1D>("nblips","Reconstructed 3D blips per event",blipBins,0,blipMax);
    h_nblips_picky  = tfs->make<TH1D>("nblips_picky","Reconstructed 3D blips per event (3-plane match, intersect #Delta < 1 cm)",blipBins,0,blipMax);
    h_blip_zy       = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy       ->SetOption("COLZ");
    h_blip_zy_picky = tfs->make<TH2D>("blip_zy_picky","3D blip location (3-plane match, intersect #Delta < 1 cm);Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy_picky ->SetOption("COLZ");
      
    art::TFileDirectory dir_diag = tfs->mkdir("Diagnostics");
    
    h_trk_length    = dir_diag.make<TH1D>("trk_length",";Track length [cm];Tracks per event per bin",500,0,500);
    h_trk_xspan     = dir_diag.make<TH1D>("trk_xspan",";Track dX [cm]",300,0,300);
    h_trk_yspan     = dir_diag.make<TH1D>("trk_yspan",";Track dY [cm]",300,0,300);
    h_trks_100cm[0]   = dir_diag.make<TH1D>("trks_100cm_data","Track length > 100 cm (data);Number of tracks per evd;Fraction of total per bin",50,0,50);
    h_trks_100cm[1]   = dir_diag.make<TH1D>("trks_100cm_mc",  "Track length > 100 cm (mc);Number of tracks per evd;Fraction of total per bin",50,0,50);
    

    if( fDoACPTrkCalib ) {
      h_ACPtrk_theta_xz   = dir_diag.make<TH1D>("trk_theta_xz","theta_xz",90,0,90);
      h_ACPtrk_theta_yz   = dir_diag.make<TH1D>("trk_theta_yz","theta_yz",90,0,90);
      h_ACPtrk_dEdx       = dir_diag.make<TH1D>("trk_acp_dEdx","Anode-cathode piercing tracks;Trajectory point dE/dx [MeV/cm]",80,0,8);
      h_ACPtrk_dEdx_near  = dir_diag.make<TH1D>("trk_acp_dEdx_near","20-40 cm from anode;dE/dx at x=20-40cm [MeV/cm]",80,0,8);
      h_ACPtrk_dEdx_far   = dir_diag.make<TH1D>("trk_acp_dEdx_far","220-240 cm from anode;dE/dx at x=220-240cm [MeV/cm]",80,0,8);
      h_ACPtrk_qratio     = dir_diag.make<TH1D>("trk_acp_qratio","Attenuation over 2m (using median dE/dx)",200,0,2.0);
      h_ACPtrk_yz         = dir_diag.make<TH2D>("trk_acp_yz",";Z [cm];Y [cm]",1037,0,1037,234,-117,117);
      h_ACPtrk_yz         ->SetOption("colz");
    }
    
    h_blip_nplanes    = dir_diag.make<TH1D>("blip_nplanes","Matched planes per blip",3,1,4);
    h_blip_charge     = dir_diag.make<TH1D>("blip_charge","3D blips;Charge [e-]",                             200,0,100e3);
    h_blip_charge_picky  = dir_diag.make<TH1D>("blip_charge_picky","3D blips (3-plane match, intersect #Delta < 1 cm);Charge [e-]",200,0,100e3);
   
    float qmax = 100;
    int   qbins = 200;
    h_blip_charge_YU = dir_diag.make<TH2D>("blip_charge_YU","3D blips (2-3 planes);Y Charge [#times 10^{3} e-];U Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_YU ->SetOption("COLZ");
    h_blip_charge_YV = dir_diag.make<TH2D>("blip_charge_YV","3D blips (2-3 planes);Y Charge [#times 10^{3} e-];V Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_YV ->SetOption("COLZ");
    h_blip_charge_UV = dir_diag.make<TH2D>("blip_charge_UV","3D blips (2-3 planes);U Charge [#times 10^{3} e-];V Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_UV ->SetOption("COLZ");
    
    h_blip_charge_YU_picky = dir_diag.make<TH2D>("blip_charge_YU_picky","3D blips (3 planes, #Delta < 1 cm);Y Charge [#times 10^{3} e-];U Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_YU_picky ->SetOption("COLZ");
    h_blip_charge_YV_picky = dir_diag.make<TH2D>("blip_charge_YV_picky","3D blips (3 planes, #Delta < 1 cm);Y Charge [#times 10^{3} e-];V Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_YV_picky ->SetOption("COLZ");
    h_blip_charge_UV_picky = dir_diag.make<TH2D>("blip_charge_UV_picky","3D blips (3 planes, #Delta < 1 cm);U Charge [#times 10^{3} e-];V Charge [#times 10^{3} e-]",qbins,0,qmax,qbins,0,qmax);
    h_blip_charge_UV_picky ->SetOption("COLZ");


    // MC histograms related to truth
    art::TFileDirectory dir_truth = dir_diag.mkdir("Truth");
    
    h_part_process    = dir_truth.make<TH1D>("part_process","MCParticle->Process()",5,0,5);
    auto xa = h_part_process->GetXaxis();
    xa->SetBinLabel(1,"primary");
    xa->SetBinLabel(2,"compt");
    xa->SetBinLabel(3,"phot");
    xa->SetBinLabel(4,"conv");
    xa->SetBinLabel(5,"other");
    
    h_nblips_tm       = dir_truth.make<TH1D>("nblips_tm","Truth-matched 3D blips per event",blipBins,0,blipMax);
    h_blip_qcomp      = dir_truth.make<TH1D>("blip_qcomp","Fraction of true charge (at anode) reconstructed into 3D blips",202,0,1.01);
    h_blip_reszy      = dir_truth.make<TH2D>("blip_res_zy","Blip position resolution;Z_{reco} - Z_{true} [cm];Y_{reco} - Y_{true} [cm]",150,-15,15,150,-15,15);
      h_blip_reszy    ->SetOption("colz");

    h_blip_resx       = dir_truth.make<TH1D>("blip_res_x","Blip position resolution;X_{reco} - X_{true} [cm]",150,-15,15);
    h_blip_resE       = dir_truth.make<TH2D>("blip_res_energy","Energy resolution of 3D blips;Energy [MeV];#deltaE/E_{true}",100,0,5,300,-1.5,1.5);
                      h_blip_resE ->SetOption("colz");
    
    h_clust_qres_vs_q       = dir_truth.make<TH2D>("qres_vs_q","Clusters on collection plane;True charge deposited [ #times 10^{3} e- ];Reco resolution",160,0,80,200,-1,1);
      h_clust_qres_vs_q     ->SetOption("colz");
    
    h_clust_qres_anode      = dir_truth.make<TH1D>("qres_anode","Reco charge vs true charge collected;( reco-true ) / true;Area-normalized entries",200,-1.,1.);
    h_clust_qres_dep        = dir_truth.make<TH1D>("qres_dep","Reco charge vs true charge deposited;( reco-true ) / true;Area-normalized entries",200,-1.,1.);
    h_qratio_vs_time_sim  = dir_truth.make<TH2D>("qratio_vs_time_sim",";Drift time [#mus]; Q_{anode} / Q_{dep}",44,100,2300, 1000,0.50,1.50);
    h_qratio_vs_time_sim  ->SetOption("colz");
    //h_adc_factor      = dir_truth.make<TH1D>("adc_per_e","Collection plane;ADC per electron;Area-normalized entries",200,0,0.01);
    //h_sim_lifetime    = dir_truth.make<TH1D>("sim_lifetime","Calculated attenuation: #tau = - t_{drift} / ln(Q'/Q_{0}));#tau_{e} [ms];Area-normalized entries",1200,0,1200);
    //h_sim_timeconst   = dir_truth.make<TH1D>("sim_timeconst","Calculated attenuation: 1/#tau_{e} = - ln(Q'/Q_{0}) / t_{drift};1/#tau_{e} [ms^{-1}];Area-normalized entries",100,0,0.2);
    
    


    //h_calib_nctrks    = dir_diag.make<TH1D>("calib_nctrks","Usable calibration tracks;Number of trks per evt",10,0,10);
    //h_calib_trkdy     = dir_diag.make<TH1D>("calib_trkdy","dY [cm]",300,0,300);
    //h_calib_trkdx     = dir_diag.make<TH1D>("calib_trkdx","dX [cm]",300,0,300);
    //h_calib_trknpts   = dir_diag.make<TH1D>("calib_trknpts","Number hits in calibration tracks",200,0,1000);
    //h_calib_lifetime  = dir_diag.make<TH1D>("calib_lifetime","Quick and dirty lifetime;Fitted lifetime per track [#mus]",1000,0,1000e3);
   
    float hitMax  = 15000;  int hitBins  = 1500;
    float ampMax  = 50;     int ampBins   = 250;
    float rmsMax  = 20;     int rmsBins   = 200;
    //float ratioMax = 5;     int ratioBins = 250;
 
    for(int i=kNplanes-1; i >= 0; i--) {
      h_nhits[i]      = dir_diag.make<TH1D>(Form("pl%i_nhits",i),  Form("Plane %i;total number of hits",i),hitBins,0,hitMax);
      //h_nhits_ut[i]   = dir_diag.make<TH1D>(Form("pl%i_nhits_untracked",i),  Form("Plane %i;total number of untracked hits",i),hitBins,0,hitMax);
      //h_nhits_m[i]    = dir_diag.make<TH1D>(Form("pl%i_nhits_planematched",i), Form("Plane %i;total number of untracked plane-matched hits",i),hitBins,0,hitMax);
      h_hitamp[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_amp",i), Form("Plane %i untracked hits;hit amplitude [ADC]",i),ampBins,0,ampMax);
      h_hitamp_mip[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_amp_mip",i),         Form("Plane %i mip hits;hit amplitude [ADC]",i),ampBins,0,ampMax);
      h_hitamp_alpha[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_amp_alpha",i),         Form("Plane %i alpha hits;hit amplitude [ADC]",i),ampBins,0,ampMax);

      h_hitrms[i]               = dir_diag.make<TH1D>(Form("pl%i_hit_rms",i),               Form("Plane %i hits (isolated);RMS [ADC time-tick]",i),                                rmsBins,0,rmsMax);
      h_hitrms_mip[i]           = dir_diag.make<TH1D>(Form("pl%i_hit_rms_mip",i),           Form("Plane %i hits in cosmic #mu tracks;RMS [ADC time-tick]",i),                       rmsBins,0,rmsMax);
      h_hitrms_iso[i]           = dir_diag.make<TH1D>(Form("pl%i_hit_rms_iso",i),           Form("Plane %i hits (isolated, no truth match);RMS [ADC time-tick]",i),rmsBins,0,rmsMax);
      h_hitrms_iso_true[i]      = dir_diag.make<TH1D>(Form("pl%i_hit_rms_iso_true",i),      Form("Plane %i hits (isolated, no plane match, truth-matched);RMS [ADC time-tick]",i), rmsBins,0,rmsMax);
      h_hitrms_isomatch[i]      = dir_diag.make<TH1D>(Form("pl%i_hit_rms_isomatch",i),      Form("Plane %i hits (isolated, plane-matched, no truth match);RMS [ADC time-tick]",i), rmsBins,0,rmsMax);
      h_hitrms_isomatch_true[i] = dir_diag.make<TH1D>(Form("pl%i_hit_rms_isomatch_true",i), Form("Plane %i hits (isolated, plane-matched, truth-matched);RMS [ADC time-tick]",i),  rmsBins,0,rmsMax);
      h_hitrms_alpha[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_rms_alpha",i),               Form("Plane %i alpha hits;hit RMS [ADC time-tick]",i),rmsBins,0,rmsMax);
      h_hitrms_electron[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_rms_electron",i),               Form("Plane %i primary electron hits;hit RMS [ADC time-tick]",i),rmsBins,0,rmsMax);
      
      //h_hitratio[i]               = dir_diag.make<TH1D>(Form("pl%i_hit_ratio",i),               Form("Plane %i hits (isolated);RMS/amp",i),                                ratioBins,0,ratioMax);
      //h_hitratio_mip[i]           = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_mip",i),           Form("Plane %i hits in cosmic #mu tracks;RMS/amp",i),                       ratioBins,0,ratioMax);
      //h_hitratio_iso[i]           = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_iso",i),           Form("Plane %i hits (isolated, no plane match, no truth match);RMS/amp",i),ratioBins,0,ratioMax);
      //h_hitratio_iso_true[i]      = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_iso_true",i),      Form("Plane %i hits (isolated, no plane match, truth-matched);RMS/amp",i), ratioBins,0,ratioMax);
      //h_hitratio_isomatch[i]      = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_isomatch",i),      Form("Plane %i hits (isolated, plane-matched, no truth match);RMS/amp",i), ratioBins,0,ratioMax);
      //h_hitratio_isomatch_true[i] = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_isomatch_true",i), Form("Plane %i hits (isolated, plane-matched, truth-matched);RMS/amp",i),  ratioBins,0,ratioMax);
      //h_hitratio_alpha[i]         = dir_diag.make<TH1D>(Form("pl%i_hit_ratio_alpha",i),         Form("Plane %i alpha hits;hit RMS/Amplitude ratio",i),                      ratioBins,0,ratioMax);

      h_hitmult[i]                = dir_diag.make<TH1D>(Form("pl%i_hit_mult",i),                Form("Plane %i hits (isolated);Fit multiplicity ",i),                                 15,0,15);
      h_hitmult_mip[i]            = dir_diag.make<TH1D>(Form("pl%i_hit_mult_mip",i),            Form("Plane %i hits in cosmic #mu tracks;Fit multiplicity ",i),                        15,0,15);
      h_hitmult_iso[i]            = dir_diag.make<TH1D>(Form("pl%i_hit_mult_iso",i),            Form("Plane %i hits (isolated, no plane match, no truth match);Fit multiplicity ",i), 15,0,15);
      h_hitmult_iso_true[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_mult_iso_true",i),       Form("Plane %i hits (isolated, no plane match, truth-matched);Fit multiplicity ",i),  15,0,15);
      h_hitmult_isomatch[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_mult_isomatch",i),       Form("Plane %i hits (isolated, plane-matched, no truth match);Fit multiplicity ",i),  15,0,15);
      h_hitmult_isomatch_true[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_mult_isomatch_true",i),  Form("Plane %i hits (isolated, plane-matched, truth-matched);Fit multiplicity ",i),   15,0,15);
      h_hitmult_alpha[i]          = dir_diag.make<TH1D>(Form("pl%i_hit_mult_alpha",i),          Form("Plane %i alpha hits;Fit multiplicity ",i),                                       15,0,15);
      

      h_hitgof_mip[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_gof_mip",i),     Form("Plane %i hits in cosmic #mu tracks (mult=1);log10(GOF/ndf)",i),200,-10,10);
      h_hitgof_iso[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_gof_iso",i),     Form("Plane %i hits (mult=1, isolated, unmatched);log_{10}(GOF/ndf)",i),200,-10,10);       
      h_hitgof_isomatch[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_gof_isomatch",i),Form("Plane %i hits (mult=1, isolated, no plane match, no truth match);log_{10}(GOF/ndf)",i),200,-10,10);
      h_hitgof_isomatch_true[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_gof_isomatch_true",i),Form("Plane %i hits (mult=1, isolated, truth-matched);log_{10}(GOF/ndf)",i),200,-10,10);
      h_hitgof_alpha[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_gof_alpha",i),         Form("Plane %i alpha hits;log_{10}(GOF/ndf)",i),200,-10,10);
      
      //h_hitfit[i]                = dir_diag.make<TH1D>(Form("pl%i_hit_isfit",i),                Form("Plane %i hits (isolated);GOF>0",i),                                 2,0,2);
      //h_hitfit_mip[i]            = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_mip",i),            Form("Plane %i hits in cosmic #mu tracks;GOF>0",i),                        2,0,2);
      //h_hitfit_iso[i]            = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_iso",i),            Form("Plane %i hits (isolated, no plane match, no truth match);GOF>0",i), 2,0,2);
      //h_hitfit_iso_true[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_iso_true",i),       Form("Plane %i hits (isolated, no plane match, truth-matched);GOF>0",i),  2,0,2);
      //h_hitfit_isomatch[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_isomatch",i),       Form("Plane %i hits (isolated, plane-matched, no truth match);GOF>0",i),  2,0,2);
      //h_hitfit_isomatch_true[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_isomatch_true",i),  Form("Plane %i hits (isolated, plane-matched, truth-matched);GOF>0",i),   2,0,2);
      //h_hitfit_alpha[i]          = dir_diag.make<TH1D>(Form("pl%i_hit_isfit_alpha",i),          Form("Plane %i alpha hits;GOF>0",i),                                       2,0,2);
      
      h_hit_charge_vs_rms[i]  = dir_diag.make<TH2D>(Form("pl%i_hit_charge_vs_rms",i), Form("Plane %i isolated hits;hit charge [#times 10^{3} e-];hit RMS [ticks]",i),200,0,50, 200,0,10);
      h_hit_charge_vs_rms[i] ->SetOption("colz");
      //h_hit_charge_vs_ratio[i]  = dir_diag.make<TH2D>(Form("pl%i_hit_charge_vs_ratio",i), Form("Plane %i isolated hits;hit charge [#times 10^{3} e-];hit RMS/amplitude [ADC]",i),200,0,50,200,0,5);
      //h_hit_charge_vs_ratio[i] ->SetOption("colz");
      h_hit_charge_vs_gof[i]  = dir_diag.make<TH2D>(Form("pl%i_hit_charge_vs_gof",i), Form("Plane %i isolated hits;hit charge [#times 10^{3} e-];log_{10}(GOF/ndf)",i),200,0,50,180,-6,3);
      h_hit_charge_vs_gof[i] ->SetOption("colz");
      h_hit_sigmaint[i]       =dir_diag.make<TH1D>(Form("pl%i_hit_sigmaint",i),"SigmaIntegral / Integral", 200, 0, 2);
      h_hit_adcdiff[i]       =dir_diag.make<TH1D>(Form("pl%i_adcdiff",i),"(Integral - SumADC) / SumADC", 200, -1, 1);

      h_nclusts[i] = dir_diag.make<TH1D>(Form("pl%i_nclusts",i),Form("nclusts, plane %i",i),1000,0,1000);
      h_nclusts_pm[i] = dir_diag.make<TH1D>(Form("pl%i_nclusts_planematched",i),Form("nclusts plane matched, plane %i",i),1000,0,1000);
      
      
      h_chargecomp[i] = dir_truth.make<TH1D>(Form("pl%i_hit_charge_completeness",i),Form("charge completness, plane %i",i),101,0,1.01);
      //h_hitpur[i]     = dir_truth.make<TH1D>(Form("pl%i_hit_purity",i),Form("hit purity, plane %i",i),101,0,1.01);
      h_hitqres[i] = dir_truth.make<TH1D>( Form("pl%i_hit_qres",i),Form("Plane %i;hit charge resolution: (reco-true)/true",i),300,-1.5,1.5);
      h_hitqres_scatter[i] = dir_truth.make<TH2D>( Form("pl%i_hit_qres_scatter",i),
        Form("Plane %i;true hit charge [ #times 10^{3} electrons ];Reconstructed hit charge [ #times 10^{3} electrons ]",i),160,0,80, 160,0,80);
        h_hitqres_scatter[i] ->SetOption("colz");
      h_hitqres_vs_q[i] = dir_truth.make<TH2D>( Form("pl%i_hit_qres_vs_q",i),
        Form("Plane %i;true hit charge [ #times 10^{3} electrons ];hit charge resolution: (reco-true)/true",i),320,0,80, 300,-1.5,1.5);
        h_hitqres_vs_q[i] ->SetOption("colz");
      //h_hitqres_mip[i]            = dir_truth.make<TH1D>(Form("pl%i_hit_qres_mip",i),            Form("Plane %i hits in cosmic #mu tracks;(reco-true)/true",i),                       300,-1.5,1.5);
      //h_hitqres_iso_true[i]       = dir_truth.make<TH1D>(Form("pl%i_hit_qres_iso_true",i),       Form("Plane %i hits (isolated, no plane match, truth-matched);(reco-true)/true",i),  300,-1.5,1.5);
      //h_hitqres_isomatch_true[i]  = dir_truth.make<TH1D>(Form("pl%i_hit_qres_isomatch_true",i),  Form("Plane %i hits (isolated, plane-matched, truth-matched);(reco-true)/true",i),   300,-1.5,1.5);
      //h_hitqres_alpha[i]          = dir_truth.make<TH1D>(Form("pl%i_hit_qres_alpha",i),          Form("Plane %i alpha hits;(reco-true)/true",i),                                       300,-1.5,1.5);
      

    }//endloop over planes
   
    
    // Cross-check that diffusion is being simulated
    h_hit_tdrift_vs_RMS = dir_truth.make<TH2D>("hit_tdrift_vs_RMS","Collection plane;Charge drift time [#mus];Hit RMS [ADC]",230,0,2300,100,0,10);
    h_hit_tdrift_vs_RMSratio = dir_truth.make<TH2D>("hit_tdrift_vs_RMSratio","Collection plane;Charge drift time [#mus];Hit RMS/amplitude ratio",230,0,2300,100,0,5);

 }

};//class BlipAna


//###################################################
//  BlipAna constructor and destructor
//###################################################
BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData  (nullptr)
  ,fBlipAlg ( pset.get<fhicl::ParameterSet>("BlipAlg") )
{
  // blip reconstruction algorithm class
  fhicl::ParameterSet pset_blipalg = pset.get<fhicl::ParameterSet>("BlipAlg");
  fHitProducer    = pset_blipalg.get<art::InputTag> ("HitProducer","");
  fTrkProducer    = pset_blipalg.get<std::string>   ("TrkProducer",     "pandora");
  fGeantProducer  = pset_blipalg.get<std::string>   ("GeantProducer",   "largeant");
  fSimDepProducer = pset_blipalg.get<std::string>   ("SimEDepProducer", "ionization");
  fCaloPlane      = pset_blipalg.get<int>           ("CaloPlane",       2);
  fSavePlaneInfo  = pset.get<std::vector<bool>>     ("SavePlaneInfo",   {true,true,true});
  fDebugMode      = pset.get<bool>                  ("DebugMode",       false);
  fDoACPTrkCalib  = pset.get<bool>                  ("DoACPTrkCalib",   true);
  fACPTCaliProducer = pset.get<art::InputTag>       ("ACPTCaliProducer","pandoracaliInit");

  fPFPproducer = pset.get<art::InputTag>("PFPproducer","pandora");
  fSHRproducer = pset.get<art::InputTag>("SHRproducer","shrreco3d");
  fHITproducer = pset.get<art::InputTag>("HITproducer","pandora");
  fVTXproducer = pset.get<art::InputTag>("VTXproducer","pandora");
  fPCAproducer = pset.get<art::InputTag>("PCAproducer","pandora");
  fCLSproducer = pset.get<art::InputTag>("CLSproducer","pandora");
  fSLCproducer = pset.get<art::InputTag>("SLCproducer","pandora");
  fMCTproducer = pset.get<art::InputTag>("MCTproducer","generator");
  fTRKproducer = pset.get<art::InputTag>("TRKproducer","pandora");


  // data tree object
  fData = new BlipAnaTreeDataStruct();
  fData ->treeName        = pset.get<std::string> ("EventTreeName", "anatree");
  fData ->saveTruthInfo   = pset.get<bool>        ("SaveTruthInfo", true);
  fData ->saveTrueParticles = pset.get<bool>      ("SaveTrueParticles", true);
  fData ->savePrimaries     = pset.get<bool>      ("SavePrimaries", false);
  fData ->saveTrueEDeps = pset.get<bool>          ("SaveTrueEDeps", true);
  fData ->saveTrkInfo     = pset.get<bool>        ("SaveTrkInfo",   true);
  fData ->saveHitInfo     = pset.get<bool>        ("SaveHitInfo",   true);
  fData ->saveClustInfo   = pset.get<bool>        ("SaveClustInfo", true);
  fData ->saveNuInfo      = pset.get<bool>        ("SaveNeutrinoInfo",true);
  fData ->Clear();
  fData ->MakeTree();
  if( fDoACPTrkCalib ) fData->MakeCalibTree();

  // initialize histograms
  InitializeHistograms();
    
  /*
  //==================================================
  // Map out the space charge effects
  //==================================================
  auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* SCE     = reinterpret_cast<spacecharge::SpaceChargeMicroBooNE const*>(lar::providerFrom<spacecharge::SpaceChargeService>());
  TH2D* h_xz_num  = h_efield_distortion_xz;
  TH2D* h_yz_num  = h_efield_distortion_yz;
  TH2D* h_xz      = (TH2D*)h_efield_distortion_xz->Clone("xz");
  TH2D* h_yz      = (TH2D*)h_efield_distortion_yz->Clone("yz");
  float Efield = detProp->Efield();
  int xbins = h_xz->GetYaxis()->GetNbins();
  int ybins = h_yz->GetYaxis()->GetNbins();
  int zbins = h_yz->GetXaxis()->GetNbins();
  
  for(int iz=1; iz<=zbins; iz++){
    float z = h_yz->GetXaxis()->GetBinCenter(iz);
    
    for(int ix=1; ix<=xbins; ix++){
      float x = h_xz->GetYaxis()->GetBinCenter(ix);
      


      for(int iy=1; iy<=ybins; iy++){
        float y = h_yz->GetYaxis()->GetBinCenter(iy);
      
        // fiducialize
        //if( x < 100 || x > 105 ) continue;
        //if( y < -80 || y > 80 ) continue;
        //if( z < 50  || z > 985 ) continue;

        auto const field_offset = SCE->GetCalEfieldOffsets(geo::Point_t{x,y,z});
        float EfieldMod = Efield * std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());
        h_yz_num  ->Fill(z,y, EfieldMod/Efield-1);
        h_yz      ->Fill(z,y);
        h_xz_num  ->Fill(z,x, EfieldMod/Efield-1);
        h_xz      ->Fill(z,x);
      
      }
    }
  }
  
  h_yz_num->Divide(h_yz);
  h_xz_num->Divide(h_xz);
  */
  
}

BlipAna::~BlipAna(){
  delete fData;
}



//###################################################
//  Main event-by-event analysis
//###################################################
void BlipAna::analyze(const art::Event& evt)
{ 
  //============================================
  // New event!
  //============================================
  fData            ->Clear();
  fData->event      = evt.id().event();
  fData->run        = evt.id().run();
  fData->subrun     = evt.id().subRun();
  fIsRealData       = evt.isRealData();
  fNumEvents++;

  // Get timestamp
  unsigned long long int tsval = evt.time().value();
  const unsigned long int mask32 = 0xFFFFFFFFUL;
  fData->timestamp = ( tsval >> 32 ) & mask32;

  // Retrieve lifetime
  const lariov::UBElectronLifetimeProvider& elifetime_provider 
    = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
  float electronLifetime = elifetime_provider.Lifetime() * /*convert ms->mus*/ 1e3;
  fData->lifetime = electronLifetime;
  
  auto const  detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
  auto const& SCE       = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const& tpcCalib  = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
  

  //============================================
  // Run blip reconstruction: 
  //============================================
  //  
  //  In this step, we pass the entire art::Event to the algorithm, 
  //  and it creates a single collection of blip 'objects', a special data
  //  struct in the 'blip' namespace defined in BlipUtils.h.
  //  
  //  We can then retrieve these blips and incorporate them into
  //  our analysis however we like:
  //
  //    std::vector<blipobj::Blip> blipVec = fBlipAlg.blips;
  //
  //  The alg also creates collections of 'HitInfo' and 'HitClust'
  //  structs used in the blip reconstruction process, which can be
  //  accessed in the same way as blips. 
  //    
  //    * HitInfo simply saves some calculations for each hit that aren't 
  //      present in the native recob::Hit object, like drift time, associated 
  //      G4 particle IDs, etc.
  //
  //    * HitClust is just a cluster of hits on a specific plane; these are 
  //      used to create 3D blips by plane-matching.
  //
  fBlipAlg.RunBlipReco(evt);


  //=======================================
  // Get data products for this event
  //========================================

  // -- MCTruth 
  art::Handle< std::vector<simb::MCTruth> > truthHandle;
  std::vector<art::Ptr<simb::MCTruth> > truthlist;
  if (evt.getByLabel("generator",truthHandle))
    art::fill_ptr_vector(truthlist, truthHandle); 
 
  // -- G4 particles
  art::Handle< std::vector<simb::MCParticle> > pHandle;
  std::vector<art::Ptr<simb::MCParticle> > plist;
  if (evt.getByLabel("largeant",pHandle))
    art::fill_ptr_vector(plist, pHandle);
  
  // -- hits (from input module)
  fHitProducer = fBlipAlg.fHitProducer;
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitProducer,hitHandle))
    art::fill_ptr_vector(hitlist, hitHandle);

  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrkProducer,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
 
  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = std::min((int)plist.size(),(int)kMaxG4);
  fData->ntrks      = (int)tracklist.size();
  fData->badchans   = fBlipAlg.EvtBadChanCount;
  fData->Resize();
 
  // flag this data as MC
  fIsMC = ( plist.size()>0 );
 


  //=======================================================
  // Check if the event contains neutrino information
  // (much of this copied from ubana/searchingfornues)
  //=======================================================
  // -- PFPs
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfplist;
  if (evt.getByLabel(fPFPproducer,pfpHandle))
    art::fill_ptr_vector(pfplist, pfpHandle);

  // -- associated tracks/vertex
  art::FindManyP<recob::Track> fmtrk_from_pfp(pfpHandle,evt,fTrkProducer);
  if( pfplist.size() ) {
    // grab PFParticles in event
    ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(evt, fPFPproducer,
      proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
      proxy::withAssociated<recob::Cluster>(fCLSproducer),
      proxy::withAssociated<recob::Slice>(fSLCproducer),
      proxy::withAssociated<recob::Track>(fTRKproducer),
      proxy::withAssociated<recob::Vertex>(fVTXproducer),
      proxy::withAssociated<recob::PCAxis>(fPCAproducer),
      proxy::withAssociated<recob::Shower>(fSHRproducer),
      proxy::withAssociated<recob::SpacePoint>(fPFPproducer));
    BuildPFPMap(pfp_proxy);
    // loopthrough PFParticles
    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
    {
      // get metadata for this PFP
      const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

      //  find neutrino candidate
      if (pfp_pxy->IsPrimary() == false) continue;
      auto PDG = fabs(pfp_pxy->PdgCode());
      fData->nu_pfp_pdg=PDG;
      if ( (PDG == 12) || (PDG == 14) )
      {
        //std::cout<<"Found a neutrino PFP\n";
        if (pfParticleMetadataList.size() != 0)
        {
          for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
          {
            const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
            auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
            if (!pfParticlePropertiesMap.empty())
            {
              for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
              {
                if( it->first == "IsNeutrino" ) fData->nu_isNeutrino  = it->second;
                if( it->first == "NuScore"    ) fData->nu_nuscore.push_back(it->second);
              }
            }
          }
        } // if PFP metadata exists!
  

        // Get vertex info
        double xyz[3] = {};
        auto vtx = pfp_pxy.get<recob::Vertex>();
        if (vtx.size() == 1)
        {
          // save vertex to array
          vtx.at(0)->XYZ(xyz);
          auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);
          float _reco_nu_vtx_sce[3];
          nuselection::ApplySCECorrectionXYZ(nuvtx.X(),nuvtx.Y(),nuvtx.Z(), _reco_nu_vtx_sce);
          fData->nu_reco_vtx_x = _reco_nu_vtx_sce[0];
          fData->nu_reco_vtx_y = _reco_nu_vtx_sce[1];
          fData->nu_reco_vtx_z = _reco_nu_vtx_sce[2];
          //std::cout<<"Vertex:  "<<nuvtx.X()<<"  "<<nuvtx.Y()<<"  "<<nuvtx.Z()<<"\n";
        }
        else
        {
          std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
        }

        // collect PFParticle hierarchy originating from this neutrino candidate
        std::vector<ProxyPfpElem_t> slice_pfp_v;
        AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
        //std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl;
        // create list of tracks and showers associated to this slice
        std::vector<float>  trkscore_v;
        std::vector<int>    trkid_v;
        std::vector<float>  shwrscore_v;

        for (auto pfp : slice_pfp_v)
        { 
          auto const &ass_trk_v = pfp.get<recob::Track>();
          auto const &ass_shr_v = pfp.get<recob::Shower>();
          float score = nuselection::GetTrackShowerScore(pfp);
          if (ass_trk_v.size() == 1) {
            trkid_v   .push_back(ass_trk_v.at(0)->ID());
            trkscore_v.push_back(score);
          }
          if (ass_shr_v.size() == 1) {
            shwrscore_v.push_back(score);
          }
        } // for all PFParticles in the slice
        //std::cout<<"  - "<<trkscore_v.size()<<" tracks\n";
        //std::cout<<"  - "<<shwrscore_v.size()<<" showers\n";
        for(size_t itrk = 0; itrk < trkscore_v.size(); itrk++){
          fData->nu_trk_id    .push_back(trkid_v[itrk]);
          fData->nu_trk_score .push_back(trkscore_v[itrk]);
        }
        for(size_t ishwr = 0; ishwr < shwrscore_v.size(); ishwr++){
          fData->nu_shwr_score .push_back(shwrscore_v[ishwr]);
        }
      
      }//if PDG of neutrino
    }//end loop over PFPs
  }
  
  
  // Tell us what's going on!
  if( fNumEvents < 200 || (fNumEvents % 100) == 0 ) {
  std::cout<<"\n"
  <<"=========== BlipAna =========================\n"
  <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"; total: "<<fNumEvents<<"\n";
  std::cout
  <<"found "<<hitlist.size()<<" hits from "<<fBlipAlg.fHitProducer<<"\n"
  <<"found "<<tracklist.size()<<" tracks from "<<fTrkProducer<<"\n"
  <<"found "<<pfplist.size()<<" PFPs from "<<fPFPproducer<<"\n"
  <<"reconstructed "<<fBlipAlg.hitclust.size()<<" 2D hit clusters\n"
  <<"reconstructed "<<fBlipAlg.blips.size()<<" 3D blips\n"
  ;
  }





  //====================================
  // Keep tabs on total energy, charge,
  // and electrons arriving at anode
  //====================================
  float total_depEnergy         = 0;
  float total_depElectrons      = 0;
  float total_numElectrons      = 0;

  //===================================
  // Check neutrinos in MCTruth
  // (NuanceOffset found through simb::kNuanceOffset)
  //===================================
  if( truthlist.size() > 0 ) {
    auto& nu    = truthlist[0]->GetNeutrino();
    auto& part  = truthlist[0]->GetParticle(0);
    auto PDG    = fabs(part.PdgCode());
    if ( (PDG == 12) || (PDG == 14) ) {
      fData->mctruth_nu_pdg   = part.PdgCode();
      fData->mctruth_nu_ccnc  = nu.CCNC();
      fData->mctruth_nu_mode  = nu.Mode();
      fData->mctruth_nu_vtx_x = part.EndPosition()[0];
      fData->mctruth_nu_vtx_y = part.EndPosition()[1];
      fData->mctruth_nu_vtx_z = part.EndPosition()[2];
      fData->mctruth_nu_KE    =  /*GeV->MeV*/1e3 * (part.E()-part.Mass());
    }
  }



  //====================================
  // Save MCParticle information
  //====================================
  std::map<int,int> map_g4trkid_index;
  if( plist.size() ) {
    
    std::vector<blipobj::ParticleInfo>& pinfo = fBlipAlg.pinfo;
    
    // Loop through the MCParticles
    if( fDebugMode ) std::cout<<"\nLooping over "<< plist.size()<<" G4 MCParticles: \n";
    int printed = 0;
    for(size_t i = 0; i<plist.size(); i++){
      auto& pPart = plist[i];
      map_g4trkid_index[pPart->TrackId()] = i;
      total_depEnergy       += pinfo[i].depEnergy;
      total_depElectrons    += pinfo[i].depElectrons;
      
      // Save to TTree object
      if(i<kMaxG4){
        
        if( fData->savePrimaries && pinfo[i].isPrimary ) {
          fData->nprimaries++;
          fData->primary_g4id[i] = pPart->TrackId();
          fData->primary_pdg[i]     = pPart->PdgCode();
          fData->primary_Px[i]      = pinfo[i].Px;
          fData->primary_Py[i]      = pinfo[i].Py;
          fData->primary_Pz[i]      = pinfo[i].Pz;
          fData->primary_x0[i]      = pPart->Vx();
          fData->primary_y0[i]      = pPart->Vy();
          fData->primary_z0[i]      = pPart->Vz();
          fData->primary_T0[i]      = pinfo[i].time;
          if( pinfo[i].pathLength ) {
            fData->primary_xAV[i]     = pinfo[i].startPoint.X();
            fData->primary_yAV[i]     = pinfo[i].startPoint.Y();
            fData->primary_zAV[i]     = pinfo[i].startPoint.Z();
          }
        }

        fData->part_g4id[i]         = pPart->TrackId();
        fData->part_pdg[i]             = pPart->PdgCode();
        fData->part_nDaughters[i]      = pPart->NumberDaughters();
        fData->part_mother[i]          = pPart->Mother();
        fData->part_E[i]               = pinfo[i].E;
        fData->part_endE[i]            = pinfo[i].endE;
        fData->part_mass[i]            = pinfo[i].mass;
        fData->part_KE[i]              = pinfo[i].KE;
        fData->part_endKE[i]           = pinfo[i].endKE;
        fData->part_P[i]               = pinfo[i].P;
        fData->part_Px[i]              = pinfo[i].Px;
        fData->part_Py[i]              = pinfo[i].Py;
        fData->part_Pz[i]              = pinfo[i].Pz;
        fData->part_startPointx[i]     = pPart->Vx();
        fData->part_startPointy[i]     = pPart->Vy();
        fData->part_startPointz[i]     = pPart->Vz();
        fData->part_endPointx[i]       = pPart->EndPosition()[0];
        fData->part_endPointy[i]       = pPart->EndPosition()[1];
        fData->part_endPointz[i]       = pPart->EndPosition()[2];
        fData->part_startT[i]          = pinfo[i].time;
        fData->part_endT[i]            = pinfo[i].endtime;
        fData->part_pathlen[i]         = pinfo[i].pathLength;
        fData->part_numTrajPts[i]      = pinfo[i].numTrajPts;
        fData->part_process[i]         = pPart->Process();
        fData->part_depEnergy[i]       = pinfo[i].depEnergy;
        fData->part_depElectrons[i]    = pinfo[i].depElectrons;
        fData->part_isPrimary[i]       = pinfo[i].isPrimary;
        // check containment
        bool startInAV  = BlipUtils::IsPointInAV(pPart->Vx(),pPart->Vy(),pPart->Vz());
        bool endInAV    = BlipUtils::IsPointInAV(pPart->EndPosition()[0],pPart->EndPosition()[1],pPart->EndPosition()[2]);
        if( startInAV && endInAV ) fData->part_isContained[i] = true;
 
        if( fDebugMode && pPart->Process() != "muIoni" ) {
          PrintParticleInfo(i);
          printed++;
        }
      }
    } // endloop over G4 particles
  
    if( fDebugMode ) std::cout<<"True total energy deposited: "<<total_depEnergy<<" MeV \n";
  
  }//endif particles found in event
  

  //====================================
  // Save TrueBlip information
  //====================================
  std::vector<blipobj::TrueBlip>& trueblips = fBlipAlg.trueblips;
  fData->nedeps = (int)trueblips.size();
  if( trueblips.size() ) {
    if( fDebugMode ) std::cout<<"\nLooping over "<<trueblips.size()<<" true blips:\n";
    for(auto& trueblip : trueblips ) {
      if( fDebugMode ) PrintTrueBlipInfo(trueblip);
      int i     = trueblip.ID;
      int ig4   = trueblip.LeadG4Index;
      auto& pPart = plist[ig4];
      fData->edep_tpc[i]      = trueblip.TPC;
      fData->edep_energy[i]   = trueblip.Energy;
      fData->edep_electrons[i]= trueblip.DepElectrons;
      //fData->edep_charge[i]   = trueblip.NumElectrons;
      fData->edep_x[i]        = trueblip.Position.X();
      fData->edep_y[i]        = trueblip.Position.Y();
      fData->edep_z[i]        = trueblip.Position.Z();
      fData->edep_tdrift[i]   = trueblip.DriftTime;
      fData->edep_pdg[i]      = trueblip.LeadG4PDG;
      fData->edep_g4id[i]     = trueblip.LeadG4ID;
      fData->edep_g4qfrac[i]  = trueblip.G4ChargeMap[trueblip.LeadG4ID] / trueblip.DepElectrons;
      fData->edep_dz[i]       = fabs(pPart->EndPosition()[2]-pPart->Vz());
      fData->edep_dx[i]       = fabs(pPart->EndPosition()[0]-pPart->Vx());
      fData->edep_allchansgood[i] = trueblip.AllChansGood; 

      int proc_code = -9;
      std::string proc = fData->part_process[ig4];
      if(     proc == "primary") { h_part_process->Fill("primary",1); proc_code = 0; }
      else if(proc == "compt"  ) { h_part_process->Fill("compt",  1); proc_code = 1; }
      else if(proc == "phot"   ) { h_part_process->Fill("phot",   1); proc_code = 2; }
      else if(proc == "conv"   ) { h_part_process->Fill("conv",   1); proc_code = 3; }
      else                       { h_part_process->Fill("other",  1); proc_code = 4; }
      fData->edep_proc[i] = proc_code;
      //float ne_dep  = trueblip.DepElectrons;
      float ne      = trueblip.NumElectrons;
      total_numElectrons += ne;
      // calculate simulated lifetime
      //if( ne>10000 && ne<ne_dep && trueblip.DriftTime > 1000. ) {
      //  float tau = trueblip.DriftTime / log(ne_dep/ne);
      //  h_sim_lifetime->Fill(tau);
      //}
    }
  }//endif trueblips were made
  
  float driftVelocity = detProp.DriftVelocity(detProp.Efield(),detProp.Temperature());
  
  //====================================
  // Save track information
  //====================================
  if( fDebugMode ) std::cout<<"Looping over "<<tracklist.size()<<" tracks...\n";
  std::map<int,bool> map_trkid_isMIP;
  std::map<int,float> map_trkid_length;
  int trks_100cm[2] = {0, 0};
  fData->longtrks=0;
  for(size_t i=0; i<tracklist.size(); i++){
    auto& trk = tracklist[i];
    const auto& startPt = trk->Vertex();
    const auto& endPt   = trk->End();
    bool isMC = fBlipAlg.map_trkid_isMC[trk->ID()];
    fData->trk_isMC[i]  = isMC;
    fData->trk_g4id[i]  = fBlipAlg.map_trkid_g4id[trk->ID()];
    fData->trk_id[i]    = trk->ID();
    fData->trk_npts[i]  = trk->NumberTrajectoryPoints();
    fData->trk_length[i]= trk->Length();
    fData->trk_startx[i]= startPt.X();
    fData->trk_starty[i]= startPt.Y();
    fData->trk_startz[i]= startPt.Z();
    fData->trk_endx[i]  = endPt.X();
    fData->trk_endy[i]  = endPt.Y();
    fData->trk_endz[i]  = endPt.Z();
    fData->trk_startd[i]= BlipUtils::DistToBoundary(startPt);
    fData->trk_endd[i]  = BlipUtils::DistToBoundary(endPt);
    h_trk_length  ->Fill(trk->Length());
    float dX = fabs(startPt.X() - endPt.X());
    float dY = fabs(startPt.Y() - endPt.Y());
    h_trk_xspan   ->Fill( dX );
    h_trk_yspan   ->Fill( dY );
    map_trkid_length[trk->ID()] = trk->Length();
    map_trkid_isMIP[trk->ID()]  = (trk->Length()>100) ? true : false;
    // count the number of non-blippy tracks to use
    // as a metric for cosmic activity in event
    if( trk->Length() > 5 ) fData->longtrks++;
    if( trk->Length() > 100 ) trks_100cm[(int)isMC]++;
    // TODO:
    // identify "cosmic"-looking tracks that pierce
    // the top of the TPC ceiling at Y = + 117cm
    //float ytop = 117;
    //bool isStartAtBnd = ( fabs(startPt.Y()-ytop) < 2. );
    //bool isEndAtTop   = ( fabs(endPt.Y()  -ytop) < 2. );
    //if( (isStartAtTop || isEndAtTop) && dY > 20. ) {
    //}
  }
  
  h_trks_100cm[0]->Fill(trks_100cm[0]);
  h_trks_100cm[1]->Fill(trks_100cm[1]);

    

  //================================================ 
  // In-house lifetime calibration using anode-to-
  // cathode-piercing (ACP) tracks. 
  //================================================
  if( fDoACPTrkCalib ) {
    // retrieve track calo data product
    // - pandoracalo = no corrections; 
    // - pandoracali = YZ transparency corrections only
    art::FindManyP<anab::Calorimetry> fmcal(tracklistHandle, evt, fACPTCaliProducer);
    
    if( fmcal.isValid() ) {
      
      //std::cout<<"We are doing ACPT\n";
      // set the required 'dX' that would indicate a
      // track crossed the full drift distance
      float dx_min = 250;
      float dx_max = 270;
      //================================================
      // loop over the tracklist 
      //================================================
      for(size_t i=0; i<tracklist.size(); i++){
        auto& trk = tracklist[i];
        auto& startPt = trk->Vertex();
        auto& endPt   = trk->End();
        float dX      = fabs(startPt.X() - endPt.X());
        float dY      = fabs(startPt.Y() - endPt.Y());
        float dZ      = fabs(startPt.Z() - endPt.Z());
        //std::cout<<"trk dX "<<dX<<"\n";
        // skip tracks that don't look like anode-cathode-piercing
        if( dX < dx_min || dX > dx_max ) continue;
        // select only good angles (exclude 90 degree +/- 5 deg)
        float theta_xz = 180*atan(dX/dZ)/3.14159;
        float theta_yz = 180*atan(dY/dZ)/3.14159;
        h_ACPtrk_theta_xz->Fill(theta_xz);
        h_ACPtrk_theta_yz->Fill(theta_yz);
        //std::cout<<"trk angle "<<theta_xz<<"\n";
        if( theta_xz > 80 ) continue;
      
        // Make sure this isn't just a really long wiggly track
        //   1 = impossibly wiggly
        //   0 = perfectly straight
        float diff = sqrt( pow(dX,2) + pow(dY,2) + pow(dZ,2) );
        float wigglyness = 1.-diff/trk->Length();
        if(wigglyness > 0.01 ) continue;
        //std::cout<<"wigglyness "<<wigglyness<<"\n";

        //==============================================
        // loop the calo objects and find the one for collection plane
        // also require at least 100 dE/dx points
        //==============================================
        //std::cout<<"LOOKING FOR TRACK CALO .....................................................................\n";
        std::vector<art::Ptr<anab::Calorimetry> > caloObjs = fmcal.at(i);
        for(auto& caloObj : caloObjs ) {
          if( caloObj->PlaneID().Plane != 2 )       continue;
          if( caloObj->dEdx().size() < 100 )        continue;
          size_t Npts = caloObj->dEdx().size();
          //std::cout<<"calo NPts "<<Npts<<"\n";

          // Calculate the X-coordinate offset that we'll use to shift
          // the track so that its start/end correspond to anode and cathode
          //  - in MicroBooNE, X starts at 0 and increases toward cathode,
          //    so we take the MIN of the two track endpoints as 0
          float X0 = std::min( startPt.X(), endPt.X() );
         
          // First, apply SCE spatial corrections at every 
          // 3D point along the track...
          std::vector<geo::Point_t> vec_XYZ_orig;
          std::vector<geo::Point_t> vec_XYZ_corr;
          for(size_t j=0; j<Npts; j++){
            geo::Point_t point{ caloObj->XYZ().at(j).X() - X0, caloObj->XYZ().at(j).Y(), caloObj->XYZ().at(j).Z() };
            vec_XYZ_orig.push_back( point );
            geo::Vector_t loc_offset = SCE->GetCalPosOffsets( point, 0 );
            point.SetXYZ(point.X()-loc_offset.X(),point.Y()+loc_offset.Y(),point.Z()+loc_offset.Z());
            vec_XYZ_corr.push_back( point );
          }
          
          // record the "dx correction factor" applied at each
          // trajectory point along the track
          std::vector<float> dx_corr_factor( Npts, 1. );
          for(size_t j=1; j<Npts-1; j++){
            float ds_prev = sqrt( ( vec_XYZ_orig[j-1] - vec_XYZ_orig[j+1] ).Mag2() );
            float ds_new  = sqrt( ( vec_XYZ_corr[j-1] - vec_XYZ_corr[j+1] ).Mag2() );
            float corr_factor = ds_new/ds_prev;
            dx_corr_factor[j] = corr_factor;
          }
          
          // ------------------------------------
          // loop over track points
          // ------------------------------------
          size_t nSavedPts = 0;
          //FillWith(fData->acptrk_dEdx,   0);
          //FillWith(fData->acptrk_tdrift, 0);
          std::vector<float> vec_dEdx_near;
          std::vector<float> vec_dEdx_far;
          
          for(size_t j=1; j<Npts-1; j++){
            
            // make sure there are neighboring wires
            float locZ  = vec_XYZ_corr.at(j).Z();
            float prevZ = vec_XYZ_corr.at(j-1).Z();
            float nextZ = vec_XYZ_corr.at(j+1).Z();
            if( fabs(locZ-prevZ) > 0.5 || fabs(locZ-nextZ) > 0.5 ) continue;

            // correct for SCE field offset
            auto& point = vec_XYZ_corr.at(j);
            float Efield = detProp.Efield();
            if( SCE->EnableCalEfieldSCE() ) {
              auto const field_offset = SCE->GetCalEfieldOffsets(point, 0); 
              Efield = detProp.Efield()*std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());;
            }
            if( Efield <= 0 ) continue;
            
            // Calculate new dE/dx
            float dQdx_ADC  = caloObj->dQdx().at(j) / dx_corr_factor.at(j);
            float dQdx_e    = fBlipAlg.fCaloAlg.ElectronsFromADCArea(dQdx_ADC,2);
            float dEdx      = fBlipAlg.dQdx_to_dEdx( dQdx_e, Efield );
            float tdrift    = point.X() / driftVelocity;

            // Quality control
            if( dEdx < 1.0 || dEdx > 20. || tdrift > 2400 ) continue;
            
            // Save this point!
            nSavedPts++;
            h_ACPtrk_yz   ->Fill( point.Z(), point.Y() );
            h_ACPtrk_dEdx ->Fill(dEdx);
          
            // Divide dEdx into those near anode or cathode
            if( fabs(point.X()-30.) < 10. ) {
              vec_dEdx_near       .push_back( dEdx );
              h_ACPtrk_dEdx_near  ->Fill( dEdx );
            } else if( fabs(point.X()-230) < 10. ) {
              vec_dEdx_far        .push_back( dEdx );
              h_ACPtrk_dEdx_far   ->Fill( dEdx );
            }
          
          }//endloop over track points

          // Median of near and far zones (only calculate if there are
          // at least 15 data points in each of them)
          float qratio_median = -9; 
          float qratio_mean   = -9;
          if( vec_dEdx_near.size() >= 15 && vec_dEdx_far.size() >= 15 ) {
            float median_near = BlipUtils::FindMedian(vec_dEdx_near);
            float median_far  = BlipUtils::FindMedian(vec_dEdx_far);
            float mean_near   = BlipUtils::FindMean(vec_dEdx_near);
            float mean_far    = BlipUtils::FindMean(vec_dEdx_far);
            qratio_median     = median_far / median_near;
            qratio_mean       = mean_far / mean_near;
            h_ACPtrk_qratio   ->Fill(qratio_median);
          
            // Fill calib tree
            fData->acptrk_npts            = nSavedPts;
            fData->acptrk_qratio_median   = qratio_median;
            fData->acptrk_qratio_mean     = qratio_mean;
            fData->acptrk_theta_xz        = theta_xz;
            fData->acptrk_theta_yz        = theta_yz;
            fData->calibTree              ->Fill();
          }
        
        }//endloop on calo objs
      }//end tracklist loop
    }//endIF calo objects are valid
  }//endIF do ACPTrkCalib
  // DONE with ACP track lifetime calibrations
  
  
  


  //====================================
  // Save hit information
  //====================================
  if( fDebugMode ) std::cout<<"Looping over the hits...\n";
  int   num_hits[kNplanes]            ={0};
  int   num_hits_untracked[kNplanes]  ={0};
  int   num_hits_true[kNplanes]       ={0};
  int   num_hits_pmatch[kNplanes]     ={0};
  float total_hit_charge[kNplanes]    ={0};
  
  for(size_t i=0; i<hitlist.size(); i++){
   
    int     plane   = hitlist[i]->WireID().Plane;
    int     ndf     = hitlist[i]->DegreesOfFreedom();
    double  gof     = (ndf>0) ? hitlist[i]->GoodnessOfFit()/ndf : -9;
    int     isFit   = int(gof>=0);
    int     mult    = (isFit) ? hitlist[i]->Multiplicity() : -9;
    float   logGOF  = (isFit && mult == 1) ? log10(gof) : 99;
    float   amp     = (isFit && mult == 1) ? hitlist[i]->PeakAmplitude() : -9;
    float   rms     = (isFit && mult == 1) ? hitlist[i]->RMS() : -9; 
    float sumADC = hitlist[i]->ROISummedADC();
    float integral = hitlist[i]->Integral();
    
    auto const& hinfo = fBlipAlg.hitinfo[i];
    bool    isMC    = (hinfo.g4trkid >= 0 );
    bool    isTrked = (hinfo.trkid >= 0 && map_trkid_length[hinfo.trkid] > 5 ); 
    bool    isMIP     = map_trkid_isMIP[hinfo.trkid];
    bool    isMatched = hinfo.ismatch; 
    bool    isAlpha     = (hinfo.g4pdg == 1000020040);
    bool    isElectron  = (hinfo.g4pdg == 11);
    bool    isPrimary = true;
    int     g4index = -9;
    if( hinfo.g4trkid >= 0 ) {
      g4index = map_g4trkid_index[hinfo.g4trkid];
      isPrimary = fBlipAlg.pinfo[g4index].isPrimary;
    }

    fNumHits[plane]++;
    num_hits[plane]++;
    
    // calculate reco-true resolution
    float qres      = -9;
    if( isMC ) {
      float qtrue = hinfo.g4charge;
      float qcoll = hinfo.charge;
      if(qcoll && qtrue){
        qres = (qcoll-qtrue)/qtrue;
        fNumHitsTrue[plane]++;
        if( isMatched ) fNumHitsMatchedTrue[plane]++;
        num_hits_true[plane]++;
        total_hit_charge[plane] += qtrue;
        h_hitqres_scatter[plane]->Fill(qtrue/1e3,qcoll/1e3);
        h_hitqres_vs_q[plane]->Fill(qtrue/1e3,(qcoll-qtrue)/qtrue);
      }

      // Find associated EDep
      // blah
      if( plane==2 && rms > 0 && amp > 0 ) {
        for(auto& trueblip : fBlipAlg.trueblips ) {
          if( trueblip.LeadG4ID == hinfo.g4trkid ) {
            h_hit_tdrift_vs_RMS->Fill(trueblip.DriftTime,rms);
            h_hit_tdrift_vs_RMSratio->Fill(trueblip.DriftTime,rms/amp);
          }
        }
      }
    }

    // *** untracked hits (not in tracks > 5cm in length) ***
    if( !isTrked ) {

      fNumHitsUntracked[plane]++;
      num_hits_untracked[hinfo.plane]++;
      h_hitamp[plane]   ->Fill(amp);
      h_hitrms[plane]   ->Fill(rms);
      h_hitmult[plane]  ->Fill(mult);
      h_hitqres[plane]  ->Fill(qres);
      h_hit_charge_vs_rms[plane]  ->Fill(hinfo.charge/1e3,rms);
      h_hit_charge_vs_gof[plane]->Fill(hinfo.charge/1e3,logGOF);
     
      // -- untracked and plane-matched --
      if( isMatched ) {
        fNumHitsMatched[plane]++;
        num_hits_pmatch[hinfo.plane]++;
        if( isMC  ) { 
          h_hitmult_isomatch_true[plane]  ->Fill(mult);
          h_hitrms_isomatch_true[plane]   ->Fill(rms); 
        } else { 
          h_hitmult_isomatch[plane]   ->Fill(mult);
          h_hitgof_isomatch[plane]    ->Fill(logGOF);
          h_hitrms_isomatch[plane]    ->Fill(rms); 
        }
      
      // -- untracked and non-plane-matched --
      } else {
        h_hitgof_iso[plane]->Fill(logGOF);
        if( isMC  ) { 
          h_hitmult_iso_true[plane]   ->Fill(mult);
            h_hitrms_iso_true[plane]    ->Fill(rms);
        } else { 
          h_hitmult_iso[plane]  ->Fill(mult);
          h_hitrms_iso[plane]   ->Fill(rms); 
        }
      }

    //*** MIP muon-like tracks
    } else if ( isMIP ) {
      h_hitmult_mip[plane]  ->Fill(mult);
        h_hitamp_mip[plane]   ->Fill(amp);
        h_hitrms_mip[plane]   ->Fill(rms);
        h_hitgof_mip[plane]   ->Fill(logGOF);
    }
    
    // alphas
    if( isAlpha && hinfo.g4frac == 1.0 ) {
      h_hitmult_alpha[plane]  ->Fill(mult);
        h_hitamp_alpha[plane]   ->Fill(amp);
        h_hitrms_alpha[plane]   ->Fill(rms);
        h_hitgof_alpha[plane]   ->Fill(logGOF);
    }
  
    if( isElectron && isPrimary && hinfo.g4frac == 1.0 ) {
      h_hitrms_electron[plane]->Fill(rms);
    }

    if( hitlist[i]->Integral() != 0 ) h_hit_sigmaint[plane] -> Fill( fabs(hitlist[i]->SigmaIntegral()/hitlist[i]->Integral()) );
    if( sumADC != 0 ) h_hit_adcdiff[plane]->Fill( (integral-sumADC)/sumADC );

    // fill data to be saved to event tree
    if( i < kMaxHits && fData->saveHitInfo ){
      fData->hit_plane[i]     = hitlist[i]->WireID().Plane;
      fData->hit_wire[i]      = hitlist[i]->WireID().Wire;
      fData->hit_tpc[i]       = hitlist[i]->WireID().TPC;
      fData->hit_channel[i]   = hitlist[i]->Channel();
      fData->hit_peakT[i]     = hitlist[i]->PeakTime();
      fData->hit_rms[i]       = hitlist[i]->RMS();
      fData->hit_amp[i]	      = hitlist[i]->PeakAmplitude();
      fData->hit_area[i]      = hitlist[i]->Integral();
      fData->hit_sumadc[i]    = hitlist[i]->ROISummedADC();
      fData->hit_mult[i]      = hitlist[i]->Multiplicity();
      fData->hit_gof[i]       = gof;
      fData->hit_trkid[i]     = hinfo.trkid;
      fData->hit_time[i]      = hinfo.driftTime;
      fData->hit_charge[i]    = hinfo.charge;
      fData->hit_ismatch[i]   = hinfo.ismatch;
      fData->hit_g4trkid[i]   = hinfo.g4trkid;
      fData->hit_g4frac[i]    = hinfo.g4frac;
      fData->hit_g4energy[i]  = hinfo.g4energy;
      fData->hit_g4charge[i]  = hinfo.g4charge;
      fData->hit_blipid[i]    = hinfo.blipid;
      fData->hit_clustid[i]   = hinfo.clustid;
    }
  
  }//endloop over hits

  // Now that we've looped all the hits, calculate some
  // plane-specific variables and fill histograms
  for(size_t ip=0; ip<kNplanes; ip++){
    h_nhits[ip]   ->Fill(num_hits[ip]);
    float qcomp     = -9;
    if( num_hits_true[ip] ) {
      if(total_numElectrons )  qcomp = total_hit_charge[ip]/total_depElectrons;
      h_chargecomp[ip]->Fill( qcomp );
    }
  }//endloop over planes
    

 





  /*
  //=======================================================
  // Temporary test of SCE fix implementation
  //=======================================================
  art::Handle< std::vector<anab::Calorimetry> > caloHandle;
  std::vector<art::Ptr<anab::Calorimetry> > calolist;
  if (evt.getByLabel("pandoracaliSCE",caloHandle))
    art::fill_ptr_vector(calolist, caloHandle);
  
  float Efield        = detProp->Efield(0);
  float driftVelocity = detProp->DriftVelocity(Efield,detProp->Temperature());
  
  // loop over the anab::Calorimetry objects
  for(auto& caloObj : calolist ) {
    int plane = caloObj->PlaneID().Plane;
    int Npts  = caloObj->dEdx().size();
    if( plane != 2 ) continue;
    auto& dEdx_vec = caloObj->dEdx();
    auto& dQdx_vec = caloObj->dQdx();
    auto& xyz_vec  = caloObj->XYZ();
    
    // go through the vectors of dQ/dx, dE/dx for this track
    //std::cout<<"  looping over "<<Npts<<" points...\n";
    for(int i=0; i<Npts; i++){
      float x = xyz_vec[i].X();
      float y = xyz_vec[i].Y();
      float z = xyz_vec[i].Z();
	    //float yzcorr = energyCalibProvider.YZdqdxCorrection(plane,y,z);
	    //float xcorr  = energyCalibProvider.XdqdxCorrection(plane,x);
      
      // find E-field offset; skip if 0 (i.e., if outside TPC volume)
      auto const field_offset = SCE->GetCalEfieldOffsets(geo::Point_t{x,y,z});
      if( field_offset.R() == 0 ) continue;
      
      // calcualte recombination using Mod Box, inputting the dE/dx from the calo object
      // and the SCE-modified E-field at this XYZ location
      float EfieldMod = Efield * std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());
      float recomb    = fBlipAlg.ModBoxRecomb(dEdx_vec[i],EfieldMod);

      // what's the predicted dQ/dx?
      float dQdx_pred = (dEdx_vec[i]/(23.6e-6)) * recomb;
    
      // convert dQdx from ADC to electrons
      float dQdx_e = dQdx_vec[i] / (4.1e-3);
      float res = (dQdx_e - dQdx_pred) / dQdx_pred;
     
      h_dQdx_diff->Fill(res);
      h_dEdx_calo->Fill(dEdx_vec[i]);
      if( res > 0.1 ) {
      std::cout <<"XYZ: "<<x<<", "<<y<<", "<<z<<"    Efield "<<EfieldMod
                <<", recomb "<<recomb<<",   lifetime "<<electronLifetime
                <<",   attenuation "<<exp(-T/electronLifetime)<<"\n";
      std::cout <<"tdrift = "<<T<<"    dEdx = "<<dEdx_vec[i]
                <<"   dQdx = "<<dQdx_e<<"   predicted dQdx: "<<dQdx_pred
                <<" --> "<<res<<" fractional diff\n";
      }
    }
  }
  */


      
   
  //=============================================
  // Save hit cluster info
  //=============================================
  fData->nclusts = (int)fBlipAlg.hitclust.size();
  int num_clusts[kNplanes]     ={0};
  int num_clusts_pm[kNplanes]   ={0};
  if( fDebugMode ) std::cout<<"\nLooping over clusters...\n";
  for(size_t i=0; i < fBlipAlg.hitclust.size(); i++){
    auto const& clust = fBlipAlg.hitclust[i];
    num_clusts[clust.Plane]++;
    if( clust.isMatched ) num_clusts_pm[clust.Plane]++;
    if( !fSavePlaneInfo[clust.Plane] ) continue;
    if( i < kMaxClusts ) {
      fData->clust_id[i]        = clust.ID;
      fData->clust_tpc[i]       = clust.TPC;
      fData->clust_plane[i]     = clust.Plane;
      fData->clust_wire[i]        = clust.CenterWire;
      fData->clust_startwire[i]   = clust.StartWire;
      fData->clust_endwire[i]     = clust.EndWire;
      fData->clust_nwires[i]      = clust.NWires;
      fData->clust_nwiresNoisy[i] = clust.NWiresNoisy;
      fData->clust_bydeadwire[i]  = (clust.DeadWireSep==0);
      fData->clust_nhits[i]       = clust.NHits;
      fData->clust_pulsetrain[i] = (clust.NPulseTrainHits>0);
      //fData->clust_time[i]      = clust.Time;
      //fData->clust_charge[i]    = clust.Charge;
      // Truncate precision to reduce file size after ROOT compression
      // (we don't need to know these to the Nth decimal place)
      fData->clust_charge[i]    = Truncate(clust.Charge,    10);
      fData->clust_chargeErr[i] = Truncate(clust.SigmaCharge, 10);
      fData->clust_amp[i]       = Truncate(clust.Amplitude, 0.01);
      fData->clust_time[i]      = Truncate(clust.Time,      0.1);
      fData->clust_timespan[i]  = Truncate(clust.Timespan,  0.01);
      fData->clust_starttime[i] = Truncate(clust.StartTime, 0.1);
      fData->clust_endtime[i]   = Truncate(clust.EndTime,   0.1);
      fData->clust_ismatch[i]   = clust.isMatched;
      fData->clust_blipid[i]    = clust.BlipID;
      fData->clust_touchtrk[i]  = (clust.TouchTrkID >= 0 );
      fData->clust_touchtrkid[i]= clust.TouchTrkID;
    }
    // if this clust has an associated "trueblip" ID, find it
    // and figure out the true G4 charge, energy, etc
    int tbi = clust.EdepID;
    if( tbi >= 0 && tbi < (int)trueblips.size() ) {
      auto const& trueBlip = trueblips[tbi];
      int g4index = trueBlip.LeadG4Index;
      if( clust.Plane==2 ){
        fData->part_madeClustCol[g4index]  = true;
        fData->edep_madeClustCol[tbi]      = true;
      }

      fData->clust_edepid[i]   = trueBlip.ID;
      // fill histograms of electron/alpha charge resolution,
      // also derive the electron-to-ADC factor (this ~should~ 
      // match up with CalAreaConstants!)
      float q_reco  = clust.Charge;
      float q_anode = trueBlip.NumElectrons;
      float q_dep   = trueBlip.DepElectrons;
      float tdrift  = trueBlip.DriftTime;
      int   pdg     = trueBlip.LeadG4PDG;

      // fill diagnostic histograms for energy deposits from electrons
      if( clust.Plane==fCaloPlane && abs(pdg) == 11 && q_dep > 2000 ) { //&& tdrift > 100 ) {
        h_clust_qres_anode   ->Fill( (q_reco-q_anode)/q_anode );
        h_clust_qres_dep     ->Fill( (q_reco-q_dep)/q_dep );
        h_clust_qres_vs_q    ->Fill( q_dep/1e3, (q_reco-q_dep)/q_dep );
        h_qratio_vs_time_sim ->Fill( tdrift, q_anode/q_dep );
      }
    
    }
      
    if( fDebugMode ) {
      PrintClusterInfo(clust);
      if( clust.TouchTrkID >= 0 ) {
        std::cout<<"!!! This cluster touches track "<<clust.TouchTrkID<<" (blipID "<<clust.BlipID<<")\n";
      }
    }
    
  }//endloop over 2D hit clusters
 
  for(size_t ip=0; ip<kNplanes; ip++){
    h_nclusts[ip]   ->Fill(num_clusts[ip]);
    h_nclusts_pm[ip]->Fill(num_clusts_pm[ip]);
  }//endloop over planes

  //====================================
  // Save blip info to tree
  //===================================
  fData->nblips             = fBlipAlg.blips.size();
  int nblips_matched        = 0;
  int nblips_total          = 0;
  int nblips_picky          = 0;
  float true_blip_charge    = 0;
  for(size_t i=0; i<fBlipAlg.blips.size(); i++){
    if( i > kMaxBlips ) break;
    auto& blp = fBlipAlg.blips[i];
    nblips_total++;
    fNum3DBlips++;
    if( blp.NPlanes >= 3 ) fNum3DBlips3Plane++;
    fData->blip_id[i]         = i;
    fData->blip_tpc[i]        = blp.TPC;
    fData->blip_nplanes[i]    = blp.NPlanes;
    fData->blip_x[i]          = blp.Position.X();
    fData->blip_y[i]          = blp.Position.Y();
    fData->blip_z[i]          = blp.Position.Z();
    fData->blip_sigmayz[i]    = blp.SigmaYZ;
    fData->blip_dx[i]         = blp.dX;
    fData->blip_dw[i]        = blp.dYZ;
    fData->blip_size[i]       = sqrt( pow(blp.dX,2) + pow(blp.dYZ,2) );
    fData->blip_proxtrkdist[i]= blp.ProxTrkDist;
    fData->blip_proxtrkid[i]  = blp.ProxTrkID;
    fData->blip_touchtrk[i]   = (blp.TouchTrkID >= 0 );
    fData->blip_touchtrkid[i] = blp.TouchTrkID;
    fData->blip_incylinder[i] = blp.inCylinder;
    fData->blip_charge[i]     = blp.Charge;
    fData->blip_energy[i]     = blp.Energy;
    fData->blip_yzcorr[i]     = tpcCalib.YZdqdxCorrection(fCaloPlane,blp.Position.Y(),blp.Position.Z());
    
    // Fill cluster charge 2D histograms
    h_blip_charge   ->Fill(blp.Charge);
    h_blip_charge_YU->Fill( 0.001*blp.clusters[2].Charge, 0.001*blp.clusters[0].Charge );
    h_blip_charge_YV->Fill( 0.001*blp.clusters[2].Charge, 0.001*blp.clusters[1].Charge );
    h_blip_charge_UV->Fill( 0.001*blp.clusters[0].Charge, 0.001*blp.clusters[1].Charge );
    for(size_t ipl = 0; ipl<kNplanes; ipl++){
      if( blp.clusters[ipl].NHits <= 0 ) continue;
      fData->blip_clustid[ipl][i] = blp.clusters[ipl].ID;
    }

    // Select picky (high-quality) blips:
    if(blp.NPlanes == 3 && blp.SigmaYZ < 1.) {
      nblips_picky++;
      fNum3DBlipsPicky++;
      h_blip_charge_picky ->Fill(blp.clusters[fCaloPlane].Charge);
       h_blip_zy_picky     ->Fill(blp.Position.Z(), blp.Position.Y());
      h_blip_charge_YU_picky->Fill( 0.001*blp.clusters[2].Charge, 0.001*blp.clusters[0].Charge );
      h_blip_charge_YV_picky->Fill( 0.001*blp.clusters[2].Charge, 0.001*blp.clusters[1].Charge );
      h_blip_charge_UV_picky->Fill( 0.001*blp.clusters[0].Charge, 0.001*blp.clusters[1].Charge );
    }

    h_blip_zy     ->Fill(blp.Position.Z(), blp.Position.Y());
    h_blip_nplanes->Fill(blp.NPlanes);
   
    // -----------------------------------------------
    // save the clustIDs and true energy deposits to the blip
    // (use the association between clust <--> edep)
    // -----------------------------------------------
    if( blp.truth.ID >= 0 && blp.truth.Energy > 0 ) {
      fData->blip_g4id[i]             = blp.truth.LeadG4ID;
      fData->blip_edepid[i]           = blp.truth.ID;
      fData->blip_energyTrue[i]       = blp.truth.Energy;
      fData->edep_blipid[blp.truth.ID]  = blp.ID;
      fNum3DBlipsTrue++;
      nblips_matched++;
      true_blip_charge += blp.truth.NumElectrons;
      h_blip_reszy->Fill( blp.Position.Z()-blp.truth.Position.Z(), blp.Position.Y()-blp.truth.Position.Y() );
      h_blip_resx->Fill( blp.Position.X()-blp.truth.Position.X() );
      h_blip_resE->Fill( blp.truth.Energy, (blp.Energy - blp.truth.Energy) / blp.truth.Energy );
    }
 

  }//endloop over 3D blips
 
  // Fill some more histograms...
  h_nblips->Fill(nblips_total);
  h_nblips_picky->Fill(nblips_picky);
  if( fIsMC ) {
    h_nblips_tm->Fill(nblips_matched);
    if( total_numElectrons ) h_blip_qcomp ->Fill(true_blip_charge / total_numElectrons );
  }
  
  if( fDebugMode ) {
    std::cout<<"\nLooping over "<<fBlipAlg.blips.size()<<" 3D blips:\n";
    for(auto const& b : fBlipAlg.blips ) PrintBlipInfo(b);
  }
 
  
  //====================================
  // Fill TTree
  //====================================
  fData->evtTree->Fill();
  
}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipAna::endJob(){
 
  /*
  fBlipAlg.h_recoWireEff_num->Divide(fBlipAlg.h_recoWireEff_denom);
  fBlipAlg.h_recoWireEff_num->SetOption("hist");
  fBlipAlg.h_recoWireEff_num->SetBit(TH1::kIsAverage);
  fBlipAlg.h_recoWireEffQ_num->Divide(fBlipAlg.h_recoWireEffQ_denom);
  fBlipAlg.h_recoWireEffQ_num->SetOption("hist");
  fBlipAlg.h_recoWireEffQ_num->SetBit(TH1::kIsAverage);
  */

  float nEvents   =  float(fNumEvents);
  float scalefac  = 1/nEvents;

  h_trk_length->Scale(scalefac);
  h_trk_xspan->Scale(scalefac);
  h_trks_100cm[0]->Scale(scalefac);
  h_trks_100cm[1]->Scale(scalefac);
  
  h_trk_length->SetOption("hist");
  h_trk_xspan->SetOption("hist");
  h_trks_100cm[0]->SetOption("hist");
  h_trks_100cm[1]->SetOption("hist");

  BlipUtils::NormalizeHist(h_clust_qres_anode);
  BlipUtils::NormalizeHist(h_clust_qres_dep);
  for(size_t i=0; i<kNplanes; i++){
    BlipUtils::NormalizeHist(h_hit_sigmaint[i]);
    BlipUtils::NormalizeHist(h_hit_adcdiff[i]);
    BlipUtils::NormalizeHist(h_hitamp[i]);
    BlipUtils::NormalizeHist(h_hitamp_mip[i]);
    BlipUtils::NormalizeHist(h_hitamp_alpha[i]);
    BlipUtils::NormalizeHist(h_hitgof_mip[i]);
    BlipUtils::NormalizeHist(h_hitgof_iso[i]);
    BlipUtils::NormalizeHist(h_hitgof_isomatch[i]);
    BlipUtils::NormalizeHist(h_hitgof_alpha[i]);
    BlipUtils::NormalizeHist(h_hitmult[i]);
    BlipUtils::NormalizeHist(h_hitmult_mip[i]);
    BlipUtils::NormalizeHist(h_hitmult_iso[i]);
    BlipUtils::NormalizeHist(h_hitmult_isomatch[i]);
    BlipUtils::NormalizeHist(h_hitmult_iso_true[i]);
    BlipUtils::NormalizeHist(h_hitmult_isomatch_true[i]);
    BlipUtils::NormalizeHist(h_hitmult_alpha[i]);
    BlipUtils::NormalizeHist(h_hitrms[i]);
    BlipUtils::NormalizeHist(h_hitrms_mip[i]);
    BlipUtils::NormalizeHist(h_hitrms_iso[i]);
    BlipUtils::NormalizeHist(h_hitrms_isomatch[i]);
    BlipUtils::NormalizeHist(h_hitrms_iso_true[i]);
    BlipUtils::NormalizeHist(h_hitrms_isomatch_true[i]);
    BlipUtils::NormalizeHist(h_hitrms_alpha[i]);
  }

  printf("\n***********************************************\n");
  fBlipAlg.PrintConfig();
  printf("BlipAna Summary\n\n");
  printf("  Total events                : %i\n",        fNumEvents);
  printf("  Blips per evt, total        : %.3f\n",      fNum3DBlips/nEvents);
  printf("                 3 planes     : %.3f\n",      fNum3DBlips3Plane/nEvents);
  printf("                 picky        : %.3f\n",      fNum3DBlipsPicky/nEvents);
  printf("                 picky frac   : %5.3f\n",     fNum3DBlipsPicky/float(fNum3DBlips));
  
  if(fIsMC){
  printf("  MC-matched blips per evt    : %.3f\n",       fNum3DBlipsTrue/nEvents);
  if( h_blip_qcomp->GetMean() > 0 ) 
  printf("  Charge completeness, total  : %.4f +/- %.4f\n", h_blip_qcomp->GetMean(), h_blip_qcomp->GetStdDev()/sqrt(fNumEvents));
  //printf("                       < 2MeV : %.4f +/- %.4f\n", h_blip_qcomp_2MeV->GetMean(), h_blip_qcomp_2MeV->GetStdDev()/sqrt(fNumEvents));
  //printf("  Blip purity                 : %.4f\n",       h_blip_pur->GetMean());
  }
  printf("  Mean blip charge            : %.0f\n",      h_blip_charge->GetMean());
  printf("\n");
  for(size_t i=0; i<kNplanes; i++){
  printf("  Plane %lu -------------------------\n",i);
  printf("   * total hits/evt           : %.2f\n",fNumHits[i]/(float)fNumEvents);
  printf("   * untracked hits/evt       : %.2f (%.2f plane-matched)\n",fNumHitsUntracked[i]/(float)fNumEvents, fNumHitsMatched[i]/(float)fNumEvents);
  //printf("   * plane-matched hits/evt   : %.2f\n",fNumHitsMatched[i]/(float)fNumEvents);
  if(fIsMC) {
  printf("   * true-matched hits/evt    : %.2f (%.2f plane-matched)\n",fNumHitsTrue[i]/(float)fNumEvents, fNumHitsMatchedTrue[i]/(float)fNumEvents);
  if( h_chargecomp[i]->GetMean() > 0 ) 
  printf("   * charge completeness      : %.4f\n",h_chargecomp[i]->GetMean());
  //printf("   * hit purity               : %.4f\n",h_hitpur[i]->GetMean());
  }
  } 
  printf("\n***********************************************\n");

  if( fDoACPTrkCalib ) {
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    float driftVelocity = detProp.DriftVelocity();
    float calcLifetime = 999999.;
    if( h_ACPtrk_qratio->GetMean() < 1. ){ 
      //calcLifetime = (200./detProp->DriftVelocity()) / std::log( 1./h_ACPtrk_qratio->GetMean() );
      calcLifetime = (200./driftVelocity) / std::log( 1./h_ACPtrk_qratio->GetMean() );
    }
    printf("\n***********************************************\n");
    printf("Anode-cathode piercing track calibrations\n");
    printf("  - tracks used     : %i\n",(int)h_ACPtrk_qratio->GetEntries());
    printf("  - 2m attenuation  : %f +/- %f\n",h_ACPtrk_qratio->GetMean(),h_ACPtrk_qratio->GetMeanError());
    printf("  - lifetime        : %f ms\n",calcLifetime/1000.);
    printf("************************************************\n");
  }

}

//#################################################
// Neutrino selection functions
//#################################################
void BlipAna::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
{
  _pfpmap.clear();
  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
  { _pfpmap[pfp_pxy->Self()] = p; p++; }
  return;
} // BuildPFPMap

void BlipAna::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                           const ProxyPfpColl_t &pfp_pxy_col,
                           std::vector<ProxyPfpElem_t> &slice_v)
{
  auto daughters = pfp_pxy->Daughters();
  slice_v.push_back(pfp_pxy);
  //std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  for (auto const &daughterid : daughters){
    if (_pfpmap.find(daughterid) == _pfpmap.end())
    {
      //std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }

    // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
      ++pfp_pxy2;
    // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));
    AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
  } // for all daughters

  return;
} // AddDaughters



//###################################################
//  Printouts for debugging
//###################################################

void BlipAna::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-10i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, Npts=%4i, KE0=%8.3f, Edep=%8.3f, T=%10.2f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->part_g4id[i],
   fData->part_pdg[i],
   fData->part_startPointx[i],
   fData->part_startPointy[i],
   fData->part_startPointz[i],
   fData->part_pathlen[i],
   fData->part_numTrajPts[i],
   fData->part_KE[i],
   fData->part_depEnergy[i],
   fData->part_startT[i]/1e3,
   fData->part_mother[i],
   fData->part_process[i].c_str(),
   fData->part_nDaughters[i]
  ); 
}

void BlipAna::PrintTrueBlipInfo(const blipobj::TrueBlip& tb){
  printf("  edepID: %5i  G4ID: %-6i PDG: %-10i XYZ: %7.2f, %7.2f, %7.2f, %8.3f MeV, %8i e- deposited, %8i e- @anode,  %12s\n",
   tb.ID,
   tb.LeadG4ID,
   tb.LeadG4PDG,
   tb.Position.X(),
   tb.Position.Y(),
   tb.Position.Z(),
   tb.Energy,
   tb.DepElectrons,
   tb.NumElectrons,
   fData->part_process[tb.LeadG4Index].c_str()
  ); 
}

void BlipAna::PrintHitInfo(const blipobj::HitInfo& hi){
  printf("  hitID: %4i, TPC: %i, plane: %i, driftTicks: %7.2f, leadWire: %3i, G4ID: %4i, recoTrack: %4i\n",
    hi.hitid,
    hi.tpc,
    hi.plane,
    hi.driftTime,
    hi.wire,
    hi.g4trkid,
    hi.trkid
  );
}

void BlipAna::PrintClusterInfo(const blipobj::HitClust& hc){
  printf("  clustID: %4i, TPC: %i, plane: %i, time range: %7.2f - %7.2f, timespan: %6.2f, leadWire: %3i, nwires: %3i, nhits: %3i, edepid: %i, isMatched: %i, blipID: %i\n",
    hc.ID,
    hc.TPC,
    hc.Plane,
    hc.StartTime,
    hc.EndTime,
    hc.Timespan,
    hc.CenterWire,
    hc.NWires,
    hc.NHits,
    hc.EdepID,
    hc.isMatched,
    hc.BlipID
  );
}

void BlipAna::PrintBlipInfo(const blipobj::Blip& bl){
  printf("  blipID: %4i, TPC: %i, charge: %8.0i,  recoEnergy: %8.3f MeV, XYZ: %6.1f, %6.1f, %6.1f,   EdepID: %i\n",
  bl.ID,
  bl.TPC,
  (int)bl.clusters[2].Charge,
  bl.Energy,
  bl.Position.X(),bl.Position.Y(),bl.Position.Z(),
  bl.truth.ID
  );
}

float BlipAna::Truncate(float input, double base){
  if( base > 0 )  return roundf( input / base ) * base;
  else            return input;
  return input;
}


DEFINE_ART_MODULE(BlipAna)

#endif
