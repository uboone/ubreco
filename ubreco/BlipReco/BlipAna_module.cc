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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "cetlib/search_path.h"

// MicroBooNE-specific includes
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

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
#include "TF1.h"

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
const int kMaxHits  = 25000;
const int kMaxTrks  = 10000;
const int kMaxShwrs = 1000;
const int kMaxBlips = 5000;
const int kMaxG4    = 10000;
const int kMaxEDeps = 10000;
   
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
  bool  saveTrkInfo         = true;
  bool  saveHitInfo         = true;
  bool  saveClustInfo       = true;

  // --- Event information ---   
  int           event;                    // event number
  int           run;                      // run number
  unsigned int  timestamp;                // unix time of event
  float         lifetime;                 // electron lifetime

  // --- G4 information ---
  int   nparticles;               // number of G4 particles
  bool  part_isPrimary[kMaxG4];        // is primary particle
  int   part_trackID[kMaxG4];          // G4 track ID
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
  float part_depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  int   part_depElectrons[kMaxG4];     // electrons deposited
  //float part_numElectrons[kMaxG4];     // electrons reaching anode wires
  std::vector<std::string> part_process;// process name
  //float total_depEnergy;          // total deposited energy in AV
  //int   total_depElectrons;       // total deposited ionization electrons in AV
  //float total_numElectrons;       // total electrons reaching anode wires

  // --- True energy deposit info (derived from SimChannels and SimEnergyDeposits) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 track ID
  int   edep_g4index[kMaxEDeps];  // leading G4 track index
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_clustid[kMaxEDeps];  // hitclust ID
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
  int   edep_electrons[kMaxEDeps];// total ionization electrons deposited
  int   edep_charge[kMaxEDeps];   // total electrons reaching anode wires
  float edep_tdrift[kMaxEDeps];   // drift time for this energy dep
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)

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
  float hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
  int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
  float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
  float hit_g4energy[kMaxHits];   // true energy
  float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
  int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
  int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
  float hit_gof[kMaxHits];        // goodness of fit (default -1)

  // --- Track information ---
  int   ntrks;                    // number tracks
  int   trk_id[kMaxTrks];         // trackID
  int   trk_npts[kMaxTrks];       // number 3D trajectory points
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
  int   clust_id[kMaxHits];           // cluster ID (index)
  int   clust_tpc[kMaxHits];          // cluster TPC ID
  int   clust_plane[kMaxHits];        // cluster plane
  int   clust_wire[kMaxHits];         // central-most wire of cluster
  int   clust_startwire[kMaxHits];    // starting wire
  int   clust_endwire[kMaxHits];      // ending wire
  int   clust_nwires[kMaxHits];       // number of wires in this cluster
  int   clust_nhits[kMaxHits];        // number of hits
  float clust_time[kMaxHits];         // charge-weighted time
  float clust_timespan[kMaxHits];     // cluster timespan
  float clust_rms[kMaxHits];          // charge-weighted RMS
  float clust_starttime[kMaxHits];    // cluster start tick
  float clust_endtime[kMaxHits];      // cluster end tick
  float clust_amp[kMaxHits];          // maximum hit amplitude [ADC]
  float clust_charge[kMaxHits];       // cluster charge at anode [e-]
  float clust_g4charge[kMaxHits];     // true cluster charge at anode
  float clust_g4energy[kMaxHits];     // true cluster energy from G4
  int   clust_g4id[kMaxHits];         // true MCParticle ID (index for particle branches)
  int   clust_blipid[kMaxHits];       // blip ID for this nlusteer (if it was made into one)
  int   clust_edepid[kMaxHits];       // true energy dep ID
  bool  clust_ismatch[kMaxHits];      // was this cluster plane-matched?
  //int   clust_lhit_wire[kMaxHits];    // cluster wire (lead hit wire)
  //int   clust_lhit_chan[kMaxHits];    // cluster channel (lead hit wire)
  //int   clust_lhit_id[kMaxHits];      // lead hit ID (index for hit_X[i] branches)
  //float clust_lhit_amp[kMaxHits];     // lead hit peak amplitude [ADC]
  //float clust_lhit_rms[kMaxHits];     // lead hit RMS [ADC]
  //float clust_lhit_time[kMaxHits];    // lead hit time [ticks]
  //float clust_lhit_gof[kMaxHits];     // lead hit goodness-of-fit; pulse train = -1
  //bool  clust_lhit_isfit[kMaxHits];   // is there a valid goodness of fit for lead hit?

  // --- 3D Blip information ---
  int   nblips;                       // number of blips in event
  int   blip_id[kMaxBlips];           // blip ID / index
  int   blip_tpc[kMaxBlips];          // blip TPC
  int   blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_sigmayz[kMaxBlips];      // difference in wire intersection points
  float blip_dx[kMaxBlips];           // dX [cm]
  float blip_dyz[kMaxBlips];
  int   blip_trkid[kMaxBlips];
  float blip_trklength[kMaxBlips];       // if this blip was tracked, get 3D length
  float blip_size[kMaxBlips];       // rough size estimation based on time-tick extent and wire span
  float blip_sumadc[kMaxBlips];          // integrated ADCs 
  float blip_charge[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip energy [MeV]
  int   blip_edepid[kMaxBlips];       // true energy dep ID
  float blip_proxtrkdist[kMaxBlips];  // distance to nearest track
  int   blip_proxtrkid[kMaxBlips];    // index of nearest trk
  bool  blip_incylinder[kMaxBlips];   // is blip within a cylinder near a track
  float blip_time[kMaxBlips];                   // blip drift time [ticks]
  int   blip_clustid[kNplanes][kMaxBlips];     // cluster ID per plane
  /*
  int   blip_clust_wire[kNplanes][kMaxBlips];
  int   blip_clust_startwire[kNplanes][kMaxBlips];
  int   blip_clust_endwire[kNplanes][kMaxBlips];
  int   blip_clust_nhits[kNplanes][kMaxBlips];  // hits per cluster on each plane
  int   blip_clust_nwires[kNplanes][kMaxBlips]; // wires per cluster on each plane
  float blip_clust_charge[kNplanes][kMaxBlips];
  float blip_clust_lhit_rms[kNplanes][kMaxBlips];
  float blip_clust_lhit_amp[kNplanes][kMaxBlips];
  */

  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    lifetime              = -999;
    timestamp             = -999;
    //total_depEnergy       = -999;
    //total_depElectrons    = -999;
    //total_numElectrons    = -999;
    nparticles            = 0;    // --- G4 particles ---
    FillWith(part_isPrimary,   false);
    FillWith(part_trackID,     -999);
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
    FillWith(part_depElectrons,-999);
    //FillWith(part_numElectrons,-999);
    FillWith(part_depEnergy,   -999.);
    FillWith(part_process,     "");
    nedeps                = 0;    // --- EDeps ---
    FillWith(edep_tpc,    -9);
    FillWith(edep_energy, -999);
    FillWith(edep_electrons,  -999);
    FillWith(edep_charge, -999);
    FillWith(edep_tdrift, -999);
    FillWith(edep_x,      -99999.);
    FillWith(edep_y,      -99999.);
    FillWith(edep_z,      -99999.);
    FillWith(edep_g4id,   -9);
    FillWith(edep_g4index,   -9);
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
    FillWith(hit_ismatch, -9);
    FillWith(hit_trkid,   -9);
    FillWith(hit_g4id,    -999);
    FillWith(hit_g4frac,  -9);
    FillWith(hit_g4energy,-999);
    FillWith(hit_g4charge,-999);
    FillWith(hit_clustid, -9);
    FillWith(hit_blipid,  -9);
    FillWith(hit_gof,     -9);
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
    nclusts                   = 0;    // --- Hit Clusters ---
    FillWith(clust_id,        -9);
    FillWith(clust_tpc,       -9);
    FillWith(clust_plane,     -9);
    FillWith(clust_nwires,    -9);
    FillWith(clust_nhits,     -9);
    FillWith(clust_wire,      -9);
    FillWith(clust_startwire, -9);
    FillWith(clust_endwire,   -9);
    FillWith(clust_time,      -999);
    FillWith(clust_timespan,  -9);
    FillWith(clust_starttime, -999);
    FillWith(clust_endtime,   -999);
    FillWith(clust_rms,       -9);
    FillWith(clust_amp,       -9);
    FillWith(clust_charge,    -999);
    FillWith(clust_g4id,      -9);
    FillWith(clust_g4charge,  -9999);
    FillWith(clust_g4energy,  -9);
    FillWith(clust_edepid,    -9);
    FillWith(clust_blipid,    -9);
    FillWith(clust_ismatch,   false);
    nblips                    = 0;
    FillWith(blip_id,         -9);
    FillWith(blip_tpc,        -9);
    FillWith(blip_nplanes,    -9);
    FillWith(blip_x,          -9999);
    FillWith(blip_y,          -9999);
    FillWith(blip_z,          -9999);
    FillWith(blip_sigmayz,    -9);
    FillWith(blip_dx,         -9);
    FillWith(blip_dyz,        -9);
    FillWith(blip_trklength,     -9);
    FillWith(blip_trkid,      -9);
    FillWith(blip_size,     -9);
    FillWith(blip_sumadc,     -999);
    FillWith(blip_charge,     -999);
    FillWith(blip_energy,     -999);
    FillWith(blip_proxtrkdist,  -99);
    FillWith(blip_proxtrkid,   -9);
    FillWith(blip_incylinder, false);
    FillWith(blip_edepid,     -9);
    for(int i=0; i<kNplanes; i++) {
      FillWith(blip_clustid[i],        -9);
      /*
      FillWith(blip_clust_wire[i],      -9);
      FillWith(blip_clust_startwire[i], -9);
      FillWith(blip_clust_endwire[i],   -9);
      FillWith(blip_clust_nhits[i],     -9);
      FillWith(blip_clust_nwires[i],    -9);
      FillWith(blip_clust_charge[i],    -999);
      FillWith(blip_clust_lhit_rms[i],  -9);
      FillWith(blip_clust_lhit_amp[i],  -9);
      */
    }
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
    evtTree->Branch("timestamp",&timestamp,"timestamp/i");
    evtTree->Branch("lifetime",&lifetime,"lifetime/F");
  
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
      evtTree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
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
      evtTree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
      evtTree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
      evtTree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
      evtTree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
      evtTree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
      evtTree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
    }

    if( saveClustInfo ) {
      evtTree->Branch("nclusts",        &nclusts,       "nclusts/I");
      evtTree->Branch("clust_plane",    clust_plane,    "clust_plane[nclusts]/I");
      evtTree->Branch("clust_wire",     clust_wire,     "clust_wire[nclusts]/I");
      evtTree->Branch("clust_startwire",clust_startwire,"clust_startwire[nclusts]/I");
      evtTree->Branch("clust_endwire",  clust_endwire,  "clust_endwire[nclusts]/I");
      evtTree->Branch("clust_nwires",   clust_nwires,   "clust_nwires[nclusts]/I");
      evtTree->Branch("clust_nhits",    clust_nhits,    "clust_nhits[nclusts]/I");
      evtTree->Branch("clust_time",     clust_time,     "clust_time[nclusts]/F");
      evtTree->Branch("clust_rms",      clust_rms,      "clust_rms[nclusts]/F");
      evtTree->Branch("clust_amp",      clust_amp,      "clust_amp[nclusts]/F");
      evtTree->Branch("clust_charge",   clust_charge,   "clust_charge[nclusts]/F");
      evtTree->Branch("clust_ismatch",  clust_ismatch,  "clust_ismatch[nclusts]/O");
      evtTree->Branch("clust_blipid",   clust_blipid,   "clust_blipid[nclusts]/I");
      if( saveTruthInfo ) 
      evtTree->Branch("clust_edepid",   clust_edepid,   "clust_edepid[nclusts]/I");
    }

    evtTree->Branch("nblips",&nblips,"nblips/I");
    evtTree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/I");
    evtTree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
    evtTree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
    evtTree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
    evtTree->Branch("blip_sigmayz",blip_sigmayz,"blip_sigmayz[nblips]/F");
    evtTree->Branch("blip_dx",blip_dx,"blip_dx[nblips]/F");
    evtTree->Branch("blip_dyz",blip_dyz,"blip_dyz[nblips]/F");
    evtTree->Branch("blip_size",blip_size,"blip_size[nblips]/F");
    evtTree->Branch("blip_trkid",blip_trkid,"blip_trkid[nblips]/I");
    evtTree->Branch("blip_trklength",blip_trklength,"blip_trklength[nblips]/F");
    evtTree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
    evtTree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
    evtTree->Branch("blip_incylinder",blip_incylinder,"blip_incylinder[nblips]/O");
    evtTree->Branch("blip_proxtrkdist",blip_proxtrkdist,"blip_proxtrkdist[nblips]/F");
    if( saveTrkInfo ) evtTree->Branch("blip_proxtrkid",blip_proxtrkid,"blip_proxtrkid[nblips]/I");
    if( saveTruthInfo ) evtTree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
   
    for(int i=0;i<kNplanes;i++) {
      if( saveClustInfo ) evtTree->Branch(Form("blip_pl%i_clustid",i),blip_clustid[i],Form("blip_pl%i_clustid[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_wire",i),     blip_clust_wire[i],Form("blip_pl%i_wire[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_startwire",i),blip_clust_startwire[i],Form("blip_pl%i_startwire[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_endwire",i),  blip_clust_endwire[i],Form("blip_pl%i_endwire[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_nhits",i),    blip_clust_nhits[i],Form("blip_pl%i_nhits[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_nwires",i),   blip_clust_nwires[i],Form("blip_pl%i_nwires[nblips]/I",i));
      //evtTree->Branch(Form("blip_pl%i_charge",i),   blip_clust_charge[i],Form("blip_pl%i_charge[nblips]/F",i));
      //evtTree->Branch(Form("blip_pl%i_lhit_rms",i),   blip_clust_lhit_rms[i],Form("blip_pl%i_lhit_rms[nblips]/F",i));
      //evtTree->Branch(Form("blip_pl%i_lhit_amp",i),   blip_clust_lhit_amp[i],Form("blip_pl%i_lhit_amp[nblips]/F",i));
    }
    
    if( saveTruthInfo ) {
      evtTree->Branch("nparticles",&nparticles,"nparticles/I");
      evtTree->Branch("part_isPrimary",part_isPrimary,"part_isPrimary[nparticles]/O");
      evtTree->Branch("part_trackID",part_trackID,"part_trackID[nparticles]/I");
      evtTree->Branch("part_pdg",part_pdg,"part_pdg[nparticles]/I");
      evtTree->Branch("part_nDaughters",part_nDaughters,"part_nDaughters[nparticles]/I");
      evtTree->Branch("part_mother",part_mother,"part_mother[nparticles]/I");
      evtTree->Branch("part_KE",part_KE,"part_KE[nparticles]/F");
      //evtTree->Branch("part_endKE",part_endKE,"part_endKE[nparticles]/F");
      //evtTree->Branch("part_mass",part_mass,"part_mass[nparticles]/F");
      //evtTree->Branch("part_P",part_P,"part_P[nparticles]/F");
      evtTree->Branch("part_Px",part_Px,"part_Px[nparticles]/F");
      evtTree->Branch("part_Py",part_Py,"part_Py[nparticles]/F");
      evtTree->Branch("part_Pz",part_Pz,"part_Pz[nparticles]/F");
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
      //evtTree->Branch("part_process",&part_process);
      
      evtTree->Branch("nedeps",&nedeps,"nedeps/I");
      evtTree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
      //evtTree->Branch("edep_g4index",edep_g4index,"edep_g4index[nedeps]/I"); 
      evtTree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
      evtTree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
      evtTree->Branch("edep_clustid",edep_clustid,"edep_clustid[nedeps]/I"); 
      evtTree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
      evtTree->Branch("edep_electrons",edep_electrons,"edep_electrons[nedeps]/I"); 
      evtTree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/I"); 
      evtTree->Branch("edep_tdrift",edep_tdrift,"edep_tdrift[nedeps]/F"); 
      evtTree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
      evtTree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
      evtTree->Branch("edep_z",edep_z,"edep_z[nedeps]/F"); 
    }
    
    
  }//end MakeTree

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
  void    PrintTrueBlipInfo(const blip::TrueBlip&);
  void    PrintClusterInfo(const blip::HitClust&);
  void    PrintHitInfo(const blip::HitInfo&);
  //float   Truncate(float, double = 0.1);
  //double  Truncate(double, double = 0.1);

  // --- Data and calo objects ---
  BlipAnaTreeDataStruct*  fData;
  blip::BlipRecoAlg*      fBlipAlg;

  //TGraph2D*               ESTAR;

  // --- FCL configs ---
  bool                fDebugMode;
  std::string         fHitProducer;
  std::string         fTrkProducer;
  std::string         fGeantProducer;
  std::string         fSimDepProducer;
  int                 fCaloPlane;
  std::vector<bool>   fSavePlaneInfo;

  // --- Counters and such ---
  bool  fIsRealData         = false;
  bool  fIsMC               = false;
  int   fNumEvents          = 0;
  int   fNumHits[3]         = {};
  int   fNumHitsUntracked[3]= {};
  int   fNumHitsMatched[3]  = {};
  int   fNumHitsTrue[3]     = {};
  int   fNum3DBlips         = 0;
  int   fNum3DBlips3Plane   = 0;
  int   fNum3DBlipsPicky    = 0;
  int   fNum3DBlipsTrue     = 0;

  // --- Histograms ---
  TH1D*   h_nhits[kNplanes];
  TH1D*   h_nhits_ut[kNplanes];
  TH1D*   h_nhits_m[kNplanes];
  TH1D*   h_nhits_tm[kNplanes];
  TH1D*   h_hitamp[kNplanes];
  //TH1D*   h_hitsigt[kNplanes];
  TH1D*   h_hitrms[kNplanes];
  TH1D*   h_hitratio[kNplanes];
  //TH1D*   h_hitint[kNplanes];
  TH1D*   h_hitmult[kNplanes];
  TH2D*   h_nelec_TrueVsReco[kNplanes];
  TH1D*   h_nelec_Resolution[kNplanes];
  TH1D*   h_chargecomp[kNplanes];
  TH1D*   h_hitpur[kNplanes];
  TH1D*   h_clust_nwires;
  TH1D*   h_clust_timespan;
  
  TH1D*   h_ntrks;
  TH1D*   h_trk_length;
  TH1D*   h_trk_xspan;
  
  TH1D*   h_nblips;
  TH1D*   h_nblips_picky;
  TH1D*   h_nblips_tm;
  TH1D*   h_blip_nplanes;
  TH1D*   h_blip_qcomp;
  TH1D*   h_blip_qcomp_2MeV;
  TH1D*   h_blip_pur;
  TH2D*   h_blip_reszy;
  TH1D*   h_blip_resx;
  TH2D*   h_blip_zy;
  TH2D*   h_blip_zy_picky;
  TH1D*   h_blip_sumadc;
  TH1D*   h_blip_charge;
  TH1D*   h_blip_charge_picky;
 
  TH2D*   h_blip_charge_YU;
  TH2D*   h_blip_charge_YV;
  TH2D*   h_blip_charge_UV;
  
  TH2D*   h_blip_charge_YU_picky;
  TH2D*   h_blip_charge_YV_picky;
  TH2D*   h_blip_charge_UV_picky;

  // Some truth metrics for debugging
  TH1D*   h_true_lifetime;
  TH1D*   h_qres_electrons;
  TH1D*   h_qres_alphas;
  TH1D*   h_adc_factor;


  // Initialize histograms
  void InitializeHistograms(){

    art::ServiceHandle<art::TFileService> tfs;
  
    float blipMax = 500;    int blipBins  = 500;
    h_nblips          = tfs->make<TH1D>("nblips","Reconstructed 3D blips per event",blipBins,0,blipMax);
    h_blip_zy         = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy         ->SetOption("COLZ");
    h_blip_zy_picky         = tfs->make<TH2D>("blip_zy_picky","3D blip location (3-plane match, intersect #Delta < 1 cm);Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy_picky         ->SetOption("COLZ");
      
    art::TFileDirectory dir_diag = tfs->mkdir("Diagnostics");
    h_nblips_picky    = dir_diag.make<TH1D>("nblips_picky","Reconstructed 3D blips per event (3-plane match, intersect #Delta < 1 cm)",blipBins,0,blipMax);
    h_blip_nplanes    = dir_diag.make<TH1D>("blip_nplanes","Matched planes per blip",3,1,4);
    h_blip_sumadc     = dir_diag.make<TH1D>("blip_sumadc","3D blips;Integral on coll plane [ADC]",            200,0,500);
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

    float hitMax  = 15000; int hitBins  = 1500;
    float ampMax  = 50;   int ampBins   = 250;
    float rmsMax  = 30;   int rmsBins   = 300;
    float ratioMax = 5;   int ratioBins = 250;

    // MC histograms related to truth
    art::TFileDirectory dir_truth = dir_diag.mkdir("Truth");
    h_nblips_tm     = dir_truth.make<TH1D>("nblips_tm","Truth-matched 3D blips per event",blipBins,0,blipMax);
    h_blip_qcomp    = dir_truth.make<TH1D>("blip_qcomp","Fraction of true charge (at anode) reconstructed into 3D blips",202,0,1.01);
    h_blip_qcomp_2MeV = dir_truth.make<TH1D>("blip_qcomp_2MeV","Fraction of true charge (at anode) reconstructed into 3D blips, dE < 2 MeV",202,0,1.01);
    h_blip_pur      = dir_truth.make<TH1D>("blip_pur","Fraction of truth-matched blips",202,0,1.01);
    h_blip_reszy    = dir_truth.make<TH2D>("blip_res_zy","Blip position resolution;Z_{reco} - Z_{true} [cm];Y_{reco} - Y_{true} [cm]",150,-15,15,150,-15,15);
      h_blip_reszy  ->SetOption("colz");
    h_blip_resx     = dir_truth.make<TH1D>("blip_res_x","Blip position resolution;X_{reco} - X_{true} [cm]",150,-15,15);
    h_true_lifetime   = dir_truth.make<TH1D>("blip_true_lifetime","Calculated charge attenuation;Lifetime [#mus]",200,0,200e3);
    h_qres_electrons= dir_truth.make<TH1D>("qres_electrons","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_qres_alphas   = dir_truth.make<TH1D>("qres_alphas","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_adc_factor    = dir_truth.make<TH1D>("adc_per_e","Collection plane;ADC per electron",200,0,0.01);
    
    h_ntrks           = dir_diag.make<TH1D>("ntrks","Number of reconstructed tracks per event",150,0,150);
    h_trk_length      = dir_diag.make<TH1D>("trk_length",";Track length [cm]",1000,0,500);
    h_trk_xspan       = dir_diag.make<TH1D>("trk_xspan",";Track dX [cm]",300,0,300);
    
    for(int i=0; i<kNplanes; i++) {
      h_nhits[i]      = dir_diag.make<TH1D>(Form("pl%i_nhits",i),  Form("Plane %i;total number of hits",i),hitBins,0,hitMax);
      h_nhits_ut[i]   = dir_diag.make<TH1D>(Form("pl%i_nhits_untracked",i),  Form("Plane %i;total number of untracked hits",i),hitBins,0,hitMax);
      h_nhits_m[i]    = dir_diag.make<TH1D>(Form("pl%i_nhits_planematched",i), Form("Plane %i;total number of untracked plane-matched hits",i),hitBins,0,hitMax);
      h_hitamp[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_amp",i), Form("Plane %i untracked hits;hit amplitude [ADC]",i),ampBins,0,ampMax);
      h_hitrms[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_rms",i), Form("Plane %i untracked hits;hit RMS [ADC]",i),rmsBins,0,rmsMax);
      h_hitratio[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_ratio",i), Form("Plane %i untracked hits;hit RMS/Amplitude ratio",i),ratioBins,0,ratioMax);
      //h_hitint[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_integral",i), Form("Plane %i hits;hit integral [ADC]",i),areaBins,0,areaMax);
      h_hitmult[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_multiplicity",i), Form("Plane %i untracked hits;multiplicity [ADC]",i),30,0,30);
      //h_hitsigt[i]    = dir_diag.make<TH1D>(Form("pl%i_hit_sigt",i), Form("Plane %i untracked hits;hit time uncertainty [ADC]",i),200,0,5);
      h_nhits_tm[i]   = dir_truth.make<TH1D>(Form("pl%i_nhits_truthmatched",i), Form("Plane %i;number of untracked truth-matched hits",i),hitBins,0,hitMax);
      h_chargecomp[i] = dir_truth.make<TH1D>(Form("pl%i_hit_charge_completeness",i),Form("charge completness, plane %i",i),101,0,1.01);
      h_hitpur[i]     = dir_truth.make<TH1D>(Form("pl%i_hit_purity",i),Form("hit purity, plane %i",i),101,0,1.01);
      h_nelec_TrueVsReco[i] = dir_truth.make<TH2D>( Form("pl%i_nelec_TrueVsReco",i),
        Form("Plane %i;true hit charge [ #times 10^{3} electrons ];Reconstructed hit charge [ #times 10^{3} electrons ]",i),100,0,50, 100,0,50);
        h_nelec_TrueVsReco[i] ->SetOption("colz");
      h_nelec_Resolution[i] = dir_truth.make<TH1D>( Form("pl%i_nelec_res",i),Form("Plane %i;hit charge resolution: (reco-true)/true",i),300,-3,3);
    }//endloop over planes
    
 }

};//class BlipAna


//###################################################
//  BlipAna constructor and destructor
//###################################################
BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData  (nullptr)
{
  // blip reconstruction algorithm class
  fhicl::ParameterSet pset_blipalg = pset.get<fhicl::ParameterSet>("BlipAlg");
  fBlipAlg        = new blip::BlipRecoAlg( pset_blipalg );
  fHitProducer    = pset_blipalg.get<std::string>   ("HitProducer",     "gaushit");
  fTrkProducer    = pset_blipalg.get<std::string>   ("TrkProducer",     "pandora");
  fGeantProducer  = pset_blipalg.get<std::string>   ("GeantProducer",   "largeant");
  fSimDepProducer = pset_blipalg.get<std::string>   ("SimEDepProducer", "ionization");
  fCaloPlane      = pset_blipalg.get<int>           ("CaloPlane",       2);
  fDebugMode      = pset.get<bool>                  ("DebugMode",       false);
  fSavePlaneInfo  = pset.get<std::vector<bool>>     ("SavePlaneInfo",   {true,true,true});

  // data tree object
  fData = new BlipAnaTreeDataStruct();
  fData ->treeName        = pset.get<std::string> ("EventTreeName", "anatree");
  fData ->saveTruthInfo   = pset.get<bool>        ("SaveTruthInfo", true);
  fData ->saveTrkInfo     = pset.get<bool>        ("SaveTrkInfo",   true);
  fData ->saveHitInfo     = pset.get<bool>        ("SaveHitInfo",   true);
  fData ->saveClustInfo   = pset.get<bool>        ("SaveClustInfo", true);
  fData ->Clear();
  fData ->MakeTree();

  // initialize histograms
  InitializeHistograms();
}
BlipAna::~BlipAna(){}



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
  fIsRealData       = evt.isRealData();
  fNumEvents++;

  // Get timestamp
  unsigned long long int tsval = evt.time().value();
  const unsigned long int mask32 = 0xFFFFFFFFUL;
  fData->timestamp = ( tsval >> 32 ) & mask32;
  
  // Retrieve lifetime
  const lariov::UBElectronLifetimeProvider& elifetime_provider = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
  float electronLifetime = elifetime_provider.Lifetime() * /*convert ms->mus*/ 1e3;
  fData->lifetime = electronLifetime;

  //============================================
  // Run blip reconstruction: 
  //============================================
  
  fBlipAlg->RunBlipReco(evt);
  
  //  
  //  In the above step, we pass the entire art::Event to the algorithm, 
  //  and it creates a single collection of blip 'objects', a special data
  //  struct in the 'blip' namespace defined in BlipUtils.h.
  //  
  //  We can then retrieve these blips and incorporate them into
  //  our analysis however we like:
  //
  //    std::vector<blip::Blip> = fBlipAlg->blips;
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


  // Tell us what's going on!
  std::cout<<"\n"
  <<"=========== BlipAna =========================\n"
  <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"; total: "<<fNumEvents<<"\n";
  std::cout<<"Lifetime is "<<electronLifetime<<" microseconds\n";
  
  

  //=======================================
  // Get data products for this event
  //========================================
  
  // -- G4 particles
  art::Handle< std::vector<simb::MCParticle> > pHandle;
  std::vector<art::Ptr<simb::MCParticle> > plist;
  if (evt.getByLabel("largeant",pHandle))
    art::fill_ptr_vector(plist, pHandle);
  
  // -- hits (from input module)
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitProducer,hitHandle))
    art::fill_ptr_vector(hitlist, hitHandle);
  
  // -- hits (from gaushit)
  art::Handle< std::vector<recob::Hit> > hitHandleGH;
  std::vector<art::Ptr<recob::Hit> > hitlistGH;
  if (evt.getByLabel("gaushit",hitHandleGH))
    art::fill_ptr_vector(hitlistGH, hitHandleGH);
  

  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrkProducer,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  
  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = (int)plist.size();
  fData->ntrks      = (int)tracklist.size();
  fData->Resize();
 
  // flag this data as MC
  if( plist.size() ) fIsMC = true;
  
  std::cout<<"Found "<<hitlist.size()<<" hits from "<<fHitProducer<<"\n";
  std::cout<<"Found "<<tracklist.size()<<" tracks from "<<fTrkProducer<<"\n";


  //====================================
  // Keep tabs on total energy, charge,
  // and electrons arriving at anode
  //====================================
  float total_depEnergy         = 0;
  float total_depElectrons      = 0;
  float total_numElectrons      = 0;
  float total_numElectrons_2MeV = 0;


  //====================================
  // Save MCParticle information
  //====================================
  std::map<int,int> map_g4id_index;
  if( plist.size() ) {
    
    std::vector<blip::ParticleInfo> pinfo = fBlipAlg->pinfo;
    
    //fData->total_depEnergy    = 0;
    //fData->total_depElectrons = 0;
    //fData->total_numElectrons = 0;

    // Loop through the MCParticles
    if( fDebugMode ) std::cout<<"\nLooping over G4 MCParticles: \n";
    for(size_t i = 0; i<plist.size(); i++){
      auto pPart = plist[i];
      
      map_g4id_index[pPart->TrackId()] = i;

      total_depEnergy       += pinfo[i].depEnergy;
      total_depElectrons    += pinfo[i].depElectrons;
      ////total_numElectrons    += pinfo[i].numElectrons;
      //fData->total_depEnergy     += pinfo[i].depEnergy;
      //fData->total_depElectrons  += pinfo[i].depElectrons;
      //fData->total_numElectrons  += pinfo[i].numElectrons;
      
      //std::cout<<"  G4 ID "<<i<<", PDG: "<<pPart->PdgCode()<<", mother: "<<pPart->Mother()<<", KE0: "<<pinfo[i].KE<<", eDep: "<<pinfo[i].depEnergy<<"\n";
      
      // Save to TTree object
      if(i<kMaxG4){
        fData->part_trackID[i]         = pPart->TrackId();
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
        fData->part_process[i]         = pPart->Process();
        fData->part_depEnergy[i]       = pinfo[i].depEnergy;
        fData->part_depElectrons[i]    = pinfo[i].depElectrons;
        //fData->part_numElectrons[i]    = pinfo[i].numElectrons;
        fData->part_isPrimary[i]       = pinfo[i].isPrimary;
        if( fDebugMode ) PrintParticleInfo(i);
      }
    
    } // endloop over G4 particles
    
    if( fDebugMode ) std::cout<<"True total energy deposited: "<<total_depEnergy<<" MeV \n";
  
  }//endif particles found in event




  
  //====================================
  // Save TrueBlip information
  //====================================
  std::vector<blip::TrueBlip> trueblips = fBlipAlg->trueblips;
  fData->nedeps = (int)trueblips.size();
  if( trueblips.size() ) {
    //std::cout<<"Found "<<trueblips.size()<<" true blips\n";

    if( fDebugMode ) std::cout<<"\nLooping over true blips:\n";
    //for(size_t i=0; i<trueblips.size(); i++ ) {
    for(auto& trueblip : trueblips ) {
      int i = trueblip.ID;
      fData->edep_tpc[i]      = trueblip.TPC;
      fData->edep_energy[i]   = trueblip.Energy;
      fData->edep_electrons[i]= trueblip.DepElectrons;
      fData->edep_charge[i]   = trueblip.NumElectrons;
      fData->edep_x[i]        = trueblip.Position.X();
      fData->edep_y[i]        = trueblip.Position.Y();
      fData->edep_z[i]        = trueblip.Position.Z();
      fData->edep_tdrift[i]   = trueblip.DriftTime;
      fData->edep_g4id[i]     = trueblip.LeadG4ID;
      fData->edep_g4index[i]  = trueblip.LeadG4Index;
      fData->edep_pdg[i]      = trueblip.LeadG4PDG;
      
      float ne_dep  = trueblip.DepElectrons;
      float ne      = trueblip.NumElectrons;
      total_numElectrons += ne;
      if( trueblip.Energy < 2 ) total_numElectrons_2MeV += ne;
      if( fDebugMode ) PrintTrueBlipInfo(trueblip);

      // calculate simulated lifetime
      if( ne>1000 && ne<ne_dep ) {
        float tau = trueblip.DriftTime * 1./log(ne_dep/float(ne));
        h_true_lifetime->Fill(tau);
      }
    }
  }//endif trueblips were made
  
  


/*
      // ------------------------------------------------------------ 
      // Look at the dEdx of this track
      // Association between Calorimetry objects and Tracks
      art::FindManyP<anab::Calorimetry> fmcal(TrackHandle, e, fTrackCalModule);

      if( fmcal.isValid() ){

        std::vector<art::Ptr<anab::Calorimetry> > calos_mu = fmcal.at(fMuTrackIndex);

        size_t N;
        size_t plane = 0;

        for(size_t i=0; i<calos_mu.size(); i++){
          plane = calos_mu[i]->PlaneID().Plane;
          N     = calos_mu[i]->dEdx().size();
          if( plane == 1 ) fMuTrackEnergy = calos_mu[i]->KineticEnergy();
          if( plane >= 0 && plane < 2 && N > 0 ) {
            fvMuTrkResRange[plane]  = calos_mu[i]->ResidualRange();
            fvMuTrkdEdx[plane]      = calos_mu[i]->dEdx();
            fvMuTrkXYZ[plane]       = calos_mu[i]->XYZ();
            fvMuTrkPitch[plane]     = calos_mu[i]->TrkPitchVec();
            for(size_t j=0; j<N; j++){
              float T = calos_mu[i]->XYZ().at(j).X() / fDriftVelocity[0];
              float eLifetimeCorr = exp( T / fElectronLifetime );
              float dADCdx = eLifetimeCorr*calos_mu[i]->dQdx().at(j);
              float dQdx_e = fCaloAlg.ElectronsFromADCArea(dADCdx,plane);
              fvMuTrkdADCdx[plane] .push_back( dADCdx );
              fvMuTrkdQdx[plane]   .push_back( dQdx_e);
            }
          }
        }
*/


  //====================================
  // Save hit information
  //====================================
  //std::cout<<"Looping over the hits...\n";
  int   num_hits[kNplanes]            ={0};
  int   num_hits_untracked[kNplanes]  ={0};
  int   num_hits_true[kNplanes]       ={0};
  int   num_hits_pmatch[kNplanes]     ={0};
  float total_hit_charge[kNplanes]    ={0};
  
  for(size_t i=0; i<hitlist.size(); i++){
   
    int plane = hitlist[i]->WireID().Plane;
    auto const& hinfo = fBlipAlg->hitinfo[i];
    
    fNumHits[plane]++;
    num_hits[plane]++;
    if( hinfo.ismatch ) {
      fNumHitsMatched[plane]++;
      num_hits_pmatch[hinfo.plane]++;
    }
      
    if( hinfo.trkid < 0 ) {
      fNumHitsUntracked[plane]++;
      num_hits_untracked[hinfo.plane]++;
      

      // TEMP 6/8/22
      //if( ! (plane == 0 && hinfo.wire > 690 && hinfo.wire < 695 ) ) {
      
      // fill diagnostic histograms
        h_hitamp[plane]    ->Fill(hitlist[i]->PeakAmplitude());
        //h_hitsigt[plane]   ->Fill(hitlist[i]->SigmaPeakTime());
        h_hitrms[plane]    ->Fill(hitlist[i]->RMS());
        h_hitratio[plane]  ->Fill(hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
        //h_hitint[plane]    ->Fill(hitlist[i]->Integral());
        h_hitmult[plane]    ->Fill(hitlist[i]->Multiplicity());
      //}
    } 
    
    // calculate reco-true resolution
    float qcoll = hinfo.charge;
    float qtrue = hinfo.g4charge;
    if( hinfo.g4id >= 0 && qcoll>0 && qtrue>0 ) {
      fNumHitsTrue[plane]++;
      num_hits_true[plane]++;
      total_hit_charge[plane] += qtrue;
      h_nelec_TrueVsReco[plane]->Fill(qtrue/1e3,qcoll/1e3);
      h_nelec_Resolution[plane]->Fill((qcoll-qtrue)/qtrue);
    }
    
    // fill data to be saved to event tree
    if( i < kMaxHits ){
      fData->hit_channel[i]   = hitlist[i]->Channel();
      fData->hit_peakT[i]     = hitlist[i]->PeakTime();
      fData->hit_gof[i]       = hitlist[i]->GoodnessOfFit();
      fData->hit_rms[i]       = hitlist[i]->RMS();
      fData->hit_amp[i]	      = hitlist[i]->PeakAmplitude();
      fData->hit_area[i]      = hitlist[i]->Integral();
      fData->hit_sumadc[i]    = hitlist[i]->SummedADC();
      fData->hit_mult[i]      = hitlist[i]->Multiplicity();
      fData->hit_plane[i]     = hitlist[i]->WireID().Plane;
      fData->hit_wire[i]      = hitlist[i]->WireID().Wire;
      fData->hit_tpc[i]       = hitlist[i]->WireID().TPC;
      fData->hit_trkid[i]     = hinfo.trkid;
      fData->hit_time[i]      = hinfo.driftTime;
      fData->hit_charge[i]    = hinfo.charge;
      fData->hit_ismatch[i]   = hinfo.ismatch;
      fData->hit_g4id[i]      = hinfo.g4id;
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
    
    // fill nhit diagnostic histograms
    h_nhits[ip]   ->Fill(num_hits[ip]);
    h_nhits_ut[ip]->Fill(num_hits_untracked[ip]);
    h_nhits_m[ip] ->Fill(num_hits_pmatch[ip]);
    h_nhits_tm[ip]->Fill(num_hits_true[ip]);
    //std::cout<<"* plane "<<ip<<": "<<num_hits[ip]<<" hits";
  
    // calculate overall hit purity/completeness per plane
    float qcomp     = -9;
    float pur       = -9;
    if( num_hits_true[ip] ) {
      if(total_numElectrons )  qcomp = total_hit_charge[ip]/total_numElectrons;
      if(num_hits[ip]       )  pur   = num_hits_true[ip]/float(num_hits[ip]);
      h_chargecomp[ip]->Fill( qcomp );
      h_hitpur[ip]    ->Fill( pur );
      //std::cout<<"  ("<<num_hits_true[ip]<<" truth-matched) -- completenes: "<<qcomp<<", purity "<<pur;
    }
    //std::cout<<"\n";
  }//endloop over planes
    
  //if( hitlist.size() ){
  //  std::cout<<nhits_untracked<<" hits not in 3D tracks \n";
    //std::cout<<"("<<100*nhits_untracked/float(hitlist.size())<<"% of hits from "<<fHitProducer<<")\n";
  //}
 

  
  //====================================
  // Save track information
  //====================================
  //std::cout<<"Looping over tracks...\n";
  h_ntrks->Fill(tracklist.size());
  for(size_t i=0; i<tracklist.size(); i++){
    auto& trk = tracklist[i];
    const auto& startPt = trk->Vertex();
    const auto& endPt   = trk->End();
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
    h_trk_length        ->Fill(trk->Length());
    h_trk_xspan         ->Fill( fabs(startPt.X()-endPt.X()) );
  }//endloop over trks
  
    
   
  //=============================================
  // Save hit cluster info
  //=============================================
  fData->nclusts = (int)fBlipAlg->hitclust.size();
  //std::cout<<"Looping over the clusts...\n";
  if( fDebugMode ) std::cout<<"\nLooping over clusters...\n";
  for(size_t i=0; i < fBlipAlg->hitclust.size(); i++){
    auto const& clust = fBlipAlg->hitclust[i];
    fData->clust_id[i]        = clust.ID;
    fData->clust_tpc[i]       = clust.TPC;
    fData->clust_plane[i]     = clust.Plane;
    fData->clust_wire[i]      = clust.CentHitWire;
    fData->clust_startwire[i] = clust.StartWire;
    fData->clust_endwire[i]   = clust.EndWire;
    
    if( !fSavePlaneInfo[clust.Plane] ) continue;
    fData->clust_nwires[i]    = clust.NWires;
    fData->clust_nhits[i]     = clust.NHits;
    
    //fData->clust_time[i]      = Truncate(clust.Time);
    //fData->clust_timespan[i]  = Truncate(clust.Timespan);
    //fData->clust_rms[i]       = Truncate(clust.TimeErr,0.01);
    //fData->clust_starttime[i] = Truncate(clust.StartTime);
    //fData->clust_endtime[i]   = Truncate(clust.EndTime);
    //fData->clust_charge[i]    = Truncate(clust.Charge,10);
    //fData->clust_ismatch[i]   = clust.isMatched;
    //fData->clust_blipid[i]    = clust.BlipID;
    //fData->clust_amp[i]       = Truncate(clust.Amplitude,0.01);
    fData->clust_time[i]      = clust.Time;
    fData->clust_timespan[i]  = clust.Timespan;
    fData->clust_rms[i]       = clust.TimeErr;
    fData->clust_starttime[i] = clust.StartTime;
    fData->clust_endtime[i]   = clust.EndTime;
    fData->clust_charge[i]    = clust.Charge;
    fData->clust_amp[i]       = clust.Amplitude;
    fData->clust_ismatch[i]   = clust.isMatched;
    fData->clust_blipid[i]    = clust.BlipID;

    //std::cout<<"hit clust on plane "<<clust.Plane<<"   NWires "<<clust.NWires<<"    NHits "<<clust.NHits<<"   lHitRMS "<<clust.LeadHit->RMS()<<"   Span "<<clust.Timespan<<"   TimeErr "<<clust.TimeErr<<"  "<<clust.TimeErr/clust.Timespan<<"\n";
      
    // if this clust has an associated "trueblip" ID, find it
    // and figure out the true G4 charge, energy, etc
    int tbi = clust.EdepID;
    if( tbi >= 0 && tbi < (int)fBlipAlg->trueblips.size() ) {
      auto const& trueBlip = fBlipAlg->trueblips[tbi];
      fData->edep_clustid[tbi] = clust.ID;
      fData->clust_edepid[i]   = trueBlip.ID;
      fData->clust_g4energy[i] = trueBlip.Energy; 
      fData->clust_g4charge[i] = trueBlip.NumElectrons;
      // fill histograms of electron/alpha charge resolution,
      // also derive the electron-to-ADC factor (this ~should~ 
      // match up with CalAreaConstants!)
      float q   = clust.Charge;
      float qt  = trueBlip.NumElectrons;
      int pdg   = trueBlip.LeadG4PDG;
      if( q && qt && clust.Plane == 2 ) {
          if( fabs(pdg) == 11 )         h_qres_electrons->Fill( (q-qt)/qt );
          if( fabs(pdg) == 1000020040 ) h_qres_alphas   ->Fill( (q-qt)/qt );
          h_adc_factor->Fill( clust.ADCs / qt );
      }
    
    }
      
    //if( fDebugMode ) PrintClusterInfo(clust);
    
  }//endloop over 2D hit clusters

  
//====================================
  // Save blip info to tree
  //===================================
  fData->nblips             = fBlipAlg->blips.size();
  int nblips_matched        = 0;
  int nblips_total          = 0;
  int nblips_picky          = 0;
  float true_blip_charge    = 0;
  float true_blip_charge_2MeV = 0;
  //fData->total_blip_energy  = 0;
  //std::cout<<"Looping over the blips...\n";
  for(size_t i=0; i<fBlipAlg->blips.size(); i++){
    auto& b = fBlipAlg->blips[i];
  
    nblips_total++;
    fNum3DBlips++;
    if( b.NPlanes >= 3 ) fNum3DBlips3Plane++;


    fData->blip_id[i]         = i;
    fData->blip_tpc[i]        = b.TPC;
    fData->blip_nplanes[i]    = b.NPlanes;
    fData->blip_x[i]          = b.X;
    fData->blip_y[i]          = b.Y;
    fData->blip_z[i]          = b.Z;
    fData->blip_sigmayz[i]    = b.SigmaYZ;
    fData->blip_dx[i]         = b.dX;
    fData->blip_dyz[i]        = b.dYZ;
    fData->blip_size[i]       = sqrt( pow(b.dX,2) + pow(b.dYZ,2) );
    fData->blip_trkid[i]      = b.TrkID;
    fData->blip_trklength[i]  = b.Length;
    fData->blip_time[i]       = b.Time;
    fData->blip_proxtrkdist[i]= b.ProxTrkDist;
    fData->blip_proxtrkid[i]  = b.ProxTrkID;
    fData->blip_incylinder[i] = b.inCylinder;
   
    for(size_t ipl = 0; ipl<kNplanes; ipl++){
      if( b.clusters[ipl].NHits <= 0 ) continue;
      fData->blip_clustid[ipl][i]        = b.clusters[ipl].ID;
      /*
      fData->blip_clust_wire[ipl][i]      = b.clusters[ipl].LeadHitWire;
      fData->blip_clust_startwire[ipl][i] = b.clusters[ipl].StartWire;
      fData->blip_clust_endwire[ipl][i]   = b.clusters[ipl].StartWire;
      fData->blip_clust_nhits[ipl][i]     = b.clusters[ipl].NHits;
      fData->blip_clust_nwires[ipl][i]    = b.clusters[ipl].NWires;
      fData->blip_clust_charge[ipl][i]    = Truncate(b.clusters[ipl].Charge,10);
      fData->blip_clust_lhit_rms[ipl][i]  = Truncate(b.clusters[ipl].LeadHit->RMS());
      fData->blip_clust_lhit_amp[ipl][i]  = Truncate(b.clusters[ipl].LeadHit->PeakAmplitude(),0.001);
      */
    }
    
    if( b.clusters[fCaloPlane].Charge > 0 ) {
      fData->blip_charge[i]     = b.clusters[fCaloPlane].Charge;
      fData->blip_sumadc[i]     = b.clusters[fCaloPlane].ADCs;
      fData->blip_energy[i]     = b.Energy;
      h_blip_charge             ->Fill(b.clusters[fCaloPlane].Charge);
      h_blip_sumadc             ->Fill(b.clusters[fCaloPlane].ADCs);
    }

    // Fill cluster charge 2D histograms
    h_blip_charge_YU->Fill( 0.001*b.clusters[2].Charge, 0.001*b.clusters[0].Charge );
    h_blip_charge_YV->Fill( 0.001*b.clusters[2].Charge, 0.001*b.clusters[1].Charge );
    h_blip_charge_UV->Fill( 0.001*b.clusters[0].Charge, 0.001*b.clusters[1].Charge );

    // Select picky (high-quality) blips:
    if(b.NPlanes == 3 && b.SigmaYZ < 1.) {
      nblips_picky++;
      fNum3DBlipsPicky++;
      h_blip_charge_picky ->Fill(b.clusters[fCaloPlane].Charge);
       h_blip_zy_picky     ->Fill(b.Z, b.Y);
      
      h_blip_charge_YU_picky->Fill( 0.001*b.clusters[2].Charge, 0.001*b.clusters[0].Charge );
      h_blip_charge_YV_picky->Fill( 0.001*b.clusters[2].Charge, 0.001*b.clusters[1].Charge );
      h_blip_charge_UV_picky->Fill( 0.001*b.clusters[0].Charge, 0.001*b.clusters[1].Charge );
    }

    h_blip_zy     ->Fill(b.Z, b.Y);
    h_blip_nplanes->Fill(b.NPlanes);
   
    // -----------------------------------------------
    // save the clustIDs and true energy deposits to the blip
    // (use the association between clust <--> edep)
    // -----------------------------------------------
    if( b.truth.ID >= 0 ) {
     
      fData->blip_edepid[i]           = b.truth.ID;
      fData->edep_blipid[b.truth.ID]  = b.ID;
      
      fNum3DBlipsTrue++;
      nblips_matched++;
      true_blip_charge+= b.truth.NumElectrons;
      if( b.truth.Energy < 2 ) true_blip_charge_2MeV += b.truth.NumElectrons;
      h_blip_reszy->Fill( b.Position.Z()-b.truth.Position.Z(), b.Position.Y()-b.truth.Position.Y() );
      h_blip_resx->Fill( b.Position.X()-b.truth.Position.X() );
    }
  
  }//endloop over 3D blips
 
  // Fill some more histograms...
  h_nblips->Fill(nblips_total);
  h_nblips_picky->Fill(nblips_picky);
  if( fIsMC ) {
    h_nblips_tm->Fill(nblips_matched);
    //if( fData->total_numElectrons ) h_blip_qcomp      ->Fill(true_blip_charge / fData->total_numElectrons);
    if( total_numElectrons        ) h_blip_qcomp      ->Fill(true_blip_charge      / total_numElectrons     );
    if( total_numElectrons_2MeV   ) h_blip_qcomp_2MeV ->Fill(true_blip_charge_2MeV / total_numElectrons_2MeV);
    if( nblips_total ) h_blip_pur    ->Fill( float(nblips_matched) / float(nblips_total) );
  }

  std::cout<<"Reconstructed "<<fBlipAlg->blips.size()<<" 3D blips";
  if( nblips_matched ) std::cout<<"; "<<nblips_matched<<" were truth-matched";
  std::cout<<"\n";

  if( fDebugMode ) {
    for(auto const& b : fBlipAlg->blips ) {
      std::cout
      <<"   -- "<<b.ID<<", TPC: "<<b.TPC
      <<"; charge: "<<b.clusters[2].Charge
      <<"; recoEnergy: "<<b.Energy<<" MeV"
      <<"; Position: "<<b.Position.X()<<", "<<b.Position.Y()<<", "<<b.Position.Z()
      //<<"; MaxIntersectDiff: "<<b.MaxIntersectDiff
      <<"; EdepID: "<<fData->blip_edepid[b.ID]
      <<"\n";
    }
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

  float nEvents = float(fNumEvents);
  
  printf("\n***********************************************\n");
  fBlipAlg->PrintConfig();
  printf("BlipAna Summary\n\n");
  printf("  Total events                : %i\n",            fNumEvents);
  printf("  Blips per evt, total        : %.3f\n",          fNum3DBlips/nEvents);
  //printf("                 3 planes   : %.3f\n",          fNum3DBlips3Plane/nEvents);
  printf("                 picky        : %.3f\n",      fNum3DBlipsPicky/nEvents);
  printf("                 picky frac   : %5.3f\n",     fNum3DBlipsPicky/float(fNum3DBlips));
  //printf("  Matched planes per blip   : %.3f\n", h_blip_nplanes->GetMean());
  if(fIsMC){
  printf("  MC-matched blips per evt    : %.3f\n",       fNum3DBlipsTrue/nEvents);
  printf("  Charge completeness, total  : %.4f +/- %.4f\n", h_blip_qcomp->GetMean(), h_blip_qcomp->GetStdDev()/sqrt(fNumEvents));
  printf("                       < 2MeV : %.4f +/- %.4f\n", h_blip_qcomp_2MeV->GetMean(), h_blip_qcomp_2MeV->GetStdDev()/sqrt(fNumEvents));
  printf("  Blip purity                 : %.4f\n",       h_blip_pur->GetMean());
  }
  printf("  Mean blip charge            : %.0f\n",      h_blip_charge->GetMean());
  printf("\n");
  for(size_t i=0; i<kNplanes; i++){
  printf("  Plane %lu -------------------------\n",i);
  printf("   * total hits/evt           : %.2f\n",fNumHits[i]/(float)fNumEvents);
  printf("   * untracked hits/evt       : %.2f\n",fNumHitsUntracked[i]/(float)fNumEvents);
  printf("   * plane-matched hits/evt   : %.2f\n",fNumHitsMatched[i]/(float)fNumEvents);
  if(fIsMC) {
  printf("   * true-matched hits/evt    : %.2f\n",fNumHitsTrue[i]/(float)fNumEvents);
  printf("   * charge completeness      : %.4f\n",h_chargecomp[i]->GetMean());
  printf("   * hit purity               : %.4f\n",h_hitpur[i]->GetMean());
  }
  } 
  printf("\n***********************************************\n");
 
}




//###################################################
//  Printouts for debugging
//###################################################

void BlipAna::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%8.3f, Edep=%8.3f, T=%10.2f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->part_trackID[i],
   fData->part_pdg[i],
   fData->part_startPointx[i],
   fData->part_startPointy[i],
   fData->part_startPointz[i],
   fData->part_pathlen[i], 
   fData->part_KE[i],
   fData->part_depEnergy[i],
   fData->part_startT[i]/1e3,
   fData->part_mother[i],
   fData->part_process[i].c_str(),
   fData->part_nDaughters[i]
  ); 
}

void BlipAna::PrintTrueBlipInfo(const blip::TrueBlip& tb){
  printf("  %5i  G4ID: %-6i PDG: %-10i XYZ: %7.2f, %7.2f, %7.2f, %8.3f MeV, %8i e- deposited, %8i e- @anode\n",
   tb.ID,
   tb.LeadG4ID,
   tb.LeadG4PDG,
   tb.Position.X(),
   tb.Position.Y(),
   tb.Position.Z(),
   tb.Energy,
   tb.DepElectrons,
   tb.NumElectrons
  ); 
}

void BlipAna::PrintHitInfo(const blip::HitInfo& hi){
  printf("  hitID: %4i, TPC: %i, plane: %i, driftTicks: %7.2f, leadWire: %3i, G4ID: %4i, recoTrack: %4i\n",
    hi.hitid,
    hi.tpc,
    hi.plane,
    hi.driftTime,
    hi.wire,
    hi.g4id,
    hi.trkid
  );
}

void BlipAna::PrintClusterInfo(const blip::HitClust& hc){
  printf("  ID: %4i, TPC: %i, plane: %i, time range: %7.2f - %7.2f, timespan: %6.2f, leadWire: %3i, nwires: %3i, nhits: %3i, edepid: %i, isMatched: %i, trkID: %i\n",
    hc.ID,
    hc.TPC,
    hc.Plane,
    hc.StartTime,
    hc.EndTime,
    hc.Timespan,
    hc.CentHitWire,
    hc.NWires,
    hc.NHits,
    hc.EdepID,
    hc.isMatched,
    hc.TrkID
  );
}

/*
float BlipAna::Truncate(float input, double base){
  if( base > 0 )  return roundf( input / base ) * base;
  else            return input;
  //return input;
}

double BlipAna::Truncate(double input, double base){
  if( base > 0 )  return roundf( input / base ) * base;
  else            return input;
  //return input;
}
*/


DEFINE_ART_MODULE(BlipAna)

#endif
