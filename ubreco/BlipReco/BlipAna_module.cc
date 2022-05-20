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
//#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
//#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
//#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "cetlib/search_path.h"

// MicroBooNE-specific includes
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

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
//#include "TGraph2D.h"
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
const int kMaxHits  = 50000;
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
  TTree* blipTree;

  blip::Blip thisBlip;

  // --- Configurations and switches ---
  std::string treeName      = "anatree";
  std::string blipTreeName  = "bliptree";
  bool  saveTruthInfo       = true;
  bool  saveTrkInfo         = true;
  bool  saveHitInfo         = true;
  bool  saveClustInfo       = true;
  bool  saveEvtTree         = true;
  bool  saveBlipTree        = true;
  
  int   caloPlane           = 2;

  // --- Event information ---   
  int           event;                    // event number
  int           run;                      // run number
  unsigned int  timestamp;                // unix time of event
  float         lifetime;                 // electron lifetime
  
  // --- G4 information ---
  int   nparticles;               // number of G4 particles
  bool  isPrimary[kMaxG4];        // is primary particle
  int   trackID[kMaxG4];          // G4 track ID
  int   pdg[kMaxG4];              // PDG
  int   nDaughters[kMaxG4];       // number of daughters
  int   mother[kMaxG4];           // mother particle
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
  float depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  int   depElectrons[kMaxG4];     // electrons deposited
  float numElectrons[kMaxG4];     // electrons reaching anode wires
  std::vector<std::string> process;// process name
  float total_depEnergy;          // total deposited energy in AV
  int   total_depElectrons;       // total deposited ionization electrons in AV
  float total_numElectrons;       // total electrons reaching anode wires

  // --- True energy deposit info (derived from SimChannels and SimEnergyDeposits) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 track ID
  int   edep_g4index[kMaxEDeps];  // leading G4 track index
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_clustid[kMaxEDeps];  // hitclust ID
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
  int   edep_depne[kMaxEDeps];    // total ionization electrons deposited
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
  float clust_lhit_gof[kMaxHits];
  bool  clust_lhit_isfit[kMaxHits];   // is there a valid goodness of fit for lead hit?
  bool  clust_ismatch[kMaxHits];      // was this cluster plane-matched?
  float clust_sumadc[kMaxHits];          // summed ADC of hits in cluster
  float clust_charge[kMaxHits];       // cluster charge at anode [e-]
  float clust_time[kMaxHits];         // cluster time-tick
  float clust_time_err[kMaxHits];     // cluster time uncertainty
  float clust_startTime[kMaxHits];    // cluster start tick
  float clust_endTime[kMaxHits];      // cluster end tick
  float clust_timespan[kMaxHits];     // cluster timespan
  float clust_g4energy[kMaxHits];     // true cluster energy from G4
  float clust_g4charge[kMaxHits];     // true cluster charge at anode
  int   clust_g4id[kMaxHits];         // true MCParticle ID (index for particle branches)
  int   clust_blipid[kMaxHits];       // blip ID for this nlusteer (if it was made into one)
  int   clust_edepid[kMaxHits];       // true energy dep ID

  // --- 3D Blip information ---
  float total_blip_energy;            // total summed blip energy in event [MeV]
  short nblips;                       // number of blips in event
  short blip_id[kMaxBlips];           // blip ID / index
  short blip_tpc[kMaxBlips];          // blip TPC
  short blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_maxdiff[kMaxBlips];      // difference in wire intersection points
  float blip_sumadc[kMaxBlips];          // integrated ADCs 
  float blip_charge[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip energy [MeV]
  int   blip_edepid[kMaxBlips];       // true energy dep ID
  float blip_trkdist[kMaxBlips];      // distance to nearest track
  short blip_trkid[kMaxBlips];        // index of nearest trk
  bool  blip_incylinder[kMaxBlips];   // is blip within a cylinder near a track
  int   blip_clustid[kNplanes][kMaxBlips];  // cluster ID per plane
  short blip_nhits[kNplanes][kMaxBlips];  // hits per cluster on each plane
  short blip_nwires[kNplanes][kMaxBlips]; // wires per cluster on each plane
  float blip_timespan[kNplanes][kMaxBlips]; // cluster timespan on each plane [ticks]
  //float blip_matchscore[kNplanes][kMaxBlips]; // cluster matchscore on each plane
  
  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    lifetime              = -999;
    timestamp             = -999;
    total_depEnergy       = -999;
    total_depElectrons    = -999;
    total_numElectrons    = -999;
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
    FillWith(depElectrons,-999);
    FillWith(numElectrons,-999);
    FillWith(depEnergy,   -999.);
    FillWith(process,     "");
    nedeps                = 0;    // --- EDeps ---
    FillWith(edep_tpc,    -9);
    FillWith(edep_energy, -999);
    FillWith(edep_depne,  -999);
    FillWith(edep_charge, -999);
    FillWith(edep_x,      -99999.);
    FillWith(edep_y,      -99999.);
    FillWith(edep_z,      -99999.);
    FillWith(edep_ds,     -999);
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
    FillWith(clust_wire,      -9);
    FillWith(clust_startwire, -9);
    FillWith(clust_endwire,   -9);
    FillWith(clust_chan,      -9);
    FillWith(clust_nwires,    -9);
    FillWith(clust_nhits,     -9);
    FillWith(clust_lhit_id,   -9);
    FillWith(clust_lhit_amp,  -9);
    FillWith(clust_lhit_rms,  -9);
    FillWith(clust_lhit_time, -9);
    FillWith(clust_lhit_peakT,-9);
    FillWith(clust_lhit_gof,  -9);
    FillWith(clust_lhit_isfit,false);
    FillWith(clust_sumadc,       -9);
    FillWith(clust_charge,    -999);
    FillWith(clust_time,      -999);
    FillWith(clust_time_err,  -999);
    FillWith(clust_startTime, -999);
    FillWith(clust_endTime,   -999);
    FillWith(clust_timespan,  -9);
    FillWith(clust_g4id,      -9);
    FillWith(clust_g4charge,  -999);
    FillWith(clust_g4energy,  -9);
    FillWith(clust_ismatch,   false);
    FillWith(clust_edepid,    -9);
    FillWith(clust_blipid,    -9);
    total_blip_energy         = -9;  // --- Blips ---
    nblips                    = 0;
    FillWith(blip_id,         -9);
    FillWith(blip_tpc,        -9);
    FillWith(blip_nplanes,    -9);
    FillWith(blip_x,          -9999);
    FillWith(blip_y,          -9999);
    FillWith(blip_z,          -9999);
    FillWith(blip_maxdiff,    -99);
    FillWith(blip_sumadc,        -999);
    FillWith(blip_charge,     -999);
    FillWith(blip_energy,     -999);
    FillWith(blip_trkdist,    -99);
    FillWith(blip_trkid,      -9);
    FillWith(blip_incylinder, false);
    FillWith(blip_edepid,     -9);
    for(int i=0; i<kNplanes; i++) {
      FillWith(blip_clustid[i],-9);
      FillWith(blip_nhits[i],-9);
      FillWith(blip_nwires[i],-9);
      FillWith(blip_timespan[i],-9);
      //FillWith(blip_matchscore[i],-9);
    }
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
   
    if( saveEvtTree ) {
    
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
        evtTree->Branch("trk_id",trk_id,"trk_id[ntrks]/S");       
        evtTree->Branch("trk_length",trk_length,"trk_length[ntrks]/F");
        evtTree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
        evtTree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
        evtTree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
        evtTree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
        evtTree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
        evtTree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
      }

      if( saveClustInfo ) {
        evtTree->Branch("nclusts",&nclusts,"nclusts/I");
        //evtTree->Branch("clust_id",&clust_id,"clust_id[nclusts]/I");
        evtTree->Branch("clust_plane",clust_plane,"clust_plane[nclusts]/S");
        evtTree->Branch("clust_wire",clust_wire,"clust_wire[nclusts]/S");
        evtTree->Branch("clust_startwire",clust_startwire,"clust_startwire[nclusts]/S");
        evtTree->Branch("clust_endwire",clust_endwire,"clust_endwire[nclusts]/S");
        evtTree->Branch("clust_nwires",clust_nwires,"clust_nwires[nclusts]/S");
        evtTree->Branch("clust_nhits",clust_nhits,"clust_nhits[nclusts]/S");
        //evtTree->Branch("clust_lhit_amp",clust_lhit_amp,"clust_lhit_amp[nclusts]/F");
        //evtTree->Branch("clust_lhit_rms",clust_lhit_rms,"clust_lhit_rms[nclusts]/F");
        //evtTree->Branch("clust_lhit_gof",clust_lhit_gof,"clust_lhit_gof[nclusts]/F");
        //evtTree->Branch("clust_lhit_isfit",clust_lhit_isfit,"clust_lhit_isfit[nclusts]/O");
        //evtTree->Branch("clust_sumadc",clust_sumadc,"clust_sumadc[nclusts]/F");
        evtTree->Branch("clust_charge",clust_charge,"clust_charge[nclusts]/F");
        evtTree->Branch("clust_time",clust_time,"clust_time[nclusts]/F");
        //evtTree->Branch("clust_timespan",clust_timespan,"clust_timespan[nclusts]/F");
        if( saveTruthInfo ) {
          evtTree->Branch("clust_g4charge",clust_g4charge,"clust_g4charge[nclusts]/F");
          evtTree->Branch("clust_g4energy",clust_g4energy,"clust_g4energy[nclusts]/F");
          evtTree->Branch("clust_edepid",clust_edepid,"clust_edepid[nclusts]/I");
        }
        evtTree->Branch("clust_ismatch",clust_ismatch,"clust_ismatch[nclusts]/O");
        evtTree->Branch("clust_blipid",clust_blipid,"clust_blipid[nclusts]/I");
      }

      evtTree->Branch("nblips",&nblips,"nblips/S");
      //evtTree->Branch("blip_tpc",blip_tpc,"blip_tpc[nblips]/S");
      //evtTree->Branch("blip_id",blip_id,"blip_id[nblips]/S");
      evtTree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/S");
      evtTree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
      evtTree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
      evtTree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
      evtTree->Branch("blip_maxdiff",blip_maxdiff,"blip_maxdiff[nblips]/F");
      //evtTree->Branch("blip_sumadc",blip_sumadc,"blip_sumadc[nblips]/F");
      evtTree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
      evtTree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
      evtTree->Branch("blip_trkdist",blip_trkdist,"blip_trkdist[nblips]/F");
      evtTree->Branch("blip_trkid",blip_trkid,"blip_trkid[nblips]/S");
      evtTree->Branch("blip_incylinder",blip_incylinder,"blip_incylinder[nblips]/O");
      if( saveTruthInfo ) evtTree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
      for(int i=0;i<kNplanes;i++) {
        evtTree->Branch(Form("blip_pl%i_clustid",i),blip_clustid[i],Form("blip_pl%i_clustid[nblips]/I",i));
        //if( i != caloPlane ) 
          //evtTree->Branch(Form("blip_pl%i_matchscore",i),blip_matchscore[i], Form("blip_pl%i_matchscore[nblips]/F",i));
        //if( i != caloPlane ) continue;
        //evtTree->Branch(Form("blip_pl%i_nhits",i),blip_nhits[i],    Form("blip_pl%i_nhits[nblips]/S",i));
        //evtTree->Branch(Form("blip_pl%i_nwires",i),blip_nwires[i],  Form("blip_pl%i_nwires[nblips]/S",i));
        //evtTree->Branch(Form("blip_pl%i_timespan",i),blip_timespan[i],  Form("blip_pl%i_timespan[nblips]/F",i));
      }
      
      //evtTree->Branch("total_blip_energy",&total_blip_energy,"total_blip_energy/F");
      
      if( saveTruthInfo ) {
        evtTree->Branch("total_depEnergy",&total_depEnergy,"total_depEnergy/F");
        //evtTree->Branch("total_depElectrons",&total_depElectrons,"total_depElectrons/I");
        evtTree->Branch("total_numElectrons",&total_numElectrons,"total_numElectrons/F");
        evtTree->Branch("nparticles",&nparticles,"nparticles/S");
        evtTree->Branch("isPrimary",isPrimary,"isPrimary[nparticles]/O");
        evtTree->Branch("trackID",trackID,"trackID[nparticles]/S");
        evtTree->Branch("pdg",pdg,"pdg[nparticles]/I");
        evtTree->Branch("nDaughters",nDaughters,"nDaughters[nparticles]/S");
        evtTree->Branch("mother",mother,"mother[nparticles]/S");
        evtTree->Branch("KE",KE,"KE[nparticles]/F");
        evtTree->Branch("endKE",endKE,"endKE[nparticles]/F");
        //evtTree->Branch("mass",mass,"mass[nparticles]/F");
        evtTree->Branch("P",P,"P[nparticles]/F");
        //evtTree->Branch("Px",Px,"Px[nparticles]/F");
        //evtTree->Branch("Py",Py,"Py[nparticles]/F");
        //evtTree->Branch("Pz",Pz,"Pz[nparticles]/F");
        evtTree->Branch("startPointx",startPointx,"startPointx[nparticles]/F");
        evtTree->Branch("startPointy",startPointy,"startPointy[nparticles]/F");
        evtTree->Branch("startPointz",startPointz,"startPointz[nparticles]/F");
        evtTree->Branch("endPointx",endPointx,"endPointx[nparticles]/F");
        evtTree->Branch("endPointy",endPointy,"endPointy[nparticles]/F");
        evtTree->Branch("endPointz",endPointz,"endPointz[nparticles]/F");
        evtTree->Branch("startT",startT,"startT[nparticles]/F");
        evtTree->Branch("endT",endT,"endT[nparticles]/F");
        evtTree->Branch("pathlen",pathlen,"pathlen[nparticles]/F");
        evtTree->Branch("depEnergy",depEnergy,"depEnergy[nparticles]/F");
        evtTree->Branch("depElectrons",depElectrons,"depElectrons[nparticles]/I");
        evtTree->Branch("numElectrons",numElectrons,"numElectrons[nparticles]/F");
        evtTree->Branch("process",&process);
        
        evtTree->Branch("nedeps",&nedeps,"nedeps/I");
        evtTree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
        evtTree->Branch("edep_g4index",edep_g4index,"edep_g4index[nedeps]/I"); 
        evtTree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
        evtTree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
        evtTree->Branch("edep_clustid",edep_clustid,"edep_clustid[nedeps]/I"); 
        evtTree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
        evtTree->Branch("edep_depne",edep_depne,"edep_depne[nedeps]/I"); 
        evtTree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/F"); 
        evtTree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
        evtTree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
        evtTree->Branch("edep_z",edep_z,"edep_z[nedeps]/F"); 
        evtTree->Branch("edep_ds",edep_ds,"edep_ds[nedeps]/F"); 
      }
    }
    
    if( saveBlipTree ) {
      blipTree = tfs->make<TTree>(blipTreeName.c_str(), "blip tree");
      blipTree->Branch("event",       &event,           "event/I");
      blipTree->Branch("run",         &run,             "run/I");
      blipTree->Branch("timestamp",   &timestamp,       "timestamp/i");
      //blipTree->Branch("lifetime",    &lifetime,        "lifetime/F");
      //blipTree->Branch("blip_tpc",    &thisBlip.TPC,    "blip_tpc/S");
      blipTree->Branch("blip_nplanes",&thisBlip.NPlanes,"blip_nplanes/S");
      blipTree->Branch("blip_maxdiff",&thisBlip.MaxIntersectDiff,"blip_maxdiff/F");
      blipTree->Branch("blip_x",      &thisBlip.X,      "blip_x/F");
      blipTree->Branch("blip_y",      &thisBlip.Y,      "blip_y/F");
      blipTree->Branch("blip_z",      &thisBlip.Z,      "blip_z/F");
      //blipTree->Branch("blip_charge", &thisBlip.Charge[2],"blip_charge/F");
      blipTree->Branch("blip_energy", &thisBlip.Energy, "blip_energy/F");
      blipTree->Branch("blip_trkdist",&thisBlip.TrkDist,"blip_trkdist/F");
      //blipTree->Branch("blip_trkid",  &thisBlip.trkid,  "blip_trkid/S");
      blipTree->Branch("blip_incylinder",&thisBlip.inCylinder,"blip_incylinder/O");
      for(int i=0; i<kNplanes; i++) {
        blipTree->Branch(Form("blip_pl%i_clustid",i), &thisBlip.ClustID[i], Form("blip_pl%i_clustid/I",i));
        if( i != caloPlane ) continue;
        blipTree->Branch(Form("blip_pl%i_nhits",i),   &thisBlip.NHits[i],   Form("blip_pl%i_nhits/S",i));
        blipTree->Branch(Form("blip_pl%i_nwires",i),  &thisBlip.NWires[i],  Form("blip_pl%i_nwires/S",i));
        blipTree->Branch(Form("blip_pl%i_timespan",i),  &thisBlip.Timespan[i],  Form("blip_pl%i_timespan/F",i));
      }
    }
   
  }//end MakeTree
  
  // === Function for filling tree ===
  void FillBlipTree( const blip::Blip& tb ) {
    if( saveBlipTree ) {
      thisBlip = tb;
      blipTree->Fill();
    }
  }
  
  // === Function for filling event tree ===
  void FillEventTree(){
    if( saveEvtTree ) evtTree->Fill();
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
  void PrintParticleInfo(size_t);
  void PrintG4Info(const simb::MCParticle&);
  void PrintClusterInfo(const blip::HitClust&);
  void PrintHitInfo(const blip::HitInfo&);

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
  int   fNum3DBlipsTrue     = 0;

  // --- Histograms ---
  TH1D*   h_nhits[kNplanes];
  TH1D*   h_nhits_m[kNplanes];
  TH1D*   h_nhits_tm[kNplanes];
  TH1D*   h_hitamp[kNplanes];
  TH1D*   h_hitsigt[kNplanes];
  TH1D*   h_hitrms[kNplanes];
  TH1D*   h_hitratio[kNplanes];
  TH1D*   h_hitint[kNplanes];
  TH2D*   h_nelec_TrueVsReco[kNplanes];
  TH1D*   h_nelec_Resolution[kNplanes];
  TH1D*   h_chargecomp[kNplanes];
  TH1D*   h_hitpur[kNplanes];
  TH1D*   h_clust_nwires;
  TH1D*   h_clust_timespan;
  TH1D*   h_hit_dt[kNplanes];
  TH1D*   h_hit_dtfrac[kNplanes];
  TH1D*   h_nmatches[3];
  //TH1D*   h_clust_matchScore[kNplanes];
  //TH1D*   h_clust_matchScore_3D[kNplanes];
  TH1D*   h_nblips;
  TH1D*   h_nblips_tm;
  TH1D*   h_blip_nplanes;
  TH1D*   h_blip_qcomp;
  TH1D*   h_blip_pur;
  TH2D*   h_blip_reszy;
  TH1D*   h_blip_resx;
  TH2D*   h_blip_zy;
  TH2D*   h_blip_zy_3p;
  TH1D*   h_blip_sumadc;
  TH1D*   h_blip_charge;
  TH1D*   h_blip_charge_3p;

  // Some truth metrics for debugging
  TH1D*   h_qres_electrons;
  TH1D*   h_qres_alphas;
  TH1D*   h_adc_factor;
  TH1D*   h_alpha_qdep;


  // Initialize histograms
  void InitializeHistograms(){

    art::ServiceHandle<art::TFileService> tfs;
  
    float blipMax = 500;    int blipBins  = 500;
    h_nblips          = tfs->make<TH1D>("nblips","Reconstructed 3D blips per event",blipBins,0,blipMax);
    h_blip_zy         = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy         ->SetOption("COLZ");
    h_blip_zy_3p         = tfs->make<TH2D>("blip_zy_3plane","3D blip location (3-plane match);Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy_3p         ->SetOption("COLZ");
      
    art::TFileDirectory dir_diag = tfs->mkdir("Diagnostics");
    h_blip_nplanes    = dir_diag.make<TH1D>("blip_nplanes","Matched planes per blip",3,1,4);
    h_blip_sumadc     = dir_diag.make<TH1D>("blip_sumadc","3D blips;Integral on coll plane [ADC]",            200,0,500);
    h_blip_charge     = dir_diag.make<TH1D>("blip_charge","3D blips;Charge [e-]",                             200,0,100e3);
    h_blip_charge_3p  = dir_diag.make<TH1D>("blip_charge_3plane","3D blips w/match on ALL planes;Charge [e-]",200,0,100e3);
    
    float hitMax  = 10000; int hitBins  = 1000;
    float ampMax  = 50;   int ampBins   = 500;
    float rmsMax  = 30;   int rmsBins   = 600;
    float areaMax = 1000; int areaBins  = 1000;
    float ratioMax= 5.0;  int ratioBins = 250;
    
    // MC histograms related to truth
    art::TFileDirectory dir_truth = dir_diag.mkdir("Truth");
    h_nblips_tm     = dir_truth.make<TH1D>("nblips_tm","Truth-matched 3D blips per event",blipBins,0,blipMax);
    h_blip_qcomp    = dir_truth.make<TH1D>("blip_qcomp","Fraction of true charge (at anode) reconstructed into 3D blips",202,0,1.01);
    h_blip_pur      = dir_truth.make<TH1D>("blip_pur","Fraction of truth-matched blips",202,0,1.01);
    h_blip_reszy    = dir_truth.make<TH2D>("blip_res_zy","Blip position resolution;Z_{reco} - Z_{true} [cm];Y_{reco} - Y_{true} [cm]",150,-15,15,150,-15,15);
      h_blip_reszy  ->SetOption("colz");
    h_blip_resx     = dir_truth.make<TH1D>("blip_res_x","Blip position resolution;X_{reco} - X_{true} [cm]",150,-15,15);
    h_qres_electrons= dir_truth.make<TH1D>("qres_electrons","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_qres_alphas   = dir_truth.make<TH1D>("qres_alphas","Collection plane;Cluster charge resolution: ( reco-true ) / true",200,-1.,1.);
    h_adc_factor    = dir_truth.make<TH1D>("adc_per_e","Collection plane;ADC per electron",200,0,0.01);
    h_alpha_qdep    = dir_truth.make<TH1D>("alpha_qdep","True charge deposited by alpha",500,0,10000);
    
    for(int i=0; i<kNplanes; i++) {
      h_nhits[i]      = dir_diag.make<TH1D>(Form("pl%i_nhits",i),  Form("Plane %i;total number of hits",i),hitBins,0,hitMax);
      h_nhits_m[i]    = dir_diag.make<TH1D>(Form("pl%i_nhits_planematched",i), Form("Plane %i;total number of untracked plane-matched hits",i),hitBins,0,hitMax);
      h_hitamp[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_amp",i), Form("Plane %i hits;hit amplitude [ADC]",i),ampBins,0,ampMax);
      h_hitrms[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_rms",i), Form("Plane %i hits;hit RMS [ADC]",i),rmsBins,0,rmsMax);
      h_hitratio[i]   = dir_diag.make<TH1D>(Form("pl%i_hit_ratio",i), Form("Plane %i hits;hit RMS/Amplitude ratio",i),ratioBins,0,ratioMax);
      h_hitint[i]     = dir_diag.make<TH1D>(Form("pl%i_hit_integral",i), Form("Plane %i hits;hit integral [ADC]",i),areaBins,0,areaMax);
      h_hitsigt[i]    = dir_diag.make<TH1D>(Form("pl%i_hit_sigt",i), Form("Plane %i hits;hit time uncertainty [ADC]",i),200,0,5);
      h_nhits_tm[i]   = dir_truth.make<TH1D>(Form("pl%i_nhits_truthmatched",i), Form("Plane %i;number of untracked truth-matched hits",i),hitBins,0,hitMax);
      h_chargecomp[i] = dir_truth.make<TH1D>(Form("pl%i_hit_charge_completeness",i),Form("charge completness, plane %i",i),101,0,1.01);
      h_hitpur[i]     = dir_truth.make<TH1D>(Form("pl%i_hit_purity",i),Form("hit purity, plane %i",i),101,0,1.01);
      h_nelec_TrueVsReco[i] = dir_truth.make<TH2D>( Form("pl%i_nelec_TrueVsReco",i),
        Form("Plane %i;True hit charge [ #times 10^{3} electrons ];reconstructed hit charge [ #times 10^{3} electrons ]",i),60,0,30, 60,0,30);
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
  //fBlipAlg = new blip::BlipRecoAlg( pset.get<fhicl::ParameterSet>("BlipAlg") );
  fhicl::ParameterSet pset_blipalg = pset.get<fhicl::ParameterSet>("BlipAlg");
  fBlipAlg        = new blip::BlipRecoAlg( pset_blipalg );
  fHitProducer    = pset_blipalg.get<std::string>   ("HitProducer");
  fTrkProducer    = pset_blipalg.get<std::string>   ("TrkProducer");
  fGeantProducer  = pset_blipalg.get<std::string>   ("GeantProducer");
  fSimDepProducer = pset_blipalg.get<std::string>   ("SimEDepProducer");
  fCaloPlane      = pset_blipalg.get<int>           ("CaloPlane");
  fDebugMode      = pset.get<bool>                  ("DebugMode",false);
   
  // data tree object
  fData = new BlipAnaTreeDataStruct();
  fData ->treeName        = pset.get<std::string> ("EventTreeName","anatree");
  fData ->blipTreeName    = pset.get<std::string> ("BlipTreeName","bliptree");
  fData ->saveTruthInfo   = pset.get<bool>        ("SaveTruthInfo",true);
  fData ->saveTrkInfo     = pset.get<bool>        ("SaveTrkInfo",true);
  fData ->saveHitInfo     = pset.get<bool>        ("SaveHitInfo",true);
  fData ->saveClustInfo   = pset.get<bool>        ("SaveClustInfo",true);
  fData ->saveEvtTree     = pset.get<bool>        ("SaveEventTree",true);
  fData ->saveBlipTree    = pset.get<bool>        ("SaveBlipTree",true);
  fData ->caloPlane       = fCaloPlane;
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




  //====================================
  // Save MCParticle information
  //====================================
  if( plist.size() ) {
    
    std::vector<blip::ParticleInfo> pinfo = fBlipAlg->pinfo;
    
    fData->total_depEnergy    = 0;
    fData->total_depElectrons = 0;
    fData->total_numElectrons = 0;

    // Loop through the MCParticles
    for(size_t i = 0; i<plist.size(); i++){
      auto pPart = plist[i];

      fData->total_depEnergy     += pinfo[i].depEnergy;
      fData->total_depElectrons  += pinfo[i].depElectrons;
      fData->total_numElectrons  += pinfo[i].numElectrons;

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
        fData->startT[i]          = pinfo[i].time;
        fData->endT[i]            = pinfo[i].endtime;
        fData->pathlen[i]         = pinfo[i].pathLength;
        fData->process[i]         = pPart->Process();
        fData->depEnergy[i]       = pinfo[i].depEnergy;
        fData->depElectrons[i]    = pinfo[i].depElectrons;
        fData->numElectrons[i]    = pinfo[i].numElectrons;
        fData->isPrimary[i]       = pinfo[i].isPrimary;
        if( fDebugMode ) PrintParticleInfo(i);
      }
    } // endloop over G4 particles
    
    std::cout<<"True total energy deposited: "<<fData->total_depEnergy<<" MeV \n";
  
  }//endif particles found in event




  
  //====================================
  // Save TrueBlip information
  //====================================
  std::vector<blip::TrueBlip> trueblips = fBlipAlg->trueblips;
  fData->nedeps = (int)trueblips.size();
  if( trueblips.size() ) {
    std::cout<<"Found "<<trueblips.size()<<" true blips\n";
    for(size_t i=0; i<trueblips.size(); i++ ) {
      fData->edep_tpc[i]    = trueblips.at(i).TPC;
      fData->edep_energy[i] = trueblips.at(i).Energy;
      fData->edep_depne[i]  = trueblips.at(i).DepElectrons;
      fData->edep_charge[i] = trueblips.at(i).NumElectrons;
      fData->edep_ds[i]     = trueblips.at(i).Length;
      fData->edep_x[i]      = trueblips.at(i).Position.X();
      fData->edep_y[i]      = trueblips.at(i).Position.Y();
      fData->edep_z[i]      = trueblips.at(i).Position.Z();
      fData->edep_g4id[i]   = trueblips.at(i).LeadG4ID;
      fData->edep_g4index[i]   = trueblips.at(i).LeadG4Index;
      fData->edep_pdg[i]    = trueblips.at(i).LeadG4PDG;
      if( fDebugMode ) {
        std::cout
        <<"   ~ "<<i<<"  "<<trueblips.at(i).ID<<"  "<<trueblips.at(i).Energy<<" MeV, "
        <<" ds= "<<trueblips.at(i).Length<<" cm, "
        <<" trkID= "<<trueblips.at(i).LeadG4ID<<", pdg "<<trueblips.at(i).LeadG4PDG<<"\n";
      }
    }
  }//endif trueblips were made
  
  std::cout<<"Found "<<hitlist.size()<<" hits from "<<fHitProducer<<"\n";
  std::cout<<"Found "<<tracklist.size()<<" tracks from "<<fTrkProducer<<"\n";
  


  //====================================
  // Save track information
  //====================================
  std::cout<<"Looping over tracks...\n";
  std::map<size_t,std::vector<size_t>> trkhitMap;
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
  }



  //====================================
  // Save hit information
  //====================================
  std::cout<<"Looping over the hits...\n";
  int   num_hits[kNplanes]        ={0};
  int   num_hits_true[kNplanes]   ={0};
  int   num_hits_pmatch[kNplanes] ={0};
  float total_hit_charge[kNplanes]={0};
  int nhits_untracked = 0;
  
  for(size_t i=0; i<hitlist.size(); i++){
   
    int plane = hitlist[i]->WireID().Plane;
    auto const& hinfo = fBlipAlg->hitinfo[i];
    
    fNumHits[plane]++;
    num_hits[plane]++;
    if( hinfo.ismatch ) {
      fNumHitsMatched[plane]++;
      num_hits_pmatch[hinfo.plane]++;
    }
      
    if( hinfo.trkid <= 0 ) {
      fNumHitsUntracked[plane]++;
      nhits_untracked++;
    }
    
    // fill diagnostic histograms
    h_hitamp[plane]    ->Fill(hitlist[i]->PeakAmplitude());
    h_hitsigt[plane]   ->Fill(hitlist[i]->SigmaPeakTime());
    h_hitrms[plane]    ->Fill(hitlist[i]->RMS());
    h_hitratio[plane]  ->Fill(hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
    h_hitint[plane]    ->Fill(hitlist[i]->Integral());
    
    // calculate reco-true resolution
    float qcoll = hinfo.qcoll;
    float qtrue = hinfo.g4charge;
    if( hinfo.g4ids.size() && qcoll>0 && qtrue>0 ) {
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
      fData->hit_time[i]      = hinfo.driftTicks;
      fData->hit_charge[i]    = hinfo.qcoll;
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
    h_nhits[ip]->Fill(num_hits[ip]);
    h_nhits_m[ip]->Fill(num_hits_pmatch[ip]);
    h_nhits_tm[ip]->Fill(num_hits_true[ip]);
    std::cout<<"* plane "<<ip<<": "<<num_hits[ip]<<" hits";
  
    // calculate overall hit purity/completeness per plane
    float qcomp     = -9;
    float pur       = -9;
    float total_ne  = fData->total_numElectrons;
    if( num_hits_true[ip] ) {
      if(total_ne)      qcomp = total_hit_charge[ip]/total_ne;
      if(num_hits[ip])  pur   = num_hits_true[ip]/float(num_hits[ip]);
      h_chargecomp[ip]->Fill( qcomp );
      h_hitpur[ip]    ->Fill( pur );
      std::cout<<"  ("<<num_hits_true[ip]<<" truth-matched) -- completenes: "<<qcomp<<", purity "<<pur;
    }
    std::cout<<"\n";
  }//endloop over planes
    
  if( hitlist.size() ){
    std::cout<<nhits_untracked<<" hits not in 3D tracks \n";
    std::cout<<"("<<100*nhits_untracked/float(hitlist.size())<<"% of hits from "<<fHitProducer<<")\n";
  }
  

  //=============================================
  // Save hit cluster info
  //=============================================
  fData->nclusts = (int)fBlipAlg->hitclust.size();
  for(size_t i=0; i<fBlipAlg->hitclust.size(); i++){
    auto const& clust = fBlipAlg->hitclust[i];
    fData->clust_id[i]        = i;
    fData->clust_tpc[i]       = clust.TPC;
    fData->clust_plane[i]     = clust.Plane;
    fData->clust_wire[i]      = clust.LeadHit->WireID().Wire;
    fData->clust_startwire[i] = clust.StartWire;
    fData->clust_endwire[i]   = clust.EndWire;
    fData->clust_chan[i]      = clust.LeadHit->Channel();
    fData->clust_nwires[i]    = (short)clust.Wires.size();
    fData->clust_nhits[i]     = (short)clust.HitIDs.size();
    fData->clust_lhit_id[i]   = clust.LeadHitID;
    fData->clust_lhit_amp[i]  = (float)clust.LeadHit->PeakAmplitude();
    fData->clust_lhit_rms[i]  = clust.LeadHit->RMS();
    fData->clust_lhit_peakT[i]= clust.LeadHit->PeakTime();
    fData->clust_lhit_gof[i]  = clust.LeadHit->GoodnessOfFit();
    fData->clust_lhit_isfit[i]= (clust.LeadHit->GoodnessOfFit()>=0);
    fData->clust_lhit_time[i] = clust.LeadHitTime;
    fData->clust_sumadc[i]    = clust.ADCs;
    fData->clust_charge[i]    = clust.Charge;
    fData->clust_time[i]      = clust.Time;
    fData->clust_time_err[i]  = clust.TimeErr;
    fData->clust_startTime[i] = clust.StartTime;
    fData->clust_endTime[i]   = clust.EndTime;
    fData->clust_timespan[i]  = clust.EndTime-clust.StartTime;
    fData->clust_ismatch[i]   = clust.isMatched;
    fData->clust_blipid[i]    = clust.BlipID;

      
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

    }//endloop over trueblips
    
    if( fDebugMode && clust.Plane == 2 ) PrintClusterInfo(clust);
  
  }//endloop over 2D hit clusters



  //====================================
  // Save blip info to tree
  //===================================
  fData->nblips             = fBlipAlg->blips.size();
  int nblips_matched        = 0;
  int nblips_total          = 0;
  float true_blip_charge    = 0;
  fData->total_blip_energy  = 0;
  for(size_t i=0; i<fBlipAlg->blips.size(); i++){
    auto& b = fBlipAlg->blips[i];
    
    fNum3DBlips++;
    if( b.NPlanes == 3 ) fNum3DBlips3Plane++;

    nblips_total++;
    fData->blip_id[i]         = i;
    fData->blip_tpc[i]        = b.TPC;
    fData->blip_nplanes[i]    = b.NPlanes;
    fData->blip_maxdiff[i]    = b.MaxIntersectDiff;
    fData->blip_x[i]          = b.X;
    fData->blip_y[i]          = b.Y;
    fData->blip_z[i]          = b.Z;
    fData->blip_trkdist[i]    = b.TrkDist;
    fData->blip_trkid[i]      = b.TrkID;
    fData->blip_incylinder[i] = b.inCylinder;
   
    for(size_t ipl = 0; ipl<kNplanes; ipl++){
      fData->blip_clustid[ipl][i]     = b.ClustID[ipl];
      //fData->blip_matchscore[ipl][i]  = b.MScore[ipl];
      fData->blip_nhits[ipl][i]       = b.NHits[ipl];
      fData->blip_nwires[ipl][i]      = b.NWires[ipl];
      fData->blip_timespan[ipl][i]    = b.Timespan[ipl];
    }
   
    if( b.Charge[fCaloPlane] > 0 ) {
      fData->blip_charge[i]     = b.Charge[fCaloPlane];
      fData->blip_sumadc[i]     = b.SumADC[fCaloPlane];
      fData->blip_energy[i]     = b.Energy;
      fData->total_blip_energy += b.Energy;
    }
    
    h_blip_zy->Fill(b.Z, b.Y);
    h_blip_charge->Fill(b.Charge[fCaloPlane]);
    h_blip_sumadc->Fill(b.SumADC[fCaloPlane]);
    h_blip_nplanes->Fill(b.NPlanes);
    if( b.NPlanes == kNplanes ) {
      h_blip_charge_3p->Fill(b.Charge[fCaloPlane]);
      h_blip_zy_3p->Fill(b.Z, b.Y);
    }
   
    // -----------------------------------------------
    // save the clustIDs and true energy deposits to the blip
    // (use the association between clust <--> edep)
    // -----------------------------------------------
    float max = 0;
    for(auto clustID : b.ClustID_set ) {
      int edepid      = fData->clust_edepid[clustID];
      if( edepid >= 0 && edepid < fData->nedeps ) {
        float E = trueblips[edepid].Energy;
        if( E > max ) {
          fData->blip_edepid[i] = edepid;
          fData->edep_blipid[edepid] = i;
          max = E;
        }
      }
    }
    
    // -----------------------------------------------
    // if a true energy dep was found, fill some histograms
    // and iterate some counters
    // -----------------------------------------------
    int b_eid = fData->blip_edepid[i];
    if( b_eid >= 0 ) {
      fNum3DBlipsTrue++;
      nblips_matched++;
      true_blip_charge+= trueblips[b_eid].NumElectrons;
      TVector3 truePos = trueblips[b_eid].Position;
      h_blip_reszy->Fill( b.Position.Z()-truePos.Z(), b.Position.Y()-truePos.Y() );
      h_blip_resx->Fill( b.Position.X()-truePos.X() );
    }
  
    // -----------------------------------
    // add this blip to the blip tree
    //------------------------------------
    fData->FillBlipTree(b);
  
  }//endloop over 3D blips
  
  // Fill some more histograms...
  h_nblips->Fill(nblips_total);
  if( fIsMC ) {
    h_nblips_tm->Fill(nblips_matched);
    if( fData->total_numElectrons ) h_blip_qcomp->Fill(true_blip_charge / fData->total_numElectrons);
    if( nblips_total ) h_blip_pur    ->Fill( float(nblips_matched) / float(nblips_total) );
  }

  std::cout<<"Reconstructed "<<fBlipAlg->blips.size()<<" 3D blips";
  if( nblips_matched ) std::cout<<" ("<<nblips_matched<<" were truth-matched)";
  std::cout<<"\n";

  if( fDebugMode ) {
    for(auto const& b : fBlipAlg->blips ) {
      std::cout
      <<"   -- "<<b.ID<<", TPC: "<<b.TPC
      <<"; charge: "<<b.Charge[2]
      <<"; recoEnergy: "<<b.Energy<<" MeV"
      <<"; Position: "<<b.Position.X()<<", "<<b.Position.Y()<<", "<<b.Position.Z()
      <<"; MaxIntersectDiff: "<<b.MaxIntersectDiff
      <<"; EdepID: "<<fData->blip_edepid[b.ID]
      <<"\n";
    }
  }
  
  
  
  //====================================
  // Fill TTree
  //====================================
  fData->FillEventTree();

}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipAna::endJob(){
  
  printf("\n***********************************************\n");
  fBlipAlg->PrintConfig();
  printf("BlipAna Summary\n\n");
  printf("  Total events              : %i\n",            fNumEvents);
  printf("  Blips per evt             : %.3f\n", fNum3DBlips/float(fNumEvents));
  printf("  Matched planes per blip   : %.3f\n", h_blip_nplanes->GetMean());
  if(fIsMC){
  printf("  MC-matched blips per evt  : %.3f\n",       fNum3DBlipsTrue/float(fNumEvents));
  printf("  Blip charge completeness  : %.4f +/- %.4f\n",       h_blip_qcomp->GetMean(), h_blip_qcomp->GetStdDev()/sqrt(fNumEvents));
  printf("  Blip purity               : %.4f\n",       h_blip_pur->GetMean());
  }
  printf("  Mean blip charge [e-]     : %.0f\n",      h_blip_charge->GetMean());
  printf("\n");
  for(size_t i=0; i<kNplanes; i++){
  printf("  Plane %lu -------------------------\n",i);
  printf("   * total hits/evt        : %.2f\n",fNumHits[i]/(float)fNumEvents);
  printf("   * untracked hits/evt    : %.2f\n",fNumHitsUntracked[i]/(float)fNumEvents);
  printf("   * plane-matched hits/evt: %.2f\n",fNumHitsMatched[i]/(float)fNumEvents);
  if(fIsMC) {
  printf("   * true-matched hits/evt : %.2f\n",fNumHitsTrue[i]/(float)fNumEvents);
  printf("   * charge completeness   : %.4f\n",h_chargecomp[i]->GetMean());
  printf("   * hit purity            : %.4f\n",h_hitpur[i]->GetMean());
  }
  } 
  printf("\n***********************************************\n");
 
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

void BlipAna::PrintHitInfo(const blip::HitInfo& hi){
  printf("  hitID: %4i, TPC: %i, plane: %i, driftTicks: %7.2f, leadWire: %3i, G4ID: %4i, recoTrack: %4i\n",
    hi.hitid,
    hi.tpc,
    hi.plane,
    hi.driftTicks,
    hi.wire,
    hi.g4id,
    hi.trkid
  );
}

void BlipAna::PrintClusterInfo(const blip::HitClust& hc){
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
