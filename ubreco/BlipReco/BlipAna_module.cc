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
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
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
#include "TTimeStamp.h"


// Set global constants and max array sizes
const int kMaxHits  = 100000;
const int kMaxTrks  = 10000;
const int kMaxShwrs = 10000;
const int kMaxBlips = 5000;
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
  std::string treeName            = "anatree";
  bool  saveTruthInfo             = true;
  bool  saveHitInfo               = true;

  // --- Event information ---   
  int   event;                    // event number
  int   run;                      // run number
  float lifetime;                 // electron lifetime
  double timestamp;                // unix time of event
  //double timestamp_aux;                // unix time of event
  
  // --- G4 information ---
  float total_depEnergy;          // total deposited energy in AV
  float total_numElectrons;       // total electrons reaching anode wires
  float gamma_depEnergy;          // total gamma-induced energy deposited
  float gamma_numElectrons;       // total electrons from gamma-induced depositions
  int   nparticles;               // number of G4 particles
  int   isPrimary[kMaxG4];        // is primary particle
  int   trackID[kMaxG4];          // G4 track ID
  int   pdg[kMaxG4];              // PDG
  int   nDaughters[kMaxG4];       // number of daughters
  int   mother[kMaxG4];           // mother particle
  float E[kMaxG4];                // initial energy (MeV)
  float endE[kMaxG4];             // final energy (MeV)
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
  int   numElectrons[kMaxG4];     // electrons reaching anode wires
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
  float edep_energyESTAR[kMaxEDeps];   // total energy deposited (MeV)
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
  float	hit_ph[kMaxHits];         // amplitude
  float	hit_area[kMaxHits];       // charge (area) in ADC units
  float	hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int	  hit_shwrid[kMaxHits];      // is this hit associated with a reco shower?
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

  // --- Hit cluster information ---
  int   nclusts;
  int   clust_tpc[kMaxHits];
  int   clust_plane[kMaxHits];
  int   clust_wire[kMaxHits];
  int   clust_nwires[kMaxHits];
  int   clust_nhits[kMaxHits];
  int   clust_lhit_id[kMaxHits];
  float   clust_lhit_ph[kMaxHits];
  float   clust_lhit_rms[kMaxHits];
  float   clust_lhit_time[kMaxHits];
  float clust_lhit_peakT[kMaxHits];
  float clust_lhit_gof[kMaxHits];
  float clust_charge[kMaxHits];
  float clust_time[kMaxHits];
  float clust_time_w[kMaxHits];
  float clust_time_err[kMaxHits];
  //float clust_time_lh[kMaxHits];
  float clust_startTime[kMaxHits];
  float clust_endTime[kMaxHits];
  float clust_timespan[kMaxHits];
  //float clust_xpos[kMaxHits];
  float clust_g4energy[kMaxHits];
  float clust_g4charge[kMaxHits];
  int   clust_g4id[kMaxHits];
  int   clust_ismatch[kMaxHits];
  int   clust_blipid[kMaxHits];
  int   clust_edepid[kMaxHits];
  
  // --- Blip information ---
  float total_blip_energy;
  int   nblips;
  int   blip_tpc[kMaxBlips];
  int   blip_nplanes[kMaxBlips];
  //int   blip_caloplane[kMaxBlips];
  float blip_x[kMaxBlips];
  float blip_y[kMaxBlips];
  float blip_z[kMaxBlips];
  float blip_maxdiff[kMaxBlips];
  //float blip_charge[kMaxBlips];
  float blip_energy[kMaxBlips];
  float blip_energyESTAR[kMaxBlips];
  int   blip_edepid[kMaxBlips];  
  //int   blip_clustid[kNplanes][kMaxBlips];
  //float blip_charge[kNplanes][kMaxBlips];
  int   blip_clustid[kNplanes][kMaxBlips];
  float blip_charge[kNplanes][kMaxBlips];

  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    lifetime              = -999;
    timestamp             = -999;
    //timestamp_aux             = -999;
    total_depEnergy       = -999;
    total_numElectrons    = -999;
    gamma_depEnergy       = -999;
    gamma_numElectrons    = -999;
    nparticles            = 0;    // --- G4 particles ---
    FillWith(isPrimary,   -9);
    FillWith(trackID,     -999);
    FillWith(pdg,         -99999);
    FillWith(nDaughters,  -999);
    FillWith(mother,      -999);
    FillWith(E,           -999.);
    FillWith(endE,        -999.);
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
    FillWith(edep_energyESTAR, -999);
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
    FillWith(hit_ph,      -999);
    FillWith(hit_area,    -999);
    FillWith(hit_sumadc,  -999);
    FillWith(hit_mult,    -999);
    FillWith(hit_charge,  -999);
    FillWith(hit_isreal,  -9);
    FillWith(hit_ismatch, -9);
    FillWith(hit_isgood, -9);
    FillWith(hit_trkid,   -9);
    FillWith(hit_shwrid,  -9);
    FillWith(hit_g4id,    -999);
    FillWith(hit_g4frac,  -9);
    FillWith(hit_g4energy,-999);
    FillWith(hit_g4charge,-999);
    FillWith(hit_clustid, -9);
    FillWith(hit_blipid,  -9);
    FillWith(hit_gof,  -9);
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
    nshwrs                = 0;    // --- Showers --- 
    FillWith(shwr_id,     -999);
    FillWith(shwr_dirx,   -999);
    FillWith(shwr_diry,   -999);
    FillWith(shwr_dirz,   -999);
    FillWith(shwr_startx, -999);
    FillWith(shwr_starty, -999);
    FillWith(shwr_startz, -999);
    FillWith(shwr_length, -999);
    FillWith(shwr_openangle, -999);
    nclusts               = 0;    // --- Hit Clusters ---
    FillWith(clust_tpc,  -9);
    FillWith(clust_plane, -9);
    FillWith(clust_wire, -9);
    FillWith(clust_nwires, -9);
    FillWith(clust_nhits, -9);
    FillWith(clust_lhit_id, -9);
    FillWith(clust_lhit_ph, -9);
    FillWith(clust_lhit_rms, -9);
    FillWith(clust_lhit_time, -9);
    FillWith(clust_lhit_peakT, -9);
    FillWith(clust_lhit_gof, -9);
    FillWith(clust_charge, -999);
    FillWith(clust_time, -999);
    //FillWith(clust_time_lh, -999);
    FillWith(clust_time_w, -999);
    FillWith(clust_time_err, -999);
    FillWith(clust_startTime, -999);
    FillWith(clust_endTime, -999);
    FillWith(clust_timespan, -9);
    //FillWith(clust_xpos, -999);
    FillWith(clust_g4id, -9);
    FillWith(clust_g4charge, -999);
    FillWith(clust_g4energy, -999);
    FillWith(clust_ismatch, -9);
    FillWith(clust_edepid, -9);
    FillWith(clust_blipid, -9);
    total_blip_energy     = -9;  // --- Blips ---
    nblips                = 0;
    FillWith(blip_tpc,  -9);
    FillWith(blip_nplanes,  -9);
    //FillWith(blip_caloplane, -9);
    FillWith(blip_x, -99999);
    FillWith(blip_y, -99999);
    FillWith(blip_z, -99999);
    FillWith(blip_maxdiff, -9);
    //FillWith(blip_charge, -999);
    FillWith(blip_energy, -999);
    FillWith(blip_energyESTAR, -999);
    FillWith(blip_edepid, -9);
    for(int i=0; i<kNplanes; i++){
      FillWith(blip_clustid[i], -9);
      FillWith(blip_charge[i], -999);
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
    tree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    tree->Branch("event",&event,"event/I");
    tree->Branch("run",&run,"run/I");
    tree->Branch("timestamp",&timestamp,"timestamp/D");
    //tree->Branch("timestamp_aux",&timestamp_aux,"timestamp_aux/D");
    tree->Branch("lifetime",&lifetime,"lifetime/F");
  
    if( saveTruthInfo ) {
    tree->Branch("total_depEnergy",&total_depEnergy,"total_depEnergy/F");
    tree->Branch("total_numElectrons",&total_numElectrons,"total_numElectrons/F");
    tree->Branch("nparticles",&nparticles,"nparticles/I");
    tree->Branch("isPrimary",isPrimary,"isPrimary[nparticles]/I");
    tree->Branch("trackID",trackID,"trackID[nparticles]/I");
    tree->Branch("pdg",pdg,"pdg[nparticles]/I");
    tree->Branch("nDaughters",nDaughters,"nDaughters[nparticles]/I");
    tree->Branch("mother",mother,"mother[nparticles]/I");
    tree->Branch("E",E,"E[nparticles]/F");
    tree->Branch("endE",endE,"endE[nparticles]/F");
    tree->Branch("mass",mass,"mass[nparticles]/F");
    tree->Branch("P",P,"P[nparticles]/F");
    tree->Branch("Px",Px,"Px[nparticles]/F");
    tree->Branch("Py",Py,"Py[nparticles]/F");
    tree->Branch("Pz",Pz,"Pz[nparticles]/F");
    tree->Branch("startPointx",startPointx,"startPointx[nparticles]/F");
    tree->Branch("startPointy",startPointy,"startPointy[nparticles]/F");
    tree->Branch("startPointz",startPointz,"startPointz[nparticles]/F");
    tree->Branch("endPointx",endPointx,"endPointx[nparticles]/F");
    tree->Branch("endPointy",endPointy,"endPointy[nparticles]/F");
    tree->Branch("endPointz",endPointz,"endPointz[nparticles]/F");
    tree->Branch("startT",startT,"startT[nparticles]/F");
    tree->Branch("endT",endT,"endT[nparticles]/F");
    tree->Branch("pathlen",pathlen,"pathlen[nparticles]/F");
    tree->Branch("numElectrons",numElectrons,"numElectrons[nparticles]/I");
    tree->Branch("depEnergy",depEnergy,"depEnergy[nparticles]/F");
    tree->Branch("process",&process);
    tree->Branch("nedeps",&nedeps,"nedeps/I");
//    tree->Branch("edep_tpc",edep_tpc,"edep_tpc[nedeps]/I"); 
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
   
    if( saveHitInfo ) {
    tree->Branch("nhits",&nhits,"nhits/I");
    tree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I"); 
    tree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I"); 
    tree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I"); 
    //tree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I"); 
    tree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F"); 
    tree->Branch("hit_time",hit_time,"hit_time[nhits]/F"); 
    tree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F"); 
    tree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/F"); 
    tree->Branch("hit_area",hit_area,"hit_area[nhits]/F"); 
    tree->Branch("hit_sumadc",hit_sumadc,"hit_sumadc[nhits]/F"); 
    tree->Branch("hit_mult",hit_mult,"hit_mult[nhits]/I"); 
    tree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
    tree->Branch("hit_ismatch",hit_ismatch,"hit_ismatch[nhits]/I");
    tree->Branch("hit_isgood",hit_isgood,"hit_isgood[nhits]/I");
    tree->Branch("hit_isreal",hit_isreal,"hit_isreal[nhits]/I"); 
    tree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
    tree->Branch("hit_shwrid",hit_shwrid,"hit_shwrid[nhits]/I"); 
    if( saveTruthInfo ) {
    tree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
    tree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F"); 
    tree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F"); 
    tree->Branch("hit_g4charge",hit_g4charge,"hit_g4charge[nhits]/F"); 
    }
    tree->Branch("hit_clustid",hit_clustid,"hit_clustid[nhits]/I"); 
    tree->Branch("hit_blipid",hit_blipid,"hit_blipid[nhits]/I");
    tree->Branch("hit_gof",hit_gof,"hit_gof[nhits]/F");
    }
    tree->Branch("ntrks",&ntrks,"ntrks/I");          
    tree->Branch("trk_id",trk_id,"trk_id[ntrks]/I");       
    tree->Branch("trk_npts",trk_npts,"trk_npts[ntrks]/I");
    tree->Branch("trk_length",trk_length,"trk_length[ntrks]/F");
    tree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
    tree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
    tree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
    tree->Branch("trk_startd",trk_startd,"trk_startd[ntrks]/F");
    tree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
    tree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
    tree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
    tree->Branch("trk_endd",trk_endd,"trk_endd[ntrks]/F");
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
    tree->Branch("nclusts",&nclusts,"nclusts/I");
    tree->Branch("clust_tpc",clust_tpc,"clust_tpc[nclusts]/I");
    tree->Branch("clust_plane",clust_plane,"clust_plane[nclusts]/I");
    tree->Branch("clust_wire",clust_wire,"clust_wire[nclusts]/I");
    tree->Branch("clust_nwires",clust_nwires,"clust_nwires[nclusts]/I");
    tree->Branch("clust_nhits",clust_nhits,"clust_nhits[nclusts]/I");
    tree->Branch("clust_lhit_id",clust_lhit_id,"clust_lhit_id[nclusts]/I");
    tree->Branch("clust_lhit_ph",clust_lhit_ph,"clust_lhit_ph[nclusts]/F");
    tree->Branch("clust_lhit_rms",clust_lhit_rms,"clust_lhit_rms[nclusts]/F");
    //tree->Branch("clust_lhit_time",clust_lhit_time,"clust_lhit_time[nclusts]/F");
    tree->Branch("clust_lhit_peakT",clust_lhit_peakT,"clust_lhit_peakT[nclusts]/F");
    tree->Branch("clust_lhit_gof",clust_lhit_gof,"clust_lhit_gof[nclusts]/F");
    tree->Branch("clust_charge",clust_charge,"clust_charge[nclusts]/F");
    tree->Branch("clust_time",clust_time,"clust_time[nclusts]/F");
    tree->Branch("clust_time_w",clust_time_w,"clust_time_w[nclusts]/F");
    tree->Branch("clust_time_err",clust_time_err,"clust_time_err[nclusts]/F");
    //tree->Branch("clust_startTime",clust_startTime,"clust_startTime[nclusts]/F");
    //tree->Branch("clust_endTime",clust_endTime,"clust_endTime[nclusts]/F");
    tree->Branch("clust_timespan",clust_timespan,"clust_timespan[nclusts]/F");
    //tree->Branch("clust_xpos",clust_xpos,"clust_xpos[nclusts]/F");
    tree->Branch("clust_g4charge",clust_g4charge,"clust_g4charge[nclusts]/F");
    tree->Branch("clust_g4energy",clust_g4energy,"clust_g4energy[nclusts]/F");
    tree->Branch("clust_ismatch",clust_ismatch,"clust_ismatch[nclusts]/I");
    tree->Branch("clust_edepid",clust_edepid,"clust_edepid[nclusts]/I");
    tree->Branch("clust_blipid",clust_blipid,"clust_blipid[nclusts]/I");
    tree->Branch("total_blip_energy",&total_blip_energy,"total_blip_energy/F");
    tree->Branch("nblips",&nblips,"nblips/I");
    tree->Branch("blip_tpc",blip_tpc,"blip_tpc[nblips]/I");
    tree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/I");
    //tree->Branch("blip_caloplane",blip_caloplane,"blip_caloplane[nblips]/I");
    tree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
    tree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
    tree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
    tree->Branch("blip_maxdiff",blip_maxdiff,"blip_maxdiff[nblips]/F");
    //tree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
    tree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
    tree->Branch("blip_energyESTAR",blip_energyESTAR,"blip_energyESTAR[nblips]/F");
    tree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
    for(int i=0; i<kNplanes; i++) 
      tree->Branch(Form("blip_clustid_pl%i",i),blip_clustid[i],Form("blip_clustid_pl%i[nblips]/I",i));
    for(int i=0; i<kNplanes; i++)
      tree->Branch(Form("blip_charge_pl%i",i),blip_charge[i],Form("blip_charge_pl%i[nblips]/F",i));
    
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

  // --- Detector and clock data ---
  //float TickOffset[kNplanes];

  // --- Data and calo objects ---
  BlipAnaTreeDataStruct*  fData;
  calo::CalorimetryAlg    fCaloAlg;
  TGraph2D*               ESTAR;

  // --- FCL configs ---
  std::string         fAnaTreeName;
  std::string         fHitModuleLabel;
  std::string         fLArG4ModuleLabel;
  std::string         fTrackModuleLabel;
  std::string         fShowerModuleLabel;
  std::string         fCaloModuleLabel;
  bool                fSaveTruthInfo;
  bool                fSaveHitInfo;
  bool                fMakeDiagHists;
//  bool                fMakeTruthHists;
  float               fTrueBlipMergeDist;
  bool                fDoHitFiltering;
  std::vector<float>  fMinHitRMS;
  std::vector<float>  fMaxHitRMS;
  std::vector<float>  fMinHitRatio;
  std::vector<float>  fMaxHitRatio;
  std::vector<float>  fMinHitGOF;
  //std::vector<float>  fMaxHitGOF;
  float               fHitClustWidthFact;
  int                 fHitClustWireRange;
  float               fHitMatchWidthFact;
  float               fHitMatchMaxTicks;
  int                 fMaxWiresInCluster;
  int                 fCaloPlane;         // use this plane for calorimetry
  int                 fMinMatchedPlanes;  
  float               fDiffTolerance;     // 1cm in MicroBooNE

  // --- Counters and such ---
  bool  fIsRealData         = 0;
  int   fNumEvents          = 0;
  int   fNumHits_pl[3][2]   = {};
  int   fNum3DBlips         = 0;

  // --- Histograms ---
  //TH1D*   h_nhits[kNplanes][2];
  //TH1D*   h_hitph[kNplanes][2];
  //TH1D*   h_hitrms[kNplanes][2];
  //TH1D*   h_hitratio[kNplanes][2];
  //TH1D*   h_hitint[kNplanes][2];
//  TH1D*   h_nhits[kNplanes];
  TH1D*   h_hitph[kNplanes];
  TH1D*   h_hitrms[kNplanes];
  TH1D*   h_hitratio[kNplanes];
  TH1D*   h_hitint[kNplanes];
  //TH2D*   h_hitrms_vs_ph[kNplanes][2];
  //TH1D*   h_sumadc_res[kNplanes][2];
  //TH2D*   h_hitarea_vs_sumadc[kNplanes][2];
  TH2D*   h_nelec_TrueVsReco[kNplanes];
  TH1D*   h_nelec_Resolution[kNplanes];
  TH1D*   h_chargecomp[kNplanes];
  TH1D*   h_hitpur[kNplanes];
  TH1D*   h_hit_dt;
  TH1D*   h_hit_dtfrac;
  TH1D*   h_nmatches[3];
  TH1D*   h_hitCuts[3][3];
  //TH1D*   h_trueBlipEw;
  TH1D*   h_hitgof_2D[3];
  TH1D*   h_hitgof_3D[3];
  TH1D*   h_clust_nwires;

  TH2D*   h_blip_zy;

  //TH1D*   h_energyDep_electrons;
  //TH1D*   h_energyDep_gammas;

  // Initialize histograms
  void InitializeHistograms(){
    art::ServiceHandle<art::TFileService> tfs;
  
    h_hit_dt      = tfs->make<TH1D>("hit_dt","Hit matching;dT [ticks]",200,0,20);
    h_hit_dtfrac  = tfs->make<TH1D>("hit_dtfrac","Hit matching;dT/RMS",200,0,5);
   
    h_blip_zy   = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy   ->SetOption("COLZ");

      //h_energyDep_electrons = tfs->make<TH1D>("energyDep_electrons","True energy deps from electrons (PDG 11);MeV",100,0,1.);
      //h_energyDep_gammas    = tfs->make<TH1D>("energyDep_gammas","True energy deps from gammas (PDG 22);MeV",100,0,1.);
      
    //h_trueBlipEw = tfs->make<TH1D>("trueBlipE_weighted","True blip energy (weighted);Energy [MeV];#sum E",500,0,5.);
    //h_trueBlipEw ->SetOption("HIST");


    if( fMakeDiagHists ) {

      art::TFileDirectory diagDir = tfs->mkdir("Diagnostics");
      //float hitMax  = 300;    int hitBins   = 300;
      float phMax = 50;       int phBins    = 250;
      float rmsMax = 15;      int rmsBins   = 300;
      float areaMax = 300;    int areaBins  = 300;
      float ratioMax = 5.0;   int ratioBins = 250;
        
      h_clust_nwires = diagDir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);

      for(int i=0; i<kNplanes; i++) {
          
        //h_nhits[i]    = tfs->make<TH1D>(Form("pl%i_nhits",i),  Form("Plane %i;Number of hits",i),hitBins,0,hitMax);
        h_hitph[i]    = diagDir.make<TH1D>(Form("pl%i_hit_ph",i), Form("Plane %i hits;Pulse height [ADC]",i),phBins,0,phMax);
        h_hitrms[i]   = diagDir.make<TH1D>(Form("pl%i_hit_rms",i), Form("Plane %i hits;Hit RMS [ADC]",i),rmsBins,0,rmsMax);
        h_hitratio[i] = diagDir.make<TH1D>(Form("pl%i_hit_ratio",i), Form("Plane %i hits;RMS/amp ratio [ADC]",i),ratioBins,0,ratioMax);
        h_hitint[i]   = diagDir.make<TH1D>(Form("pl%i_hit_integral",i), Form("Plane %i hits;Hit integral [ADC]",i),areaBins,0,areaMax);
        h_hitgof_2D[i]= diagDir.make<TH1D>(Form("pl%i_hit_gof_2d",i), Form("Plane %i hits (NOT plane-matched);Goodness of fit [ADC]",i),500,0,10);
        h_hitgof_3D[i]= diagDir.make<TH1D>(Form("pl%i_hit_gof_3d",i), Form("Plane %i hits (plane-matched);Goodness of fit [ADC]",i),500,0,10);

        /*
        for(int j=0;j<2;j++) {
          std::vector<std::string> label = { "fake", "real" };
          h_nhits[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_nhits",i,label[j].c_str()),  Form("Plane %i, %s hits;Number of hits",i,label[j].c_str()),hitBins,0,hitMax);
          h_hitph[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_hit_ph",i,label[j].c_str()), Form("Plane %i, %s hits;Pulse height [ADC]",i,label[j].c_str()),phBins,0,phMax);
          h_hitrms[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_hit_rms",i,label[j].c_str()),Form("Plane %i, %s hits;Hit RMS [ticks]",i,label[j].c_str()),rmsBins,0,rmsMax);
          //h_hitrms_vs_ph[i][j]  = diagDir.make<TH2D>(Form("pl%i_%s_hit_rms_vs_ph",i,label[j].c_str()),
           // Form("Plane %i, %s hits;Hit RMS [ticks];Pulse height [ADC]",i,label[j].c_str()),
           // rmsBins/2,0,rmsMax,
           // phBins/2,0,phMax);
          //h_hitarea_vs_sumadc[i][j]  = diagDir.make<TH2D>(Form("pl%i_%s_hit_area_vs_sumadc",i,label[j].c_str()),
          //  Form("Plane %i, %s hits;Hit area [ticks*ADC];Summed ADC [ticks*ADC]",i,label[j].c_str()),
          //  areaBins,0,areaMax,
          //  areaBins,0,areaMax);
          //h_sumadc_res[i][j] = diagDir.make<TH1D>(Form("pl%i_%s_sumadc_res",i,label[j].c_str()),
          //  Form("Plane %i, %s hits;(SummedADC - Integral) / Integral",i,label[j].c_str()),
          //  300,-3,3);
          h_hitratio[i][j]  = diagDir.make<TH1D>(Form("pl%i_%s_hit_ratio",i,label[j].c_str()),
            Form("Plane %i, %s hits;RMS/amp Ratio",i,label[j].c_str()),ratioBins, 0, ratioMax);
          h_hitint[i][j]  = diagDir.make<TH1D>(Form("pl%i_%s_hit_integral",i,label[j].c_str()),
            Form("Plane %i, %s hits;Hit integral [ADC]",i,label[j].c_str()),areaBins,0,areaMax);
        }
        */
       
          h_nelec_TrueVsReco[i] = diagDir.make<TH2D>( Form("pl%i_nelec_TrueVsReco",i),
            Form("Plane %i;True hit charge [ #times 10^{3} electrons ];Reconstructed hit charge [ #times 10^{3} electrons ]",i),60,0,30, 60,0,30);
          h_nelec_TrueVsReco[i] ->SetOption("colz");
          h_nelec_Resolution[i] = diagDir.make<TH1D>( Form("pl%i_nelec_res",i),Form("Plane %i;Hit charge resolution: (reco-true)/true",i),200,-2,2);
          h_chargecomp[i] = diagDir.make<TH1D>(Form("pl%i_charge_completeness",i),Form("Charge completness, plane %i",i),101,0,1.01);
          h_hitpur[i]     = diagDir.make<TH1D>(Form("pl%i_hit_purity",i),Form("Hit purity, plane %i",i),101,0,1.01);
          
        h_nmatches[i] = diagDir.make<TH1D>(Form("pl%i_nmatches",i),Form("Number of matched hits, plane %i",i),10,0,10);


        for(int k=0; k<3; k++){
          std::vector<std::string> label = {"all","real","fake"};
          h_hitCuts[i][k]=diagDir.make<TH1D>(Form("pl%i_hitCuts_%s",i,label[k].c_str()),Form("Hits per event, plane %i, %s",i,label[k].c_str()),6,0,6);
          h_hitCuts[i][k]->SetBit(TH1::kIsAverage);
          h_hitCuts[i][k]->SetOption("HIST TEXT");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(1,"total");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(2,"untracked");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(3,"rms");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(4,"ratio");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(5,"matched");
          h_hitCuts[i][k]->GetXaxis()->SetBinLabel(6,"blip");
        }

      }//endloop over planes
      
    
    }
 
 }

};//class BlipAna


//###################################################
//  BlipAna constructor and destructor
//###################################################
BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData              (nullptr)
  ,fCaloAlg           (pset.get< fhicl::ParameterSet >    ("CaloAlg"))
  ,fAnaTreeName       (pset.get< std::string >            ("AnaTreeName",       "anatree"))
  ,fHitModuleLabel    (pset.get< std::string >            ("HitModuleLabel",    "gaushit"))
  ,fLArG4ModuleLabel  (pset.get< std::string >            ("LArG4ModuleLabel",  "largeant"))
  ,fTrackModuleLabel  (pset.get< std::string >            ("TrackModuleLabel",  "pandoraTrack"))
  ,fShowerModuleLabel (pset.get< std::string >            ("ShowerModuleLabel", "pandoraShower"))
  ,fCaloModuleLabel   (pset.get< std::string >            ("CaloModuleLabel",   "pandoraCalo"))
  ,fSaveTruthInfo  (pset.get< bool >                      ("SaveTruthInfo",     true))
  ,fSaveHitInfo  (pset.get< bool >                        ("SaveHitInfo",     true))
  ,fMakeDiagHists    (pset.get< bool >                    ("MakeDiagnosticHists",true))
  ,fTrueBlipMergeDist     (pset.get< float >                  ("TrueBlipMergeDist",     0.3))
  ,fDoHitFiltering    (pset.get< bool >                   ("DoHitFiltering",    false))
  ,fMinHitRMS         (pset.get< std::vector< float > >   ("MinHitRMS",         {-9999,-9999,-9999}))
  ,fMaxHitRMS         (pset.get< std::vector< float > >   ("MaxHitRMS",         { 9999, 9999, 9999}))
  ,fMinHitRatio       (pset.get< std::vector< float > >   ("MinHitRatio",       {-9999,-9999,-9999}))
  ,fMaxHitRatio       (pset.get< std::vector< float > >   ("MaxHitRatio",       { 9999, 9999, 9999}))
  ,fMinHitGOF         (pset.get< std::vector< float > >   ("MinHitGOF",         {-9999,-9999,-9999}))
  ,fHitClustWidthFact (pset.get< float >                  ("HitClustWidthFact", 2.5))
  ,fHitClustWireRange (pset.get< int >                    ("HitClustWireRange", 1))
  ,fHitMatchWidthFact (pset.get< float >                  ("HitMatchWidthFact", 1.))
  ,fHitMatchMaxTicks  (pset.get< float >                  ("HitMatchMaxTicks",  5.))
  ,fMaxWiresInCluster (pset.get< int >                    ("MaxWiresInCluster", 5 ))
  ,fCaloPlane         (pset.get< int >                    ("CaloPlane",         2))
//  ,fMinMatchedPlanes  (pset.get< int >                    ("MinMatchedPlanes",  3))
//  ,fDiffTolerance     (pset.get< float >                  ("DiffTolerance",     2.0))
{
  // TODO: move above into "reconfigure" function; parameters to be set in blipana_config.fcl
  fData = new BlipAnaTreeDataStruct();
  fData ->saveTruthInfo = fSaveTruthInfo;
  fData ->saveHitInfo = fSaveHitInfo;
  fData ->treeName = fAnaTreeName;
  fData ->Clear();
  fData ->MakeTree();
  InitializeHistograms();
  fDiffTolerance = pset.get< float >("DiffTolerance", 999.);
  fMinMatchedPlanes = pset.get< int >("MinMatchedPlanes", 2);
}
BlipAna::~BlipAna(){}



//###################################################
//  beginJob: retrieve relevant detector and clock
//  data here, and save them into class variables
//###################################################
void BlipAna::beginJob() {
  // -- Detector Properties --
  //auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  //std::cout<<"Electron lifetime = "<<detProp->ElectronLifetime()<<"\n";
  
  BlipUtils::InitializeUtils();

  art::ServiceHandle<geo::Geometry> geom;
  std::cout<<"Checking geometrty\n";
  std::cout<<"  DetHalfWidth "<<geom->DetHalfWidth()<<"\n";
  std::cout<<"  DetHalfHeight "<<geom->DetHalfHeight()<<"\n";
  std::cout<<"  DetLength "<<geom->DetLength()<<"\n";
  
  // Load the ESTAR lookup table
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

}


//###################################################
//  Main event-by-event analysis
//###################################################
void BlipAna::analyze(const art::Event& evt)
{
  
  //=========================================
  // New event
  //=========================================
  fData            ->Clear();
  fData->event      = evt.id().event();
  fData->run        = evt.id().run();
  fIsRealData       = evt.isRealData();
 
  // Get event time
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fData->timestamp = tts.AsDouble();

  //fData->timestamp_aux = ev.eventAuxiliary().time().timeHigh();

  fNumEvents++;
  
  // Detector properties
  auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
  // Electron lifetime
  fData->lifetime = detProp->ElectronLifetime(); 
  
  // Plane offsets. First induction plane will be treated as reference. 
  //TickOffset[0] = 0;
  //TickOffset[1] = detProp->GetXTicksOffset(1,0,0)-detProp->GetXTicksOffset(0,0,0);
  //TickOffset[2] = detProp->GetXTicksOffset(2,0,0)-detProp->GetXTicksOffset(0,0,0);
  //TickOffset[0] = detProp->GetXTicksOffset(0,0,0);
  //TickOffset[1] = detProp->GetXTicksOffset(0,0,0);
  //TickOffset[2] = detProp->GetXTicksOffset(0,0,0);
  //TickOffset[1] = detProp->GetXTicksOffset(1,0,0)-detProp->GetXTicksOffset(0,0,0);
  //TickOffset[2] = detProp->GetXTicksOffset(2,0,0)-detProp->GetXTicksOffset(0,0,0);
  
  /*
  std::cout
  <<"detProp->GetXTicksOffset(0) = "<<detProp->GetXTicksOffset(0,0,0)<<"\n"
  <<"detProp->GetXTicksOffset(1) = "<<detProp->GetXTicksOffset(1,0,0)<<"\n"
  <<"detProp->GetXTicksOffset(2) = "<<detProp->GetXTicksOffset(2,0,0)<<"\n";
  */

  // Tell us what's going on!
  std::cout<<"\n"
  <<"=========== BlipAna =========================\n"
  <<"Processing event "<<evt.id().event()<<" in run "<<evt.id().run()<<"; total: "<<fNumEvents<<"\n";
  //<<"Electron lifetime: "<<fData->lifetime<<"\n"
  //<<"TickOffsets      : "<<TickOffset[0]<<", "<<TickOffset[1]<<", "<<TickOffset[2]<<"\n";


  //=========================================
  // Get data products for this event
  //=========================================
  
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;
  
  // -- G4 particles
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // -- hits
  art::Handle< std::vector<recob::Hit> > hitlistHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitlistHandle))
    art::fill_ptr_vector(hitlist, hitlistHandle);
  
  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  
  // -- showers
  /*
  art::Handle< std::vector<recob::Shower> > showerlistHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerlistHandle))
    art::fill_ptr_vector(showerlist, showerlistHandle);
  */

  // -- hit<->track associations
  art::FindManyP<recob::Track> fmtrk(hitlistHandle,evt,fTrackModuleLabel);

  // -- hit<->shower associations
  //art::FindManyP<recob::Shower> fmshwr(hitlistHandle,evt,fShowerModuleLabel);
  
  // -- track<->calorimetry
  //art::FindMany<anab::Calorimetry> fmcal(tracklistHandle, evt, fCaloModuleLabel);

  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = (int)plist.size();
  fData->ntrks      = (int)tracklist.size();
  //fData->nshwrs     = (int)showerlist.size();
  fData->Resize();
  
  std::cout
  <<"Found "<<fData->nparticles<<" G4 particles, "
            <<fData->nhits<<" hits from "<<fHitModuleLabel
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
  sim::ParticleList::const_iterator itPart = plist.begin();
  //for(size_t i = 0; (i<plist.size())&&(itPart!=plist.end()); i++){
  for(size_t i = 0; i<nParticles; i++){
    const simb::MCParticle* pPart = (itPart++)->second;
    
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
      PrintParticleInfo(i);
    }

    //if( pPart->PdgCode() == 11 ) h_energyDep_electrons->Fill(pinfo[i].energyDep);
    //if( pPart->PdgCode() == 22 ) h_energyDep_gammas   ->Fill(pinfo[i].energyDep);

  } // endloop over G4 particles

  // Calculate the true blips
  BlipUtils::MakeTrueBlips(pinfo, trueblips);

  // Merge and save true blip information
  MergeTrueBlips(trueblips, fTrueBlipMergeDist); 
  fData->nedeps = (int)trueblips.size();
  std::cout<<"Found "<<trueblips.size()<<" true edeps:\n";
  for(size_t i=0; i<trueblips.size(); i++ ) {
    trueblips[i].ID       = i;
    //h_trueBlipEw          ->Fill(trueblips.at(i).Energy,trueblips.at(i).Energy);
    fData->edep_tpc[i]    = trueblips.at(i).TPC;
    fData->edep_energy[i] = trueblips.at(i).Energy;
    fData->edep_charge[i] = trueblips.at(i).NumElectrons;
    fData->edep_ds[i]     = trueblips.at(i).Length;
    fData->edep_x[i]      = trueblips.at(i).Position.X();
    fData->edep_y[i]      = trueblips.at(i).Position.Y();
    fData->edep_z[i]      = trueblips.at(i).Position.Z();
    fData->edep_g4id[i]   = trueblips.at(i).LeadG4ID;
    fData->edep_pdg[i]    = trueblips.at(i).LeadG4PDG;
    std::cout
    <<"   ~ "<<i<<"  "<<trueblips.at(i).Energy<<" MeV, "
    <<" ds= "<<trueblips.at(i).Length<<" cm, "
    <<" trkID= "<<trueblips.at(i).LeadG4ID<<", pdg "<<trueblips.at(i).LeadG4PDG<<"\n";
  }


  //====================================
  // Save hit information
  //====================================
  std::vector<BlipUtils::HitInfo> hitinfo(hitlist.size());
  std::map<int,std::vector<int>> wirehitsMap;
  std::cout<<"Looping over the hits\n";
  for(size_t i=0; i<hitlist.size(); i++){
    int   wire  = hitlist[i]->WireID().Wire;
    int   plane = hitlist[i]->WireID().Plane;
    int   tpc   = hitlist[i]->WireID().TPC;
    hitinfo[i].hitid      = i;
    hitinfo[i].hit        = hitlist[i];
    hitinfo[i].wire       = wire;
    hitinfo[i].tpc        = tpc;
    hitinfo[i].plane      = plane;
    hitinfo[i].driftTicks = hitlist[i]->PeakTime() - detProp->GetXTicksOffset(plane,0,0);
    hitinfo[i].qcoll      = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),plane);

    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(hitlist[i]) ){
      BlipUtils::HitTruth( hitlist[i], hitinfo[i].g4id, hitinfo[i].g4frac, hitinfo[i].g4energy, hitinfo[i].g4charge);
      hitinfo[i].g4ids = BlipUtils::HitTruthIds(hitlist[i]);
      hitinfo[i].isreal = (hitinfo[i].g4id > 0);
    }
   
   // Find associated track
    if (fmtrk.isValid()){ 
      if (fmtrk.at(i).size())  hitinfo[i].trkid = fmtrk.at(i)[0]->ID();
    }
    
    // Find associated shower
    //if (fmshwr.isValid()){
    //  if (fmshwr.at(i).size()) hitinfo[i].shwrid = fmshwr.at(i)[0]->ID();
    //}
    
    // Add to the wire-by-wire hit map
    wirehitsMap[wire].push_back(i);
    
    //printf("  %i   tpc: %i  plane: %i  wire: %i  time: %f  rms: %f   ratio: %f\n",
    //  (int)i, tpc, plane, wire, hitinfo[i].driftTicks, hitlist[i]->RMS(), hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
    //PrintHitInfo(hitinfo[i]);

  }
 
  /*
  // Time matching
  for(size_t i=0; i<hitlist.size(); i++){
    for(size_t j=i+1; j<hitlist.size(); j++){
      if( hitinfo[i].plane == hitinfo[j].plane ) continue;
      if( hitinfo[i].tpc != hitinfo[j].tpc   ) continue;
      if( hitinfo[i].trkid >= 0 || hitinfo[j].trkid >= 0) continue;
      float dT = fabs(hitinfo[j].driftTicks - hitinfo[i].driftTicks);
      //h_hit_dT->Fill( dT );
      if( dT <= fHitMatchTolerance ) {
        hitinfo[i].ismatch = true;
        hitinfo[j].ismatch = true;
      }
    }
  }
  */

  
  

 /*
  // Flag real hits and fill hit diagnostic histograms
  int   nhits_pl[3][2]={};
  float total_chargeInHits[3] = {0};
  for(size_t i=0; i<hitlist.size(); i++){
    int   pl      = hitlist[i]->WireID().Plane; 
    int   ir      = hitinfo[i].isreal;
    float q       = hitinfo[i].qcoll;
    float qtrue   = hitinfo[i].g4charge;
    nhits_pl[pl][ir]++;
    fNumHits_pl[pl][ir]++;
    if(ir) total_chargeInHits[pl] += qtrue;
    if( fSaveDiagHistos ) {
      h_hitph[pl][ir]         ->Fill(hitlist[i]->PeakAmplitude());
      h_hitrms[pl][ir]        ->Fill(hitlist[i]->RMS());
      //h_hitrms_vs_ph[pl][ir]  ->Fill(hitlist[i]->RMS(),hitlist[i]->PeakAmplitude());
      h_hitratio[pl][ir]->Fill(hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
      h_hitint[pl][ir]->Fill(hitlist[i]->Integral());
      //h_hitarea_vs_sumadc[pl][ir]->Fill(integral, sumadc);
      //if(integral>0) h_sumadc_res[pl][ir]    ->Fill( (sumadc-integral)/integral );
      if(ir && q>0 && qtrue>0 ){
        h_nelec_TrueVsReco[pl]->Fill(qtrue/1e3,q/1e3);
        h_nelec_Resolution[pl]->Fill((q-qtrue)/qtrue);
      }
    }
    
  }
  */

  // Flag real hits and fill hit diagnostic histograms
  int   nhits_all[3]={};
  int   nhits_pl[3][2]={};
  float total_chargeInHits[3] = {0};
  for(size_t i=0; i<hitlist.size(); i++){
    int   pl      = hitlist[i]->WireID().Plane; 
    int   ir      = hitinfo[i].isreal;
    float q       = hitinfo[i].qcoll;
    float qtrue   = hitinfo[i].g4charge;
    nhits_all[pl]++;
    nhits_pl[pl][ir]++;
    fNumHits_pl[pl][ir]++;
    if(ir) total_chargeInHits[pl] += qtrue;
    h_hitph[pl]   ->Fill(hitlist[i]->PeakAmplitude());
    h_hitrms[pl]  ->Fill(hitlist[i]->RMS());
    h_hitratio[pl]->Fill(hitlist[i]->RMS()/hitlist[i]->PeakAmplitude());
    h_hitint[pl]->Fill(hitlist[i]->Integral());
    if(!fIsRealData && ir && q>0 && qtrue>0 ){
      h_nelec_TrueVsReco[pl]->Fill(qtrue/1e3,q/1e3);
      h_nelec_Resolution[pl]->Fill((q-qtrue)/qtrue);
    }
  }

  

/*
  // Fill histograms
  for(size_t i=0; i<kNplanes; i++){
    int   totHits = nhits_pl[i][0]+nhits_pl[i][1];
    //std::cout<<"  charge in hits: "<<total_chargeInHits[i]<<", total charge: "<<fData->total_numElectrons<<"\n";
    float qcomp = -9, pur = -9;
    if( fData->total_numElectrons ) qcomp = total_chargeInHits[i]/fData->total_numElectrons;
    if( totHits ) pur   = nhits_pl[i][1] / float(totHits);
    h_chargecomp[i] ->Fill( qcomp );
    h_hitpur[i]     ->Fill( pur );
    std::cout<<"Hits on plane "<<i<<": "<<nhits_pl[i][1]<<" real, "<<nhits_pl[i][0]<<" fake; compl "<<qcomp<<"; purity "<<pur<<"\n";
    //if( fMakeDiagHistos ) for(size_t j=0; j<2; j++) h_nhits[i][j]->Fill(nhits_pl[i][j]);

  }
  */
  
  for(size_t i=0; i<kNplanes; i++){
    int   totHits = nhits_all[i];
    float qcomp = -9, pur = -9;
    //h_nhits[i]->Fill(nhits_all[i]);
    if( !fIsRealData ) {
      if( fData->total_numElectrons ) qcomp = total_chargeInHits[i]/fData->total_numElectrons;
      if( totHits ) pur   = nhits_pl[i][1] / float(totHits);
      h_chargecomp[i] ->Fill( qcomp );
      h_hitpur[i]     ->Fill( pur );
      std::cout<<"Hits on plane "<<i<<": "<<nhits_pl[i][1]<<" real, "<<nhits_pl[i][0]<<" fake; compl "<<qcomp<<"; purity "<<pur<<"\n";
    }
    //if( fMakeDiagHistos ) for(size_t j=0; j<2; j++) h_nhits[i][j]->Fill(nhits_pl[i][j]);

  }
  
  
  //====================================
  // Save track information
  //====================================
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
  //std::vector<bool> hitPassesCuts(hitlist.size(), false);
  //std::vector<bool> hitIsClustered(hitlist.size(),false);
      
      
  for(size_t i=0; i<hitlist.size(); i++){
    int plane = hitlist[i]->WireID().Plane;
    
    // Exclude any hits that were part of tracks
    //if( hitinfo[i].trkid >= 0 ) continue;
    
    if( !fDoHitFiltering ) {
      hitinfo[i].isgood = true;
      continue;
    }

    // Hit goodness of fit
    if( hitlist[i]->GoodnessOfFit() < fMinHitGOF[plane] ) continue;
    
    // Check hit widths and amplitudes
    float rms = hitlist[i]->RMS();
    if( rms < fMinHitRMS[plane] || rms > fMaxHitRMS[plane] ) continue;
    hitinfo[i].rmsCut = true;
    
    // Ratio cut
    float ratio = hitlist[i]->RMS()/hitlist[i]->PeakAmplitude();
    if( ratio < fMinHitRatio[plane] || ratio > fMaxHitRatio[plane] ) continue;
    hitinfo[i].ratioCut = true;
    
    hitinfo[i].isgood = true;
    
    // since this hit looks good, also flag any overlapping hits
    // on the same wire as good too
    for(auto j : wirehitsMap[hitlist[i]->WireID().Wire] ) {
      if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) ) 
        hitinfo[j].isgood = true;
    }
  
  }
  

  // ---------------------------------------------------
  // Create collection of hit clusters on same wires
  // ---------------------------------------------------
  std::vector<BlipUtils::HitClust> hitclust;
  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
  std::cout<<"Doing first-pass cluster reco\n";
  for(auto const& wirehits : wirehitsMap){
    for(auto const& hi : wirehits.second ){
      if( !hitinfo[hi].isgood || hitinfo[hi].isclustered ) continue;
      // cluster this hit
      BlipUtils::HitClust hc = BlipUtils::MakeHitClust(hitinfo[hi]);
      if( !hc.isValid ) continue;
      hitinfo[hi].isclustered=true;
      // see if we can add other hits to it; continue until 
      // no new hits can be lumped in with this clust
      int hitsAdded;
      do{
        hitsAdded = 0;  
        for(auto const& hj : wirehits.second ) {
          if( !hitinfo[hj].isgood || hitinfo[hj].isclustered ) continue;
          if( hitlist[hj] == hitlist[hi] ) continue;
          float range = fHitClustWidthFact*hitlist[hj]->RMS();
          float t1 = hitinfo[hj].driftTicks - range;
          float t2 = hitinfo[hj].driftTicks + range;
          if( (t1 > hc.StartTime && t1 < hc.EndTime )
            ||(t2 > hc.StartTime && t2 < hc.EndTime ) ){
            BlipUtils::GrowHitClust(hitinfo[hj],hc);
            hitinfo[hj].isclustered = true;
            hitsAdded++;
          }
        }
      } while ( hitsAdded!=0 );
      hitclust.push_back(hc);
      tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(hitclust.size()-1);
    }
  }
  std::cout<<"Reconstructed "<<hitclust.size()<<" hit clusters\n";
  //for(auto const& hc : hitclust ) PrintClusterInfo(hc);

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
        h_clust_nwires->Fill(hc.Wires.size());
        if( (int)hc.Wires.size() <= fMaxWiresInCluster ) {
          hitclust_merged.push_back(hc);
          tpc_planeclustsMap_merged[hc.TPC][hc.Plane].push_back(hitclust_merged.size()-1);
        }
      }
    }
  }
  hitclust = hitclust_merged;
  tpc_planeclustsMap = tpc_planeclustsMap_merged;
  

  //--------------------------------------------
  // Save hit cluster info
  //--------------------------------------------
  fData->nclusts = (int)hitclust.size();
  for(size_t i=0; i<hitclust.size(); i++){
    hitclust[i].ID            = i;
    fData->clust_tpc[i]       = hitclust[i].TPC;
    fData->clust_plane[i]     = hitclust[i].Plane;
    fData->clust_wire[i]      = hitclust[i].LeadHit->WireID().Wire;
    fData->clust_nwires[i]    = (int)hitclust[i].Wires.size();
    fData->clust_nhits[i]     = (int)hitclust[i].HitIDs.size();
    fData->clust_lhit_id[i]   = hitclust[i].LeadHitID;
    fData->clust_lhit_ph[i]   = hitclust[i].LeadHit->PeakAmplitude();
    fData->clust_lhit_rms[i]  = hitclust[i].LeadHit->RMS();
    fData->clust_lhit_peakT[i]= hitclust[i].LeadHit->PeakTime();
    fData->clust_lhit_gof[i]  = hitclust[i].LeadHit->GoodnessOfFit();
    fData->clust_lhit_time[i] = hitclust[i].LeadHitTime;
    fData->clust_charge[i]    = hitclust[i].Charge;
    fData->clust_time[i]      = hitclust[i].Time;
    fData->clust_time_w[i]    = hitclust[i].WeightedTime;
    fData->clust_time_err[i]  = hitclust[i].TimeErr;
    fData->clust_startTime[i]=hitclust[i].StartTime;
    fData->clust_endTime[i] =hitclust[i].EndTime;
    fData->clust_timespan[i]=(hitclust[i].EndTime-hitclust[i].StartTime);
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
        break;
      }
    }
  }
  
  std::cout<<"After merging: "<<hitclust.size()<<" hit clusts:\n";
  //for(auto const& hc : hitclust ) PrintClusterInfo(hc);


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
      std::cout<<"plane "<<fCaloPlane<<" has "<<planeMap[planeA].size()<<" clusts\n";
      for(auto const& i : hitclusts_planeA ) {
       
        //std::cout<<"Looking for matches to calo plane cluster: "; PrintClusterInfo(hitclust[i]);

        // store the best-matched hitclusts
        std::vector<float> best_dT(kNplanes, 999);
        std::map<int,std::set<int>> cands;
        std::map<int,int> bestMatchedClusts;
        
        // loop over other planes
        for(auto  hitclusts_planeB : planeMap ) {
          int planeB = hitclusts_planeB.first;
          if( planeB == planeA ) continue;
          
          //std::cout
          //<<"--- checking plane "<<planeB<<": "<<hitclusts_planeB.second.size()<<" clusts\n";

          // Loop over all non-matched clusts on this plane
          for(auto const& j : hitclusts_planeB.second ) {
            if( hitclust[j].isMatched ) continue;

            //std::cout<<"    cl "<<j<<"\n";
            // Check if the two clusts overlap by any degree
            // (we don't want to waste time making calls to the geometry 
            // service if there's no chance there's a match)
            if( BlipUtils::DoHitClustsOverlap(hitclust[i],hitclust[j])) {
              
              //std::cout<<"    --> in neighborhood\n";

              // Check if the two channels actually intersect
              double y,z;
              const int chA = hitclust[i].LeadHit->Channel();
              const int chB = hitclust[j].LeadHit->Channel();
              bool doTheyCross = art::ServiceHandle<geo::Geometry>()
                                  ->ChannelsIntersect(chA,chB,y,z);
              
              if( doTheyCross ) {
                
                /*
                // Use weighted times for the match
                float tA = hitclust[i].WeightedTime;
                float tB = hitclust[j].WeightedTime;
                float dT = fabs(tA-tB);
                float maxWidth = std::max(hitclust[i].TimeErr,hitclust[j].TimeErr);
                float matchTol = std::min(fHitMatchMaxTicks, fHitMatchWidthFact*maxWidth);
                //h_hit_dt->Fill(dT);
                //h_hit_dtfrac->Fill(dT/maxWidth);
                if( dT < matchTol ){
                  cands[planeB].insert(j);
                  if( dT < best_dT[planeB] ) {
                    //std::cout<<"    --> matching hit found! "<<dT<<"\n";
                    //update vectors
                    best_dT[planeB] = dT;
                    bestMatchedClusts[planeB] = j;
                  }
                }
                */
                
                //Now check the time separation for each hit; if any are within 
                //the maximum dT threshold, the clust, is a match candidate.
                std::set<int> hitsA = hitclust[i].HitIDs;
                std::set<int> hitsB = hitclust[j].HitIDs; 
                for(auto hitA : hitsA){
                  for(auto hitB : hitsB){
                    float tA = hitinfo[hitA].driftTicks;
                    float tB = hitinfo[hitB].driftTicks;
                    float dT = fabs(tA-tB);
                    float maxWidth = std::max(hitlist[hitA]->RMS(),hitlist[hitB]->RMS());
                    float matchTol = std::min(fHitMatchMaxTicks, fHitMatchWidthFact*maxWidth);
                   
                    // only plot dT if this is the *only* hit in the cluster
                    if( hitclust[i].HitIDs.size() == 1 && hitclust[j].HitIDs.size() == 1 ) {
                      h_hit_dt->Fill( dT );
                      h_hit_dtfrac->Fill( dT / maxWidth );
                    }
                    
                    //std::cout<<"dT: "<<dT<<"  matchTolerance: "<<matchTol<<"\n";
                    if( dT < matchTol ) {
                      cands[planeB].insert(j);
                      if( dT < best_dT[planeB] ) {
                        //std::cout<<"    --> matching hit found! "<<dT<<"\n";
                        //update vectors
                        best_dT[planeB] = dT;
                        bestMatchedClusts[planeB] = j;
                      }
                    }//endif dT < hitMatchTol
                  }//endloop over hits in clusterB
                }//endloop over hits in clusterA (calo plane)
                

              }//endif valid intersection
            }//endif overlap check
          }//endloop over B clusters
        }//endloop over other planes
        
        //std::cout
        //<<"    Found "<<cands.size()<<" planes with candidate matches\n";
        if( (int)cands.size()+1 >= fMinMatchedPlanes ) {
          std::vector<BlipUtils::HitClust> hcGroup;
          hcGroup.push_back(hitclust[i]);
          hitclust[i].isMatched = true;
          for(auto hit : hitclust[i].HitIDs) hitinfo[hit].ismatch = true;
          for(auto c : cands ) {
            
            h_nmatches[c.first]->Fill(c.second.size());
            int bestID = bestMatchedClusts[c.first];
            //std::cout<<"    plane "<<c.first<<" had "<<c.second.size()<<" matches, best was "<<bestID<<"\n";
          
            // For now: pick only the match closest in dT
            // To-do: for ambiguous cases where one plane has multiple possible matches,
            // pick the match that is consistent with a match on the other (non-calo) plane
            if( bestID >= 0 ){
              hcGroup.push_back(hitclust[bestID]);
              hitclust[bestID].isMatched = true;
              for(auto hit : hitclust[bestID].HitIDs) hitinfo[hit].ismatch = true;
            }

          }
        
          //std::cout<<"    ==> Making a 3D blip out of a group of "<<hcGroup.size()<<" clusters\n";
          BlipUtils::Blip newBlip = BlipUtils::MakeBlip(hcGroup);
          if(     newBlip.isValid 
              &&  newBlip.MaxIntersectDiff <= fDiffTolerance )
            {
            //std::cout<<"New blip created with max intersection difference of "<<newBlip.MaxIntersectDiff<<"\n";
            bool isHuge = false;
            for(size_t m=0; m<kNplanes; m++){ 
              if( newBlip.Charge[m] > 500e3 ) isHuge = true; 
            }
            if (!isHuge) {
              blips.push_back(newBlip);
              for(auto hc : hcGroup ) hitclust[hc.ID].BlipID = blips.size()-1;
            }
          }
        
        }//endif check for existence of candidates on other planes

      }//endloop over caloplane clusters
    }//endif calo plane has clusters
  }//endloop over TPCs
  

  // -----------------------------------------------------
  // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
  // like associating blip with some nearby track/shower and using its tagged T0)
  //    Method 1: Assume a dE/dx = 2 MeV/cm for electrons, use that + local E-field
  //              calculate recombination.
  //    Method 2: ESTAR lookup table method ala ArgoNeuT
  for(size_t i=0; i<blips.size(); i++){


    //auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    //auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    float qColl = blips[i].Charge[fCaloPlane];
    float td    = blips[i].DriftTime;
    float depEl = qColl * exp( td / detProp->ElectronLifetime() ); 
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
    blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 

  }
 

  fData->total_blip_energy = 0;
 

  // Save blip info to tree
  fData->nblips = blips.size();
  for(size_t i=0; i<blips.size(); i++){
    fNum3DBlips++;
    auto const& b = blips[i];
    fData->blip_tpc[i]        = b.TPC;
    fData->blip_nplanes[i]    = b.NPlanes;
    //fData->blip_caloplane[i]  = fCaloPlane;
    fData->blip_charge[0][i]    = b.Charge[0];
    fData->blip_charge[1][i]    = b.Charge[1];
    fData->blip_charge[2][i]    = b.Charge[2];
    //if( fData->blip_charge[2][i]  > 1e6 ) std::cout << " watermelon "<<b.Charge[2]<<"    "<<fData->blip_charge[2][i]<<"\n";
    //fData->blip_charge[i]     = b.Charge[fCaloPlane];
    fData->blip_energy[i]     = b.Energy;
    fData->blip_energyESTAR[i]= b.EnergyESTAR;
    fData->blip_maxdiff[i]    = b.MaxIntersectDiff;
    fData->blip_x[i]          = (float)b.Position.X();
    fData->blip_y[i]          = (float)b.Position.Y();
    fData->blip_z[i]          = (float)b.Position.Z();
    
    h_blip_zy->Fill(b.Position.Z(),b.Position.Y());
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
  }
  
  std::cout<<"Reconstructed "<<blips.size()<<" blips:\n";
  /*
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
  */

  // Update clust data in Tree, so we can tie reconstructed hit clusters
  // to the blips they were grouped into
  for(size_t i=0; i<hitclust.size(); i++){
    fData->clust_ismatch[i] = hitclust[i].isMatched;
    fData->clust_blipid[i] = hitclust[i].BlipID;
  }

  // ========================================
  // Evaluate hit cuts
  // ========================================
  for(size_t i=0; i<hitlist.size(); i++){
    int plane = hitinfo[i].plane; 
    
    for(size_t j=0; j<2; j++){

      // First ALL hits
      // then only real ones
      int k = 0;
      if( j == 1 ){
        if ( hitinfo[i].isreal )  k = 1;
        else                      k = 2;
      } 
      

      h_hitCuts[plane][k]->Fill("total",1);           
      
      if( hitinfo[i].trkid < 0)   { h_hitCuts[plane][k]->Fill("untracked",1);
      if( hitinfo[i].rmsCut )     { h_hitCuts[plane][k]->Fill("rms",1);
      if( hitinfo[i].ratioCut )   { h_hitCuts[plane][k]->Fill("ratio",1);
      if( hitinfo[i].ismatch )    { h_hitCuts[plane][k]->Fill("matched",1);
      if( hitinfo[i].blipid >= 0) { h_hitCuts[plane][k]->Fill("blip",1); 
      }}}}}

    }//endloop j
   
    // Save TTree data for hits (doing this here so we can
    // incorporate blip and cluster associations that were found
    // through the blip reconstruction process)
    if( i < kMaxHits ) {
      fData->hit_plane[i]   = hitinfo[i].plane;
      fData->hit_wire[i]    = hitinfo[i].wire;
      fData->hit_tpc[i]     = hitinfo[i].tpc;
      fData->hit_trkid[i]   = hitinfo[i].trkid;
      fData->hit_shwrid[i]  = hitinfo[i].shwrid;
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_gof[i]     = hitlist[i]->GoodnessOfFit();
      fData->hit_rms[i]     = hitlist[i]->RMS();
      fData->hit_ph[i]	    = hitlist[i]->PeakAmplitude();
      fData->hit_area[i]    = hitlist[i]->Integral();
      fData->hit_sumadc[i]  = hitlist[i]->SummedADC();
      fData->hit_mult[i]    = hitlist[i]->Multiplicity();
      fData->hit_time[i]    = hitinfo[i].driftTicks;
      fData->hit_charge[i]  = hitinfo[i].qcoll;
      fData->hit_isreal[i]  = hitinfo[i].isreal;
      fData->hit_g4id[i]    = hitinfo[i].g4id;
      fData->hit_g4frac[i]  = hitinfo[i].g4frac;
      fData->hit_g4energy[i]= hitinfo[i].g4energy;
      fData->hit_g4charge[i]= hitinfo[i].g4charge;
      fData->hit_ismatch[i] = hitinfo[i].ismatch;
      fData->hit_isgood[i]  = hitinfo[i].isgood;
      fData->hit_blipid[i]  = hitinfo[i].blipid;
      fData->hit_clustid[i] = hitinfo[i].clustid;
      if( hitinfo[i].blipid >= 0 )  h_hitgof_3D[hitinfo[i].plane]->Fill(fData->hit_gof[i]);
      else                          h_hitgof_2D[hitinfo[i].plane]->Fill(fData->hit_gof[i]);
    }
  
  }


  //====================================
  // Fill TTree
  //====================================
  fData->FillTree();
}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipAna::endJob(){
 
  // Divide hitcut histos
  for(size_t i=0; i<kNplanes; i++) {
    h_hitCuts[i][0]->Scale(1./fNumEvents);
    h_hitCuts[i][1]->Scale(1./fNumEvents);
    h_hitCuts[i][2]->Scale(1./fNumEvents);
  }

  printf("\n=============================================\n");
  printf("BlipAna Summary\n\n");
  printf("  Req num matched planes  : %i\n",        fMinMatchedPlanes);
  printf("  Max match window        : %3.1f ticks\n",fHitMatchMaxTicks);
  printf("  Hit match width fact    : x%3.1f\n",    fHitMatchWidthFact);
  printf("\n");
  printf("  Total events            :   %i\n\n",fNumEvents);
  printf("  Blips per event         : %6.2f\n",fNum3DBlips/float(fNumEvents));
  printf("\n");
  for(size_t i=0; i<kNplanes; i++){
  printf("   Plane %lu -------------------------\n",i);
  printf("   - ave total fake hits  : %.2f\n",fNumHits_pl[i][0]/(float)fNumEvents);
  printf("   - ave total real hits  : %.2f\n",fNumHits_pl[i][1]/(float)fNumEvents);
  printf("   - charge completeness  : %f\n",h_chargecomp[i]->GetMean());
  printf("   - hit purity           : %f\n",h_hitpur[i]->GetMean());
  } 
  printf("\n=============================================\n");

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
  printf("  hitID: %4i, TPC: %i, plane: %i, driftTicks: %7.2f, leadWire: %3i, isreal: %2i, isgood: %2i, G4ID: %4i, recoTrack: %4i\n",
    hi.hitid,
    hi.tpc,
    hi.plane,
    hi.driftTicks,
    hi.wire,
    (int)hi.isreal,
    (int)hi.isgood,
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
