////////////////////////////////////////////////////////////////////////////////
// Class:       gammacorrelation
// Plugin Type: analyzer (art v3_01_02)
// File:        gammacorrelation_module.cc
//
// Generated at Wed Mar 20 10:42:52 2019 by Avinay Bhat using cetskelgen
// from cetlib version v3_05_01.
// Avinay Bhat (avbhat@syr.edu)
// This version  was modified by Ohana Benevides (obenevid@syr.edu) on May 2021 
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/ProcessHistoryID.h"
#include "lardata/Utilities/FindManyInChainP.h"

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "canvas/Persistency/Provenance/BranchType.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// save info associated to common optical filter
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "ubevt/Utilities/PMTRemapService.h"
#include "ubevt/Utilities/PMTRemapProvider.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include <TTree.h>

#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

// For the NuMI trigger objects
#include "lardataobj/RawData/TriggerData.h" 
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h" 

#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TRandom3.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <memory>
/*
#include "larreco/RecoAlg/TrackMomentumCalculator.h" // Momentum from range
#include "larreco/RecoAlg/TrajectoryMCSFitter.h" // Multiple Coulomb Scattering (MCS)

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

*/

class gammacorrelation : public art::EDAnalyzer {

public:
  explicit gammacorrelation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  
  gammacorrelation(gammacorrelation const&) = delete;
  gammacorrelation(gammacorrelation&&) = delete;
  gammacorrelation& operator=(gammacorrelation const&) = delete;
  gammacorrelation& operator=(gammacorrelation&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;
  void endSubRun(art::SubRun const& sr) override;

  //-------------------------------------------------//

public:

	void mudar_analyze(art::Event const& e);
  void gammacorrelation_analyze(art::Event const& e);

  void mudar_beginJob();
  void gammacorrelation_beginJob();

  void mudar_endSubRun(art::SubRun const& sr);
  void gammacorrelation_endSubRun(art::SubRun const& sr);

private:
  // Declare member data here.

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  std::string fHit_tag,fpfparticle_tag,fvertex_tag,fsps_tag,fcluster_tag,fReco_track_tag;

  //std::string fMCProducer_tag;

  Int_t evttime=0;

  Double32_t sps_x,sps_y,sps_z,sps_hit_charge,sps_cluster_charge,sps_cluster_charge10,sps_cluster_charge20,sps_cluster_charge50,sps_cluster_charge75,sps_cluster_charge100,sps_cluster_charge150,sps_cluster_charge200;
  Double_t distance, distance_smallest,Event_cluster_charge;//distance_smallest is the distance between a spacepoint and the vertex
  Double_t Vertex_x,Vertex_y,Vertex_z;
  TRandom3 rand;
  Double_t _rand_vtx_x, _rand_vtx_y, _rand_vtx_z, distance_rand_vtx, distance_smallest_rand_vtx;
  Int_t neutrinos,N_sps,N_Event,N_Run,N_SubRun,N_sps10,N_sps20,N_sps50,N_sps75,N_sps100,N_sps150,N_sps200;
  float _maxTrkLen;
  int   _neutrinoshowers;
  int   _neutrinotracks;
  float _muon_px, _muon_py, _muon_pz;
  Double_t tracklength = 0.0,X_reco=0.0,Y_reco=0.0,Z_reco=0.0;
  Double_t pointdistance_nu=0.0;
  Double_t pointdistance_nu_cosmic_smallest=0.0;
  Double_t distance_nu_cosmic_smallest=0.0;

  Double_t pointdistance_trk=0.0;
  Double_t pointdistance_trk_smallest=0.0;
  Double_t distance_trk_smallest=0.0;

  Double_t pointdistance_smallest_nu;
  Double_t track_point_length=0.0;
  Double_t track_point_length_smallest=0.0;
  Double_t distance_smallest_nu; //distance between a neutrino correlated track and a spacepoint
  Double_t X_reco_nu;
  Double_t Y_reco_nu;
  Double_t Z_reco_nu;
  Double_t X_reco_smallest_nu=0;
  Double_t Y_reco_smallest_nu=0;
  Double_t Z_reco_smallest_nu=0;
  Double_t X_reco_best_nu=0.0;
  Double_t Y_reco_best_nu=0.0;
  Double_t Z_reco_best_nu=0.0;
  Double_t R=2.0; //Radius of the cone for a neutrino correlated track
  Double_t H=4.0; //Height of the cone for a neutrino correlated track

  Int_t cosmic_trk_50=0;

  TTree *Event_Correlationtree;
  TTree *Sps_Correlationtree;

  bool isData,fData;

  /** Setup root trees  */
  TTree *potTree;
  double sr_pot = 0;
  int sr_run = 0;
  int sr_sub_run = 0;

  Double_t fidVolMinX =    0; //Fiducial Volume dimensions for MicroBooNE
  Double_t fidVolMaxX =  256;
  Double_t fidVolMinY = -116;
  Double_t fidVolMaxY =  116;
  Double_t fidVolMinZ =    0;
  Double_t fidVolMaxZ = 1030;

private:
  //*********************** From muDAR_feasibility ***********************
  //Truth information variabels
	TTree* _tree_Events;
	double _neutrino_E; //The neutrino energy 
	double _lepton_E; //The son lepton energy
  int _neutrino_mother_PDG; //The mother particle PDG (hopefully all muons)

	double _vertex_X_outgoing_lepton; //The X position of the lepton vertex 
	double _vertex_Y_outgoing_lepton; //The Y position of the lepton vertex   
	double _vertex_Z_outgoing_lepton; //The Z position of the lepton vertex 

  double _vertex_X_incoming_nu; //The X position of the incoming nu vertex 
	double _vertex_Y_incoming_nu; //The Y position of the incoming nu vertex   
	double _vertex_Z_incoming_nu; //The Z position of the incoming nu vertex 
	
  //Flux information variabels
  TTree* _flux_tree;
  int _which_decay_type; //Type of decay. For NuMI, you can find the decay codes @ http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/
  double _decay_x; //The decay x position. The origin is defined as the target. 
  double _decay_y; //The decay y position. The origin is defined as the target. 
  double _decay_z; //The decay z position. The origin is defined as the target. 

  int _mcflux_size; //Size of the MCFlux vector -- This is me trying to make sure this is all 1 and I don't need to use associations  


	// POT information variabels  
	TTree *_subrun_tree;
  int _run_sr; // The run number
	int _sub_sr; // The subRun number
  float _pot;  // The total amount of POT for the current sub run
	
	// Trigger (common optical filtervariables)                                                        
  float  _opfilter_pe_beam, _opfilter_pe_beam_tot, _opfilter_pe_veto, _opfilter_pe_veto_tot;

	// NuMI Trigger (common optical filtervariables)
	bool _swtrig_pre_ext, _swtrig_pre, _swtrig_post_ext, _swtrig_post;

  //***********************************************************************
};

gammacorrelation::gammacorrelation(fhicl::ParameterSet const& p)
: EDAnalyzer{p}  // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  //All tags get filled in the run_gammacorrelation.fcl file

  //fMCProducer_tag = p.get<std::string>("MCProducer_tag");
  fHit_tag = p.get<std::string>("hit_tag"  );
  fpfparticle_tag=p.get<std::string>("pfparticle_tag");
  fvertex_tag=p.get<std::string>("vertex_tag");
  fsps_tag=p.get<std::string>("sps_tag");
  fcluster_tag=p.get<std::string>("cluster_tag");
  fData     = p.get< bool >("IsData");
  fReco_track_tag = p.get<std::string>("recotrack_tag"  );
}

void gammacorrelation::analyze(art::Event const& e){

mudar_analyze(e);
gammacorrelation_analyze(e);
}

void gammacorrelation::gammacorrelation_analyze(art::Event const& e){
// art::EventAuxiliary const& ev_aux;

  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() << " SubRun Number: " << e.subRun() <<"**********************"<< std::endl;

  art::Timestamp ts = e.time();
  evttime = ts.timeHigh();

  // create list of tracks and showers associated to neutrino candidate
  std::vector<recob::Track  > sliceTracks;
  std::vector<recob::Shower > sliceShowers;

  N_Event=e.event();
  N_Run=e.run();
  N_SubRun=e.subRun();

  isData = fData;

  art::Handle<std::vector<recob::Hit> > hit_handle;
  e.getByLabel(fHit_tag,hit_handle);

  art::Handle<std::vector<recob::Cluster> > cluster_handle;
  e.getByLabel(fcluster_tag,cluster_handle);

  art::Handle<std::vector<recob::PFParticle> > pfparticle_handle;
  e.getByLabel(fpfparticle_tag,pfparticle_handle);

  art::Handle<std::vector<recob::Track> > recotrack_handle;
  e.getByLabel(fReco_track_tag,recotrack_handle);

  art::FindMany<recob::Track> pfp_track_assn_v  (pfparticle_handle, e, fpfparticle_tag);
  art::FindMany<recob::Shower> pfp_shower_assn_v(pfparticle_handle, e, fpfparticle_tag);

  art::Handle<std::vector<recob::Vertex> > vertex_handle;
  e.getByLabel(fvertex_tag,vertex_handle);

  art::FindMany<recob::Vertex> pfp_vertex_assn_v(pfparticle_handle, e, fvertex_tag);

  art::Handle<std::vector<recob::SpacePoint> > spacepoint_handle;
  e.getByLabel(fsps_tag,spacepoint_handle);

  art::FindMany<recob::Cluster> sps_clus_assn_v(spacepoint_handle, e, fsps_tag);
  art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_handle, e, fcluster_tag);

  // build PFParticle map  for this event
  _pfpmap.clear();
  for (unsigned int p=0; p < pfparticle_handle->size(); p++)
  _pfpmap[pfparticle_handle->at(p).Self()] = p;
  _neutrinoshowers = 0;
  _neutrinotracks = 0;
  _maxTrkLen = 0;

  recob::Vertex nuvtx;
  TVector3 rndvtx;
  neutrinos = 0;
  Event_cluster_charge=0,sps_cluster_charge50=0,sps_cluster_charge20=0,sps_cluster_charge10=0,sps_cluster_charge75=0,sps_cluster_charge100=0,sps_cluster_charge150=0,sps_cluster_charge200=0;
  for (size_t p=0; p < pfparticle_handle->size(); p++) {
    auto pfp = pfparticle_handle->at(p);

    if (pfp.IsPrimary() == false) continue;

    auto PDG = fabs(pfp.PdgCode());
    if ( (PDG == 12) || (PDG == 14) ) {

      neutrinos += 1;
      // grab daughter PFParticles
      auto daughters = pfp.Daughters();

      for(auto const& daughterid : daughters) {

        if (_pfpmap.find(daughterid) == _pfpmap.end()) {
          std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
          continue;
        }

        const auto daughter = pfparticle_handle->at( _pfpmap[daughterid] );
        // if there is a track associated to the PFParticle, add it
        auto const& ass_trk_v = pfp_track_assn_v.at( _pfpmap[daughterid] );
        if (ass_trk_v.size() == 1) {
          sliceTracks.push_back( *(ass_trk_v.at(0)) );
          if (ass_trk_v.at(0)->Length() > _maxTrkLen) {
            _maxTrkLen = ass_trk_v.at(0)->Length();
            _muon_px   = ass_trk_v.at(0)->StartDirection().X();
            _muon_py   = ass_trk_v.at(0)->StartDirection().Y();
            _muon_pz   = ass_trk_v.at(0)->StartDirection().Z();
          }// if longest track
        }
        // if there is a shower associated to the PFParticle, add it
        auto const& ass_shr_v = pfp_shower_assn_v.at( _pfpmap[daughterid] );
        if (ass_shr_v.size() == 1) sliceShowers.push_back( *(ass_shr_v.at(0)) );
      }// for all PFParticles in the slice

      auto ass_vtx_v  =pfp_vertex_assn_v.at( p );
      if (ass_vtx_v.size() != 1)
      std::cout << "ERROR. Neutrino not associated with a single vertex..." << std::endl;
      nuvtx = *(ass_vtx_v.at(0));

      _neutrinoshowers = sliceShowers.size();
      _neutrinotracks  = sliceTracks.size();

    }
  }
  distance_smallest=1e10;

  distance_smallest_rand_vtx=1e10;
  N_sps=0,N_sps10=0,N_sps20=0,N_sps50=0,N_sps75=0,N_sps100=0,N_sps150=0,N_sps200=0;
  Vertex_x=-9999.0;
  Vertex_y=-9999.0;
  Vertex_z=-9999.0;

  // random generated "neutrino" vtx in the fiducial volume
  _rand_vtx_x = rand.Uniform(fidVolMinX,fidVolMaxX);
  _rand_vtx_y = rand.Uniform(fidVolMinY,fidVolMaxY);
  _rand_vtx_z = rand.Uniform(fidVolMinZ,fidVolMaxZ);

  rndvtx = TVector3(_rand_vtx_x,_rand_vtx_y,_rand_vtx_z);
  distance_nu_cosmic_smallest=1e10;
  cosmic_trk_50=0;

  for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP (This is to determine if there is a cosmic muon track close to the neutrino vertex)

    auto const& track = recotrack_handle->at(i_t);

    if ( (sqrt((pow(nuvtx.position().x()-track.Start().X(),2))+(pow(nuvtx.position().y()-track.Start().Y(),2))+ (pow(nuvtx.position().z()-track.Start().Z(),2))) <5.0) ||  (sqrt((pow(nuvtx.position().x()-track.End().X(),2))+(pow(nuvtx.position().y()-track.End().Y(),2))+ (pow(nuvtx.position().z()-track.End().Z(),2)))<5.0))
    continue;

    if (track.Length()<20.0) //Targeting long cosmic muon tracks and ignoring the small tracks
    continue;

    pointdistance_nu_cosmic_smallest=1e10;

    for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP
      X_reco=track.LocationAtPoint(m).X();
      Y_reco=track.LocationAtPoint(m).Y();
      Z_reco=track.LocationAtPoint(m).Z();

      pointdistance_nu=sqrt((pow(nuvtx.position().x()-X_reco,2))+(pow(nuvtx.position().y()-Y_reco,2))+ (pow(nuvtx.position().z()-Z_reco,2)));

      if(pointdistance_nu<pointdistance_nu_cosmic_smallest){
        pointdistance_nu_cosmic_smallest=pointdistance_nu;

      }

    }
 
    if (pointdistance_nu_cosmic_smallest<50.0){
      cosmic_trk_50++;

    }

    if(pointdistance_nu_cosmic_smallest<distance_nu_cosmic_smallest){
      distance_nu_cosmic_smallest=pointdistance_nu_cosmic_smallest;
      tracklength= track.Length();

    }
   
  }

  for(size_t s=0;s<spacepoint_handle->size();s++){ //START SPACEPOINT LOOP
    N_sps++;
  
    auto sps = spacepoint_handle->at(s);
    auto cluster_v = sps_clus_assn_v.at(s);

    sps_x=sps.XYZ()[0];
    sps_y=sps.XYZ()[1];
    sps_z=sps.XYZ()[2];

    if (sps_y < -90.0 || sps_y > 90.0 || sps_z < 20 || sps_z > 1000)
    continue;

    distance_trk_smallest=1e10;
    distance_smallest_nu=1e10;

    for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP

      auto const& track = recotrack_handle->at(i_t);

      //The following IF loop is only meant for tracks coming out of the neutrino vertex
      if ( (sqrt((pow(nuvtx.position().x()-track.Start().X(),2))+(pow(nuvtx.position().y()-track.Start().Y(),2))+ (pow(nuvtx.position().z()-track.Start().Z(),2))) <5.0) ||  (sqrt((pow(nuvtx.position().x()-track.End().X(),2))+(pow(nuvtx.position().y()-track.End().Y(),2))+ (pow(nuvtx.position().z()-track.End().Z(),2)))<5.0)){

        pointdistance_smallest_nu=1e10;

        track_point_length=0.0;
        for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP

          X_reco_nu=track.LocationAtPoint(m).X();
          Y_reco_nu=track.LocationAtPoint(m).Y();
          Z_reco_nu=track.LocationAtPoint(m).Z();

          pointdistance_nu=sqrt((pow(sps_x-X_reco_nu,2))+(pow(sps_y-Y_reco_nu,2))+(pow(sps_z-Z_reco_nu,2)));

          if(pointdistance_nu<pointdistance_smallest_nu){//comparison IF loop

            pointdistance_smallest_nu=pointdistance_nu;

            X_reco_smallest_nu=X_reco_nu;
            Y_reco_smallest_nu=Y_reco_nu;
            Z_reco_smallest_nu=Z_reco_nu;
            track_point_length= sqrt((pow(X_reco_nu-nuvtx.position().x(),2))+(pow(Y_reco_nu-nuvtx.position().y(),2))+ (pow(Z_reco_nu-nuvtx.position().z(),2)));

          }
        }//END RECO POINT LOOP

        if(pointdistance_smallest_nu<distance_smallest_nu){
          distance_smallest_nu=pointdistance_smallest_nu;
          track_point_length_smallest=track_point_length;
          X_reco_best_nu=X_reco_smallest_nu; //variables for the coordinates of the nearest reco track
          Z_reco_best_nu=Z_reco_smallest_nu; //variables for the coordinates of the nearest reco track
        }
      }//END IF LOOP


      pointdistance_trk_smallest=1e10;

      for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP

        pointdistance_trk= sqrt((pow(track.Start().X()-sps_x,2))+(pow(track.Start().Y()-sps_y,2))+ (pow(track.Start().Z()-sps_z,2)));

        if (pointdistance_trk<pointdistance_trk_smallest){
          pointdistance_trk_smallest=pointdistance_trk;
        }

      }

      if (pointdistance_trk_smallest<distance_trk_smallest){
        distance_trk_smallest=pointdistance_trk_smallest;

      }

    }

    //The following IF CONTINUE condition is only meant for both neutrino correlated tracks and also cosmic tracks

    if ((distance_smallest_nu<((R/H)*(track_point_length_smallest)) && (track_point_length_smallest<H))||((distance_smallest_nu>R) && (track_point_length_smallest>H)))
    continue;

    Vertex_x=nuvtx.position().x();
    Vertex_y=nuvtx.position().y();
    Vertex_z=nuvtx.position().z();

    distance = sqrt((pow(nuvtx.position().x()-sps_x,2))+(pow(nuvtx.position().y()-sps_y,2))+ (pow(nuvtx.position().z()-sps_z,2)));
 
    sps_cluster_charge =0;

    for (auto const& cluster: cluster_v){
      auto plane = cluster->View();

      if (plane==2){
        sps_cluster_charge = cluster->Integral();
        Event_cluster_charge += sps_cluster_charge;
      }
    }

    if (distance<distance_smallest){
      distance_smallest=distance;
    }

    if (distance<10.0){
      N_sps10++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge10 += cluster->Integral();
        }
      }
    }

    if (distance<20.0){
      N_sps20++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge20 += cluster->Integral();
        }
      }
    }

    if (distance<50.0){
      N_sps50++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge50 += cluster->Integral();

        }
      }
    }

    if (distance<75.0){
      N_sps75++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge75 += cluster->Integral();
        }
      }
    }

    if (distance<100.0){
      N_sps100++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge100 += cluster->Integral();
        }
      }
    }

    if (distance<150.0){
      N_sps150++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge150 += cluster->Integral();
        }
      }
    }

    if (distance<200.0){
      N_sps200++;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();
        if (plane==2){
          sps_cluster_charge200 += cluster->Integral();

        }
      }
    }

    distance_rand_vtx= sqrt((pow(_rand_vtx_x-sps_x,2))+(pow(_rand_vtx_y-sps_y,2))+ (pow(_rand_vtx_z-sps_z,2)));
    if (distance_rand_vtx<distance_smallest_rand_vtx){
      distance_smallest_rand_vtx=distance_rand_vtx;
    }
   
    Sps_Correlationtree->Fill();
  }//END SPACEPOINT LOOP

  Event_Correlationtree->Fill();//Filled per event

}

void gammacorrelation::mudar_analyze(art::Event const& e){

	// ************************************************************************
	// *********************** Load MCTruth info *****************************
	// ************************************************************************

  art::Handle<art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>> MCTruthHandle;
	art::InputTag fMCTruthTag("largeant");
	std::cout << "Is true getbylabel: " << 
	e.getByLabel(fMCTruthTag, MCTruthHandle) << std::endl;
	std::cout << "Size: " << MCTruthHandle->size() << std::endl;
	const art::Ptr<simb::MCTruth> mctruth = MCTruthHandle->at(0).first;
	std::cout << "mctruth nparticles  " << mctruth->NParticles() << std::endl;


	// ************************************************************************
	// *********************** Load MCFlux info *****************************
	// ************************************************************************

  art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  e.getByLabel("generator", mcFluxHandle);
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;
  art::fill_ptr_vector(mcFluxVec, mcFluxHandle);

  //Retrieving flux information 

  _decay_x = mcFluxVec.at(0)->fvx; 
  _decay_y = mcFluxVec.at(0)->fvy; 
  _decay_z = mcFluxVec.at(0)->fvz; 
  _which_decay_type = mcFluxVec.at(0)->fndecay; 
  _mcflux_size = mcFluxVec.size(); 

  //Retrieving the neutrino information
	auto neutrino = mctruth->GetNeutrino();
  const simb::MCParticle nu = neutrino.Nu();
	_neutrino_E = nu.E();

  //Checking the PDG of the neutrino mother particle
  _neutrino_mother_PDG = nu.Mother();

  //For the nu vertex...
  const TLorentzVector vertex_nu = nu.Position(); 
	_vertex_X_incoming_nu = vertex_nu.X(); 
 	_vertex_Y_incoming_nu = vertex_nu.Y(); 
	_vertex_Z_incoming_nu = vertex_nu.Z(); 
	
	//Retrieving the lepton information
  const simb::MCParticle lepton = neutrino.Lepton();
	_lepton_E = lepton.E();

	//For the lepton vertex...
	const TLorentzVector vertex_lepton = lepton.Position(); 
	_vertex_X_outgoing_lepton = vertex_lepton.X(); 
 	_vertex_Y_outgoing_lepton = vertex_lepton.Y(); 
	_vertex_Z_outgoing_lepton = vertex_lepton.Z(); 
    
  // ***********************************************************************
  // ************************************************************************
  // ******** Trigger result output for NuMI software trigger *************** 
  // ******** Thanks to Owen who gave the code snippet for this! ************
  // ************************************************************************
  
	art::InputTag triggerTag ("swtrigger", "", "DataOverlayOpticalNuMI" );
  	const auto& triggerHandle = e.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);
  	std::vector<std::string> triggerName = triggerHandle->getListOfAlgorithms();

  for (int j=0; j!=triggerHandle->getNumberOfAlgorithms(); j++){

	if (triggerName[j] == "EXT_NUMIwin_FEMBeamTriggerAlgo"){
    	bool _swtrig_pre_ext = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_pre_ext: " << _swtrig_pre_ext << std::endl;
    }

    else if (triggerName[j] == "EXT_NUMIwin_2018May_FEMBeamTriggerAlgo"){
    	bool _swtrig_post_ext = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_post_ext: " << _swtrig_post_ext << std::endl;
    }

    else if (triggerName[j] == "NUMI_FEMBeamTriggerAlgo"){
    	bool  _swtrig_pre = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_pre: " << _swtrig_pre << std::endl;
    }

    else if (triggerName[j] == "NUMI_2018May_FEMBeamTriggerAlgo"){
    	bool _swtrig_post = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_post: " << _swtrig_post << std::endl;
    }

    else continue;

    // Print the trigger and the result
    std::cout<<triggerName[j]<<": ";
    std::cout<<triggerHandle->passedAlgo(triggerName[j])<<std::endl;
  }

	// ************************************************************************
	// *************** Load common-optical-filter output **********************
	// ************************************************************************
	art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
	art::InputTag fCommonOpFiltTag("opfiltercommon");

	e.getByLabel(fCommonOpFiltTag, CommonOpticalFilter_h);

	_opfilter_pe_beam     = CommonOpticalFilter_h->PE_Beam();
	_opfilter_pe_beam_tot = CommonOpticalFilter_h->PE_Beam_Total();
	_opfilter_pe_veto     = CommonOpticalFilter_h->PE_Veto();
	_opfilter_pe_veto_tot = CommonOpticalFilter_h->PE_Veto_Total();

	// Implementation of required member function here.
	// Software Trigger
	art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    std::string fSoftwareTriggerModuleLabel = "swtrigger"; 
	
	if (e.getByLabel(fSoftwareTriggerModuleLabel, softwareTriggerHandle)){
		
		std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
		
		for (int i = 0; i < int(algoNames.size()); i++){

			if (algoNames[i] == "BNB_FEMBeamTriggerAlgo"){ 
			
				auto EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[i]); 

				std::cout << "EventPassedSwTrigger: " << EventPassedSwTrigger << std::endl;

				std::cout << "Found BNB_FEMBeamTrigger" << std::endl; 
			}
		}
	}
	else if (e.getByLabel("daq", softwareTriggerHandle)){

		std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();

		for (int i = 0; i < int(algoNames.size()); i++){
			
			if (algoNames[i] == "EXT_BNBwin_FEMBeamTriggerAlgo"||algoNames[i] == "BNB_FEMBeamTriggerAlgo"){
		 		
				auto EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[i]); 

				std::cout << "EventPassedSwTrigger: " << EventPassedSwTrigger << std::endl;

				std::cout << "Found EXT_BNBwin_FEMBeamTrigger and BNB_FEMBeamTrigger" << std::endl; 
			}
		}
	}

	_tree_Events->Fill();
  _flux_tree->Fill();
}


void gammacorrelation::mudar_beginJob(){


	art::ServiceHandle<art::TFileService> tfs;

	_tree_Events = tfs->make<TTree>("Events", "Events");
	_subrun_tree = tfs->make<TTree>("SubRun", "SubRun");
  _flux_tree = tfs->make<TTree>("Flux_Events", "Flux_Events");

	//trigger branches
	_tree_Events -> Branch("_opfilter_pe_beam", &_opfilter_pe_beam, "_opfilter_pe_beam/F");
	_tree_Events -> Branch("_opfilter_pe_beam_tot", &_opfilter_pe_beam_tot, "_opfilter_pe_beam_tot/F");
	_tree_Events -> Branch("_opfilter_pe_veto", &_opfilter_pe_veto, "_opfilter_pe_veto/F");
	_tree_Events -> Branch("_opfilter_pe_veto_tot", &_opfilter_pe_veto_tot, "_opfilter_pe_veto_tot/F");
	
	//NuMI trigger branches
	_tree_Events->Branch("_swtrig_pre", &_swtrig_pre, "_swtrig_pre/O");
	_tree_Events->Branch("_swtrig_pre_ext", &_swtrig_pre_ext, "_swtrig_pre_ext/O");
	_tree_Events->Branch("_swtrig_post_ext", &_swtrig_post_ext, "_swtrig_post_ext/O");
	_tree_Events->Branch("_swtrig_post", &_swtrig_post, "_swtrig_post/O");

	//truth info particles branches
	_tree_Events -> Branch("neutrino_E", &_neutrino_E, "neutrino_E/D");
	_tree_Events -> Branch("lepton_E", &_lepton_E, "lepton_E/D");
  _tree_Events -> Branch("neutrino_mother_pdg", &_neutrino_mother_PDG, "neutrino_mother_pdg/I"); 
	_tree_Events -> Branch("vertex_X_lepton", &_vertex_X_outgoing_lepton, "vertex_X/D");
	_tree_Events -> Branch("vertex_Y_lepton", &_vertex_Y_outgoing_lepton, "vertex_Y/D");
	_tree_Events -> Branch("vertex_Z_lepton", &_vertex_Z_outgoing_lepton, "vertex_Z/D");
  _tree_Events -> Branch("vertex_X_nu", &_vertex_X_incoming_nu, "vertex_X/D");
	_tree_Events -> Branch("vertex_Y_nu", &_vertex_Y_incoming_nu, "vertex_Y/D");
	_tree_Events -> Branch("vertex_Z_nu", &_vertex_Z_incoming_nu, "vertex_Z/D");

	//POT counting info 
  _subrun_tree->Branch("run", &_run_sr, "run/I");
  _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");
	_subrun_tree->Branch("pot", &_pot, "pot/F");

  //Flux info
  _flux_tree->Branch("Decay_type", &_which_decay_type, "Decay_type/I");
  _flux_tree->Branch("Decay_X", &_decay_x, "Decay_X/D");
  _flux_tree->Branch("Decay_Y", &_decay_y, "Decay_Y/D");
  _flux_tree->Branch("Decay_Z", &_decay_z, "Decay_Z/D");
  _flux_tree->Branch("MCflux_size", &_mcflux_size, "MCflux_size/I");

}

void gammacorrelation::gammacorrelation_beginJob(){

   art::ServiceHandle<art::TFileService> tfs;

  Event_Correlationtree = tfs->make<TTree>("Event_Correlationtree",    "Event_Correlationtree");
  Event_Correlationtree->Branch("evttime",&evttime,"evttime/I");
  Event_Correlationtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  Event_Correlationtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  Event_Correlationtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  Event_Correlationtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  Event_Correlationtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  Event_Correlationtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  Event_Correlationtree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");// //distance_smallest is the distance between a spacepoint and the vertex
  Event_Correlationtree->Branch("distance_smallest_nu",&distance_smallest_nu,"distance_smallest_nu/D");//distance between a neutrino correlated track and a spacepoint
  Event_Correlationtree->Branch("N_sps",&N_sps,"N_sps/I");
  Event_Correlationtree->Branch("Event_cluster_charge",&Event_cluster_charge,"Event_cluster_charge/D");
  Event_Correlationtree->Branch("distance_smallest_rand_vtx",&distance_smallest_rand_vtx,"distance_smallest_rand_vtx/D");
  Event_Correlationtree->Branch("N_sps10",&N_sps10,"N_sps10/I");
  Event_Correlationtree->Branch("N_sps20",&N_sps20,"N_sps20/I");
  Event_Correlationtree->Branch("N_sps50",&N_sps50,"N_sps50/I");
  Event_Correlationtree->Branch("N_sps75",&N_sps75,"N_sps75/I");
  Event_Correlationtree->Branch("N_sps100",&N_sps100,"N_sps100/I");
  Event_Correlationtree->Branch("N_sps150",&N_sps150,"N_sps150/I");
  Event_Correlationtree->Branch("N_sps200",&N_sps200,"N_sps200/I");
  Event_Correlationtree->Branch("sps_cluster_charge10",&sps_cluster_charge10,"sps_cluster_charge10/D");
  Event_Correlationtree->Branch("sps_cluster_charge20",&sps_cluster_charge20,"sps_cluster_charge20/D");
  Event_Correlationtree->Branch("sps_cluster_charge50",&sps_cluster_charge50,"sps_cluster_charge50/D");
  Event_Correlationtree->Branch("sps_cluster_charge75",&sps_cluster_charge75,"sps_cluster_charge75/D");
  Event_Correlationtree->Branch("sps_cluster_charge100",&sps_cluster_charge100,"sps_cluster_charge100/D");
  Event_Correlationtree->Branch("sps_cluster_charge150",&sps_cluster_charge150,"sps_cluster_charge150/D");
  Event_Correlationtree->Branch("sps_cluster_charge200",&sps_cluster_charge200,"sps_cluster_charge200/D");
  Event_Correlationtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Event_Correlationtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Event_Correlationtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Event_Correlationtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Event_Correlationtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Event_Correlationtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Event_Correlationtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  Event_Correlationtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  Event_Correlationtree->Branch("distance_nu_cosmic_smallest" ,&distance_nu_cosmic_smallest ,"distance_nu_cosmic_smallest/D" );
  Event_Correlationtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  Event_Correlationtree->Branch("N_Run" ,&N_Run ,"N_Run/I" );
  Event_Correlationtree->Branch("N_SubRun" ,&N_SubRun ,"N_SubRun/I" );

  Sps_Correlationtree = tfs->make<TTree>("Sps_Correlationtree",    "Sps_Correlationtree");
  Sps_Correlationtree->Branch("evttime",&evttime,"evttime/I");
  Sps_Correlationtree->Branch("sps_x",&sps_x,"sps_x/D");
  Sps_Correlationtree->Branch("sps_y",&sps_y,"sps_y/D");
  Sps_Correlationtree->Branch("sps_z",&sps_z,"sps_z/D");
  Sps_Correlationtree->Branch("distance",&distance,"distance/D");
  Sps_Correlationtree->Branch("distance_smallest_nu",&distance_smallest_nu,"distance_smallest_nu/D");//distance between a neutrino correlated track and a spacepoint
  Sps_Correlationtree->Branch("sps_cluster_charge",&sps_cluster_charge,"sps_cluster_charge/D");
  Sps_Correlationtree->Branch("N_Event",&N_Event,"N_Event/I");
  Sps_Correlationtree->Branch("N_Run",&N_Run,"N_Run/I");
  Sps_Correlationtree->Branch("N_SubRun",&N_SubRun,"N_SubRun/I");
  Sps_Correlationtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  Sps_Correlationtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  Sps_Correlationtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  Sps_Correlationtree->Branch("distance_rand_vtx",&distance_rand_vtx,"distance_rand_vtx/D");
  Sps_Correlationtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  Sps_Correlationtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  Sps_Correlationtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  Sps_Correlationtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Sps_Correlationtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Sps_Correlationtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Sps_Correlationtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Sps_Correlationtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Sps_Correlationtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Sps_Correlationtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  Sps_Correlationtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  Sps_Correlationtree->Branch("distance_nu_cosmic_smallest" ,&distance_nu_cosmic_smallest ,"distance_nu_cosmic_smallest/D" );
  Sps_Correlationtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  Sps_Correlationtree->Branch("pointdistance_trk_smallest" ,&pointdistance_trk_smallest ,"pointdistance_trk_smallest/D" );
  Sps_Correlationtree->Branch("distance_trk_smallest" ,&distance_trk_smallest ,"distance_trk_smallest/D" );

  potTree = tfs->make<TTree>("potTree","potTree");
  potTree->Branch("sr_pot", &sr_pot, "sr_pot/D");
  potTree->Branch("sr_run", &sr_run, "sr_run/I");
  potTree->Branch("sr_sub_run", &sr_sub_run, "sr_sub_run/I");
}
void gammacorrelation::beginJob(){

  mudar_beginJob(); 
  gammacorrelation_beginJob();
}

void gammacorrelation::mudar_endSubRun(art::SubRun const &sr){

  art::InputTag MCTproducer("generator");
  _pot=0; 

  art::Handle<sumdata::POTSummary> potSummaryHandle;
  _pot = sr.getByLabel(MCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.;
  std::cout << "[muDARfeasibility::endSubRun] Storing POT info!" << std::endl;

  _run_sr = sr.run();
  _sub_sr = sr.subRun();
  _subrun_tree->Fill();

}

void gammacorrelation::gammacorrelation_endSubRun(art::SubRun const &sr){


  // Note: the entire subrun's POT is recorded in the tree for every event.
  // You must only add it once per subrun to get the correct number.

  art::Handle<sumdata::POTSummary> potsum_h;

  if (!isData) { // MC only (data is dealt with using Zarko's script)

  if(sr.getByLabel("generator", potsum_h)) {

    sr_pot = potsum_h->totpot;
  }
}

sr_run = sr.run();
sr_sub_run = sr.subRun();

potTree->Fill();
}

void gammacorrelation::endJob(){
}

void gammacorrelation::endSubRun(art::SubRun const &sr) {

mudar_endSubRun(sr);
gammacorrelation_endSubRun(sr);
}

DEFINE_ART_MODULE(gammacorrelation)
