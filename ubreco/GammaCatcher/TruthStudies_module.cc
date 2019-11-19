////////////////////////////////////////////////////////////////////////
// Class:       TruthStudies
// Plugin Type: analyzer (art v3_01_02)
// File:        TruthStudies_module.cc
//
// Generated at Thu Oct  3 11:07:53 2019 by Avinay Bhat using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/ProcessHistoryID.h"


//LARSOFT
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "canvas/Persistency/Provenance/BranchType.h"



// ROOT
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TRandom3.h"
#include <memory>


class TruthStudies;


class TruthStudies : public art::EDAnalyzer {
public:
  explicit TruthStudies(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthStudies(TruthStudies const&) = delete;
  TruthStudies(TruthStudies&&) = delete;
  TruthStudies& operator=(TruthStudies const&) = delete;
  TruthStudies& operator=(TruthStudies&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void endSubRun(art::SubRun const& sr) override;

private:

  // Declare member data here.
  std::string fMCProducer_tag;

  Int_t pdg_particle=0;
  Double_t MC_Particle_Energy=0.0;
  Double_t Event_Energy_total=0.0;
  Double_t Event_Energy_5=0.0;
  Double_t Event_Energy_10=0.0;
  Double_t Event_Energy_20=0.0;
  Double_t Event_Energy_30=0.0;
  Double_t Event_Energy_40=0.0;
  Double_t Event_Energy_50=0.0;
  Double_t Event_Energy_100=0.0;
  Double_t position_X=-9999.0;
  Double_t position_Y=-9999.0;
  Double_t position_Z=-9999.0;
  Double_t distance_center=-9999.0;
  Double_t MC_Particle_Energy_MeV=0.0;
  Double_t Event_Energy_total_MeV=0.0;
  Double_t Event_Energy_5_MeV=0.0;
  Double_t Event_Energy_10_MeV=0.0;
  Double_t Event_Energy_20_MeV=0.0;
  Double_t Event_Energy_30_MeV=0.0;
  Double_t Event_Energy_40_MeV=0.0;
  Double_t Event_Energy_50_MeV=0.0;
  Double_t Event_Energy_100_MeV=0.0;
  Int_t StatusCode=-9999;
  Double_t distance_true_reco=999999.0;
  Double_t distance_true_reco_smallest=999999.0;
  Double_t MC_Particle_Energy_MeV_smallest=99999.0;
  Double_t sps_cluster_charge_smallest=99999.0;

  TTree *Eventtree;
  TTree *MCParticletree;




  // Declare member data here.

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  std::string fHit_tag,fpfparticle_tag,fvertex_tag,fsps_tag,fcluster_tag,fReco_track_tag;

  Int_t evttime=0;

  Double32_t sps_x,sps_y,sps_z,sps_hit_charge,sps_cluster_charge,sps_cluster_charge_match,sps_cluster_charge10,sps_cluster_charge20,sps_cluster_charge50;
  Double_t distance, distance_smallest,Event_cluster_charge;
  Double_t Vertex_x,Vertex_y,Vertex_z;
  TRandom3 rand;
  Double_t _rand_vtx_x, _rand_vtx_y, _rand_vtx_z, distance_rand_vtx, distance_smallest_rand_vtx;
  Int_t neutrinos,N_sps,N_Event,N_Run,N_SubRun,N_sps10,N_sps20,N_sps50;
  float _maxTrkLen;
  int   _neutrinoshowers;
  int   _neutrinotracks;
  float _muon_px, _muon_py, _muon_pz;
  Double_t tracklength = 0.0,X_reco=0.0,Y_reco=0.0,Z_reco=0.0;
  Double_t pointdistance_nu=0.0;
  Double_t pointdistance_nu_smallest=0.0;
  Double_t distance_nu_smallest=0.0;


  Double_t pointdistance_trk=0.0;
  Double_t pointdistance_trk_smallest=0.0;
  Double_t distance_trk_smallest=0.0;


  Double_t pointdistance_smallest_nu;
  Double_t track_point_length=0.0;
  Double_t track_point_length_smallest=0.0;
  Double_t distance_smallest_nu;
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


};


TruthStudies::TruthStudies(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.


  fMCProducer_tag = p.get<std::string>("MCProducer_tag");
  fHit_tag = p.get<std::string>("hit_tag"  );
  fpfparticle_tag=p.get<std::string>("pfparticle_tag");
  fvertex_tag=p.get<std::string>("vertex_tag");
  fsps_tag=p.get<std::string>("sps_tag");
  fcluster_tag=p.get<std::string>("cluster_tag");
  fData     = p.get< bool >("IsData");
  fReco_track_tag = p.get<std::string>("recotrack_tag"  );

}


void TruthStudies::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() << " SubRun Number: " << e.subRun() <<"**********************"<< std::endl;

  art::Timestamp ts = e.time();
  //uint64_t evttime = ts.timeHigh();
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
  Event_cluster_charge=0,sps_cluster_charge50=0,sps_cluster_charge20=0,sps_cluster_charge10=0;





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

      // std::cout<<V_xyz.x()<<std::endl;

    }

  }
  distance_smallest=1e10;

  distance_smallest_rand_vtx=1e10;
  N_sps=0,N_sps10=0,N_sps20=0,N_sps50=0;
  Vertex_x=-9999.0;
  Vertex_y=-9999.0;
  Vertex_z=-9999.0;


  // random generated "neutrino" vtx in the fiducial volume
  _rand_vtx_x = rand.Uniform(fidVolMinX,fidVolMaxX);
  _rand_vtx_y = rand.Uniform(fidVolMinY,fidVolMaxY);
  _rand_vtx_z = rand.Uniform(fidVolMinZ,fidVolMaxZ);

  rndvtx = TVector3(_rand_vtx_x,_rand_vtx_y,_rand_vtx_z);
  distance_nu_smallest=1e10;
  cosmic_trk_50=0;

  for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP

    auto const& track = recotrack_handle->at(i_t);




    if ( (sqrt((pow(nuvtx.position().x()-track.Start().X(),2))+(pow(nuvtx.position().y()-track.Start().Y(),2))+ (pow(nuvtx.position().z()-track.Start().Z(),2))) <5.0) ||  (sqrt((pow(nuvtx.position().x()-track.End().X(),2))+(pow(nuvtx.position().y()-track.End().Y(),2))+ (pow(nuvtx.position().z()-track.End().Z(),2)))<5.0))
    continue;

    if (track.Length()<20.0)
    continue;

    pointdistance_nu_smallest=1e10;



    for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP
      X_reco=track.LocationAtPoint(m).X();
      Y_reco=track.LocationAtPoint(m).Y();
      Z_reco=track.LocationAtPoint(m).Z();


      pointdistance_nu=sqrt((pow(nuvtx.position().x()-X_reco,2))+(pow(nuvtx.position().y()-Y_reco,2))+ (pow(nuvtx.position().z()-Z_reco,2)));

      if(pointdistance_nu<pointdistance_nu_smallest){
        pointdistance_nu_smallest=pointdistance_nu;

      }

    }
    // std::cout<<"pointdistance_nu_smallest: "<<pointdistance_nu_smallest<<std::endl;

    if (pointdistance_nu_smallest<50.0){
      cosmic_trk_50++;

    }

    if(pointdistance_nu_smallest<distance_nu_smallest){
      distance_nu_smallest=pointdistance_nu_smallest;
      tracklength= track.Length();

    }
    // std::cout<<"distance_nu_smallest: "<<distance_nu_smallest<<std::endl;
    // std::cout<<"tracklength: "<<tracklength<<std::endl;
    // std::cout<<"cosmic_trk_50: "<<cosmic_trk_50<<std::endl;
  }


  for(size_t s=0;s<spacepoint_handle->size();s++){ //START SPACEPOINT LOOP
    N_sps++;
    // std::cout<<"SpacePoint_Number: "<<s<<std::endl;


    // std::cout<<"SpacePoint_Number: "<<s<<std::endl;

    auto sps = spacepoint_handle->at(s);
    auto cluster_v = sps_clus_assn_v.at(s);


    sps_x=sps.XYZ()[0];
    sps_y=sps.XYZ()[1];
    sps_z=sps.XYZ()[2];




    distance_trk_smallest=1e10;

    distance_smallest_nu=1e10;


    for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP




      auto const& track = recotrack_handle->at(i_t);

      /*


      if ( (sqrt((pow(nuvtx.position().x()-track.Start().X(),2))+(pow(nuvtx.position().y()-track.Start().Y(),2))+ (pow(nuvtx.position().z()-track.Start().Z(),2))) <5.0) ||  (sqrt((pow(nuvtx.position().x()-track.End().X(),2))+(pow(nuvtx.position().y()-track.End().Y(),2))+ (pow(nuvtx.position().z()-track.End().Z(),2)))<5.0)){

        // std::cout<<"************************Nu Track # "<<i_t<<" *************************"<<std::endl;
        // std::cout<<"************************Tracklength # "<<track.Length()<<" *************************"<<std::endl;


        pointdistance_smallest_nu=1e10;

        track_point_length=0.0;
        for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP

          X_reco_nu=track.LocationAtPoint(m).X();
          Y_reco_nu=track.LocationAtPoint(m).Y();
          Z_reco_nu=track.LocationAtPoint(m).Z();



          pointdistance_nu=sqrt((pow(sps_x-X_reco_nu,2))+(pow(sps_y-Y_reco_nu,2))+(pow(sps_z-Z_reco_nu,2)));
          // std::cout<<"pointdistance_nu: "<<pointdistance_nu<<std::endl;
          // std::cout<<"trackpointlength: "<<sqrt((pow(X_reco_nu-nuvtx.position().x(),2))+(pow(Y_reco_nu-nuvtx.position().y(),2))+ (pow(Z_reco_nu-nuvtx.position().z(),2)))<<std::endl;

          if(pointdistance_nu<pointdistance_smallest_nu){//comparison IF loop

            pointdistance_smallest_nu=pointdistance_nu;

            X_reco_smallest_nu=X_reco_nu;
            Y_reco_smallest_nu=Y_reco_nu;
            Z_reco_smallest_nu=Z_reco_nu;
            track_point_length= sqrt((pow(X_reco_nu-nuvtx.position().x(),2))+(pow(Y_reco_nu-nuvtx.position().y(),2))+ (pow(Z_reco_nu-nuvtx.position().z(),2)));

          }

        }//END RECO POINT LOOP

        // std::cout<<"pointdistance_smallest_nu: "<<pointdistance_smallest_nu<<std::endl;
        // std::cout<<"track_point_length: "<<track_point_length<<std::endl;



        if(pointdistance_smallest_nu<distance_smallest_nu){
          distance_smallest_nu=pointdistance_smallest_nu;
          track_point_length_smallest=track_point_length;
          X_reco_best_nu=X_reco_smallest_nu; //variables for the coordinates of the nearest reco track
          Z_reco_best_nu=Z_reco_smallest_nu; //variables for the coordinates of the nearest reco track
        }

        // std::cout<<"distance_smallest_nu: "<<distance_smallest_nu<<std::endl;
        // std::cout<<"track_point_length_smallest: "<<track_point_length_smallest<<std::endl;





      }//END IF LOOP

      */

      // std::cout<<"Track Number: "<<i_t<<std::endl;

      pointdistance_trk_smallest=1e10;

      for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP

        pointdistance_trk= sqrt((pow(track.Start().X()-sps_x,2))+(pow(track.Start().Y()-sps_y,2))+ (pow(track.Start().Z()-sps_z,2)));

        if (pointdistance_trk<pointdistance_trk_smallest){
          pointdistance_trk_smallest=pointdistance_trk;
        }

      }

      // std::cout<<"Distance to track: "<<pointdistance_trk_smallest<<std::endl;

      if (pointdistance_trk_smallest<distance_trk_smallest){
        distance_trk_smallest=pointdistance_trk_smallest;

      }

    }

    // std::cout<<"pointdistance_smallest_nu_here: "<<pointdistance_smallest_nu<<std::endl;
    // std::cout<<"distance_smallest_nu_here: "<<distance_smallest_nu<<std::endl;
    // std::cout<<"track_point_length_here: "<<track_point_length<<std::endl;
    // std::cout<<"track_point_length_smallest: "<<track_point_length_smallest<<std::endl;


    //
    // if ((distance_smallest_nu<((R/H)*(track_point_length_smallest)) && (track_point_length_smallest<H))||((distance_smallest_nu>R) && (track_point_length_smallest>H)))
    // continue;







    // if (neutrinos==0){
    //   Sps_Correlationtree->Fill();
    //   continue;
    // }

    Vertex_x=nuvtx.position().x();
    Vertex_y=nuvtx.position().y();
    Vertex_z=nuvtx.position().z();

    // std::cout<<"Vertex X: "<<nuvtx.position().x()<<"****Vertex Y: "<<nuvtx.position().y()<<"****Vertex Z: "<<nuvtx.position().z()<<std::endl;
    //
    // std::cout<<"sps_x: "<<sps_x<<"****sps_y: "<<sps_y<<"****sps_z: "<<sps_z<<std::endl;
    distance= sqrt((pow(nuvtx.position().x()-sps_x,2))+(pow(nuvtx.position().y()-sps_y,2))+ (pow(nuvtx.position().z()-sps_z,2)));
    // std::cout<<"Distance: "<<distance<<std::endl;
    sps_cluster_charge =0;

    for (auto const& cluster: cluster_v){
      auto plane = cluster->View();

      if (plane==2){
        sps_cluster_charge = cluster->Integral();
         // std::cout<<"sps_cluster_charge: "<<sps_cluster_charge<<std::endl;
      }
    }

    if (distance<distance_smallest){
      distance_smallest=distance;
      for (auto const& cluster: cluster_v){
        auto plane = cluster->View();

        if (plane==2){
          Event_cluster_charge = cluster->Integral();
        }
      }

    }

    // std::cout<<"distance: "<<distance<<std::endl;
    // std::cout<<"distance_smallest: "<<distance_smallest<<std::endl;


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


    distance_rand_vtx= sqrt((pow(_rand_vtx_x-sps_x,2))+(pow(_rand_vtx_y-sps_y,2))+ (pow(_rand_vtx_z-sps_z,2)));
    if (distance_rand_vtx<distance_smallest_rand_vtx){
      distance_smallest_rand_vtx=distance_rand_vtx;
    }
    // std::cout<<"distance_rand_vtx: "<<distance_rand_vtx<<std::endl;


    // std::cout<<"sps_cluster_charge: "<<sps_cluster_charge<<std::endl;
    Sps_Correlationtree->Fill();
  }//END SPACEPOINT LOOP


  // std::cout<<"distance_smallest: "<<distance_smallest<<std::endl;
  // std::cout<<"Event_cluster_charge: "<<Event_cluster_charge<<std::endl;









  auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCProducer_tag);

  Event_Energy_5=0.0;
  Event_Energy_10=0.0;
  Event_Energy_20=0.0;
  Event_Energy_30=0.0;
  Event_Energy_40=0.0;
  Event_Energy_50=0.0;
  Event_Energy_100=0.0;
  Event_Energy_total=0.0;

  for (size_t p = 0; p < mcp_h->size(); p++){//START MC PARTICLE FOR LOOP
    auto mcp = mcp_h->at(p);

    if (mcp.PdgCode() != 22)
    continue;

    if (mcp.Process() != "primary")
    continue;


    StatusCode=mcp.StatusCode();
    position_X=-9999.0;
    position_Y=-9999.0;
    position_Z=-9999.0;
    pdg_particle=-9999;
    MC_Particle_Energy= -9999.0;

    pdg_particle=mcp.PdgCode();


    position_X =mcp.Vx();
    position_Y =mcp.Vy();
    position_Z =mcp.Vz();

    MC_Particle_Energy= mcp.E();

    distance_center=sqrt((pow((position_X-0),2))+(pow((position_Y-0),2))+ (pow((position_Z-500),2)));


    // std::cout<<"MC Particle # "<<p<<std::endl;
    // std::cout<<"StatusCode "<<StatusCode<<std::endl;
    // std::cout<<"MCP Process "<<mcp.Process()<<std::endl;
    std::cout<<"MC Particle Energy# "<<MC_Particle_Energy<<std::endl;
    std::cout<<"pdg_particle "<<pdg_particle<<std::endl;


    /*
    Since we're using a single particle gun to generate photons at X=0,Y=0,Z=500, we want to calculate the energy deposition within 10,20,30,40,50 cm (3D distance) of these
    coordinates and also the total energy within the event.
    */


    Event_Energy_total += MC_Particle_Energy;
    std::cout<<"MC Particle Energy summed up so far:  "<<Event_Energy_total<<std::endl;

    if (distance_center <5.0){
      Event_Energy_5 += MC_Particle_Energy;

    }

    if (distance_center <10.0){
      Event_Energy_10 += MC_Particle_Energy;

    }

    if (distance_center <20.0){
      Event_Energy_20 += MC_Particle_Energy;

    }

    if (distance_center <30.0){
      Event_Energy_30 += MC_Particle_Energy;

    }

    if (distance_center <40.0){
      Event_Energy_40 += MC_Particle_Energy;

    }

    if (distance_center <50.0){
      Event_Energy_50 += MC_Particle_Energy;

    }

    if (distance_center <100.0){
      Event_Energy_100 += MC_Particle_Energy;

    }


    MC_Particle_Energy_MeV=MC_Particle_Energy*1000;


    MCParticletree->Fill();

    Event_Energy_total_MeV=Event_Energy_total*1000;
    Event_Energy_5_MeV=Event_Energy_5*1000;
    Event_Energy_10_MeV=Event_Energy_10*1000;
    Event_Energy_20_MeV=Event_Energy_20*1000;
    Event_Energy_30_MeV=Event_Energy_30*1000;
    Event_Energy_40_MeV=Event_Energy_40*1000;
    Event_Energy_50_MeV=Event_Energy_50*1000;
    Event_Energy_100_MeV=Event_Energy_100*1000;





  }//END MC PARTICLE LOOP





  // std::cout<<"MC Particle Energy per Event# "<<MC_Particle_Energy<<std::endl;
  std::cout<<"MC Particle Energy MeV per Event# "<<Event_Energy_total_MeV<<std::endl;


  Eventtree->Fill();




  distance_true_reco=999999.0;

  MC_Particle_Energy_MeV_smallest=99999.0;
  sps_cluster_charge_smallest=99999.0;

  for(size_t s=0;s<spacepoint_handle->size();s++){ //START SPACEPOINT LOOP
    auto sps = spacepoint_handle->at(s);
    auto cluster_v = sps_clus_assn_v.at(s);

    std::cout<<"SPS # "<<s<<std::endl;


    for (auto const& cluster: cluster_v){
      auto plane = cluster->View();

      if (plane==2){
        sps_cluster_charge_match = cluster->Integral();
         // std::cout<<"sps_cluster_charge: "<<sps_cluster_charge<<std::endl;
      }
    }




    distance_true_reco_smallest=999999.0;
    // std::cout<<"SpacePoint_Number: "<<s<<std::endl;
    // std::cout<<"sps_cluster_charge: "<<sps_cluster_charge_match<<std::endl;

    for (size_t p = 0; p < mcp_h->size(); p++){//START MC PARTICLE FOR LOOP

      // std::cout<<"MC Particle # "<<p<<std::endl;

      auto mcp = mcp_h->at(p);

      // // std::cout<<"MC Particle Energy # "<<mcp.E()*1000<<std::endl;
      if (mcp.PdgCode() != 11 && mcp.PdgCode()!= -11)
      continue;
      // //
      // // if (mcp.Process() != "primary")
      // // continue;
      //
      //
      // if (mcp.PdgCode() != 11)
      // continue;
      //
      // if (mcp.PdgCode() != -11)
      // continue;

      distance_true_reco=sqrt((pow((mcp.Vx()-sps.XYZ()[0]),2))+(pow((mcp.Vy()-sps.XYZ()[1]),2))+ (pow((mcp.Vz()-sps.XYZ()[2]),2)));
      // std::cout<<"distance_true_reco: "<<distance_true_reco<<std::endl;

      if (distance_true_reco < distance_true_reco_smallest){
        distance_true_reco_smallest=distance_true_reco;
        MC_Particle_Energy_MeV_smallest=mcp.E()*1000;
        sps_cluster_charge_smallest=sps_cluster_charge_match;

      }
    }

      std::cout<<"distance_true_reco_smallest: "<<distance_true_reco_smallest<<std::endl;
      std::cout<<"sps_cluster_charge_smallest: "<<sps_cluster_charge_smallest<<std::endl;
      std::cout<<"MC Particle Energy_smallest_MeV # "<<MC_Particle_Energy_MeV_smallest<<std::endl;

  }





}//END OF EVENT LOOP

void TruthStudies::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;

  Eventtree = tfs->make<TTree>("Eventtree",    "Eventtree");
  MCParticletree = tfs->make<TTree>("MCParticletree",    "MCParticletree");

  Eventtree->Branch("Event_Energy_total_MeV",&Event_Energy_total,"Event_Energy_total_MeV/D");
  Eventtree->Branch("Event_Energy_5_MeV",&Event_Energy_5_MeV,"Event_Energy_5_MeV/D");
  Eventtree->Branch("Event_Energy_10_MeV",&Event_Energy_10_MeV,"Event_Energy_10_MeV/D");
  Eventtree->Branch("Event_Energy_20_MeV",&Event_Energy_20_MeV,"Event_Energy_20_MeV/D");
  Eventtree->Branch("Event_Energy_30_MeV",&Event_Energy_30_MeV,"Event_Energy_30_MeV/D");
  Eventtree->Branch("Event_Energy_40_MeV",&Event_Energy_40_MeV,"Event_Energy_40_MeV/D");
  Eventtree->Branch("Event_Energy_50_MeV",&Event_Energy_50_MeV,"Event_Energy_50_MeV/D");
  Eventtree->Branch("Event_Energy_100_MeV",&Event_Energy_100_MeV,"Event_Energy_100_MeV/D");
  Eventtree->Branch("evttime",&evttime,"evttime/I");
  Eventtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  Eventtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  Eventtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  Eventtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  Eventtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  Eventtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  Eventtree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
  Eventtree->Branch("N_sps",&N_sps,"N_sps/I");
  Eventtree->Branch("Event_cluster_charge",&Event_cluster_charge,"Event_cluster_charge/D");
  Eventtree->Branch("distance_smallest_rand_vtx",&distance_smallest_rand_vtx,"distance_smallest_rand_vtx/D");
  Eventtree->Branch("N_sps10",&N_sps10,"N_sps10/I");
  Eventtree->Branch("N_sps20",&N_sps20,"N_sps20/I");
  Eventtree->Branch("N_sps50",&N_sps50,"N_sps50/I");
  Eventtree->Branch("sps_cluster_charge10",&sps_cluster_charge10,"sps_cluster_charge10/D");
  Eventtree->Branch("sps_cluster_charge20",&sps_cluster_charge20,"sps_cluster_charge20/D");
  Eventtree->Branch("sps_cluster_charge50",&sps_cluster_charge50,"sps_cluster_charge50/D");
  Eventtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Eventtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Eventtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Eventtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Eventtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Eventtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Eventtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  Eventtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  Eventtree->Branch("distance_nu_smallest" ,&distance_nu_smallest ,"distance_nu_smallest/D" );
  Eventtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  Eventtree->Branch("N_Run" ,&N_Run ,"N_Run/I" );
  Eventtree->Branch("N_SubRun" ,&N_SubRun ,"N_SubRun/I" );



  MCParticletree->Branch("MC_Particle_Energy",&MC_Particle_Energy,"MC_Particle_Energy/D");
  MCParticletree->Branch("position_X",&position_X,"position_X/D");
  MCParticletree->Branch("position_Y",&position_Y,"position_Y/D");
  MCParticletree->Branch("position_Z",&position_Z,"position_Z/D");
  MCParticletree->Branch("pdg_particle",&pdg_particle,"pdg_particle/D");


  Sps_Correlationtree = tfs->make<TTree>("Sps_Correlationtree",    "Sps_Correlationtree");
  Sps_Correlationtree->Branch("evttime",&evttime,"evttime/I");
  Sps_Correlationtree->Branch("sps_x",&sps_x,"sps_x/D");
  Sps_Correlationtree->Branch("sps_y",&sps_y,"sps_y/D");
  Sps_Correlationtree->Branch("sps_z",&sps_z,"sps_z/D");
  Sps_Correlationtree->Branch("distance",&distance,"distance/D");
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
  Sps_Correlationtree->Branch("distance_nu_smallest" ,&distance_nu_smallest ,"distance_nu_smallest/D" );
  Sps_Correlationtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  Sps_Correlationtree->Branch("pointdistance_trk_smallest" ,&pointdistance_trk_smallest ,"pointdistance_trk_smallest/D" );
  Sps_Correlationtree->Branch("distance_trk_smallest" ,&distance_trk_smallest ,"distance_trk_smallest/D" );



  potTree = tfs->make<TTree>("potTree","potTree");
  potTree->Branch("sr_pot", &sr_pot, "sr_pot/D");
  potTree->Branch("sr_run", &sr_run, "sr_run/I");
  potTree->Branch("sr_sub_run", &sr_sub_run, "sr_sub_run/I");



}

void TruthStudies::endJob()
{
  // Implementation of optional member function here.
}

void TruthStudies::endSubRun(art::SubRun const &sr) {
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


DEFINE_ART_MODULE(TruthStudies)
