////////////////////////////////////////////////////////////////////////
// Class:       gammacorrelation
// Plugin Type: analyzer (art v3_01_02)
// File:        gammacorrelation_module.cc
//
// Generated at Wed Mar 20 10:42:52 2019 by Avinay Bhat using cetskelgen
// from cetlib version v3_05_01.
// Avinay Bhat (avbhat@syr.edu)
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


#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"



#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


#include "TTree.h"
#include "TRandom3.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <memory>

class gammacorrelation;


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

private:

  // Declare member data here.

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  std::string fHit_tag,fpfparticle_tag,fvertex_tag,fsps_tag,fcluster_tag,fReco_track_tag;

  Double32_t sps_x,sps_y,sps_z,sps_hit_charge,sps_cluster_charge,sps_cluster_charge10,sps_cluster_charge20,sps_cluster_charge50;
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


  // art::InputTag fMCTproducer;

};


gammacorrelation::gammacorrelation(fhicl::ParameterSet const& p)
: EDAnalyzer{p}  // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  //All tags get filled in the run_gammacorrelation.fcl file


  fHit_tag = p.get<std::string>("hit_tag"  );
  fpfparticle_tag=p.get<std::string>("pfparticle_tag");
  fvertex_tag=p.get<std::string>("vertex_tag");
  fsps_tag=p.get<std::string>("sps_tag");
  fcluster_tag=p.get<std::string>("cluster_tag");
  fData     = p.get< bool >("IsData");
  fReco_track_tag = p.get<std::string>("recotrack_tag"  );
}

void gammacorrelation::analyze(art::Event const& e)
{


  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() << " SubRun Number: " << e.subRun() <<"**********************"<< std::endl;

  // create list of tracks and showers associated to neutrino candidate
  std::vector<recob::Track  > sliceTracks;
  std::vector<recob::Shower > sliceShowers;

  N_Event=e.event();
  N_Run=e.run();
  N_SubRun=e.subRun();
  // isData = e.isRealData();
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
  for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP{

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
     std::cout<<"pointdistance_nu_smallest: "<<pointdistance_nu_smallest<<std::endl;

    if (pointdistance_nu_smallest<50.0){
      cosmic_trk_50++;

    }

    if(pointdistance_nu_smallest<distance_nu_smallest){
      distance_nu_smallest=pointdistance_nu_smallest;
      tracklength= track.Length();

    }
    std::cout<<"distance_nu_smallest: "<<distance_nu_smallest<<std::endl;
    std::cout<<"tracklength: "<<tracklength<<std::endl;
    std::cout<<"cosmic_trk_50: "<<cosmic_trk_50<<std::endl;
  }


  for(size_t s=0;s<spacepoint_handle->size();s++){ //START SPACEPOINT LOOP
    N_sps++;

    auto sps = spacepoint_handle->at(s);
    auto cluster_v = sps_clus_assn_v.at(s);


    sps_x=sps.XYZ()[0];
    sps_y=sps.XYZ()[1];
    sps_z=sps.XYZ()[2];

    distance_rand_vtx= sqrt((pow(_rand_vtx_x-sps_x,2))+(pow(_rand_vtx_y-sps_y,2))+ (pow(_rand_vtx_z-sps_z,2)));
    if (distance_rand_vtx<distance_smallest_rand_vtx){
      distance_smallest_rand_vtx=distance_rand_vtx;
    }
    // std::cout<<"distance_rand_vtx: "<<distance_rand_vtx<<std::endl;



    if (neutrinos==0){
      Sps_Correlationtree->Fill();
      continue;
    }

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

    // std::cout<<"sps_cluster_charge: "<<sps_cluster_charge<<std::endl;
    Sps_Correlationtree->Fill();
  }//END SPACEPOINT LOOP


  // std::cout<<"distance_smallest: "<<distance_smallest<<std::endl;
  // std::cout<<"Event_cluster_charge: "<<Event_cluster_charge<<std::endl;
  Event_Correlationtree->Fill();//Filled per event
}

void gammacorrelation::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  Event_Correlationtree = tfs->make<TTree>("Event_Correlationtree",    "Event_Correlationtree");
  Event_Correlationtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  Event_Correlationtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  Event_Correlationtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  Event_Correlationtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  Event_Correlationtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  Event_Correlationtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  Event_Correlationtree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
  Event_Correlationtree->Branch("N_sps",&N_sps,"N_sps/I");
  Event_Correlationtree->Branch("Event_cluster_charge",&Event_cluster_charge,"Event_cluster_charge/D");
  Event_Correlationtree->Branch("distance_smallest_rand_vtx",&distance_smallest_rand_vtx,"distance_smallest_rand_vtx/D");
  Event_Correlationtree->Branch("N_sps10",&N_sps10,"N_sps10/I");
  Event_Correlationtree->Branch("N_sps20",&N_sps20,"N_sps20/I");
  Event_Correlationtree->Branch("N_sps50",&N_sps50,"N_sps50/I");
  Event_Correlationtree->Branch("sps_cluster_charge10",&sps_cluster_charge10,"sps_cluster_charge10/D");
  Event_Correlationtree->Branch("sps_cluster_charge20",&sps_cluster_charge20,"sps_cluster_charge20/D");
  Event_Correlationtree->Branch("sps_cluster_charge50",&sps_cluster_charge50,"sps_cluster_charge50/D");
  Event_Correlationtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Event_Correlationtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Event_Correlationtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Event_Correlationtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Event_Correlationtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Event_Correlationtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Event_Correlationtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  Event_Correlationtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  Event_Correlationtree->Branch("distance_nu_smallest" ,&distance_nu_smallest ,"distance_nu_smallest/D" );
  Event_Correlationtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );

  Sps_Correlationtree = tfs->make<TTree>("Sps_Correlationtree",    "Sps_Correlationtree");
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

  potTree = tfs->make<TTree>("potTree","potTree");
  potTree->Branch("sr_pot", &sr_pot, "sr_pot/D");
  potTree->Branch("sr_run", &sr_run, "sr_run/I");
  potTree->Branch("sr_sub_run", &sr_sub_run, "sr_sub_run/I");
}

void gammacorrelation::endJob(){}

void gammacorrelation::endSubRun(art::SubRun const &sr) {
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


DEFINE_ART_MODULE(gammacorrelation)
