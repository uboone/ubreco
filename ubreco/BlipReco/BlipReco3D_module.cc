////////////////////////////////////////////////////////////////////////
// Class:       BlipReco3D
// Plugin Type: producer (art v2_11_03)
// File:        BlipReco3D_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"


#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <memory>

class BlipReco3D;


class BlipReco3D : public art::EDProducer {
public:
  explicit BlipReco3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BlipReco3D(BlipReco3D const &) = delete;
  BlipReco3D(BlipReco3D &&) = delete;
  BlipReco3D & operator = (BlipReco3D const &) = delete;
  BlipReco3D & operator = (BlipReco3D &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // configurable inputs to reconstruction
  float f2DcutY;
  float f2DcutUV;
  float fdeltaY;




  std::string fMC_track_tag,fReco_track_tag,fCluster_tag,fHit_tag,fpfparticle_tag,fVertex_tag;
  // std::string fAssnProducer;
  // std::string fAssn_tag;

  UInt_t  Y_Wire;


  Double_t YV_iou;
  Double_t YU_iou;



  Double_t X_reco=0.0;
  Double_t Y_reco=0.0;
  Double_t Z_reco=0.0;

  //Double_t X_reco_nu=0.0;
  //Double_t Y_reco_nu=0.0;
  //Double_t Z_reco_nu=0.0;

  //Double_t X_reco3d=0.0;
  //Double_t Y_reco3d=0.0;
  //Double_t Z_reco3d=0.0;


  Double_t cluster_hit_z=0.0;
  Double_t cluster_hit_x=0.0;
  //Double_t cluster_hit_z_nu=0.0;
  //Double_t cluster_hit_x_nu=0.0;

  // Double_t Y_cluster_hit_z=0.0;
  // Double_t Y_cluster_hit_x=0.0;

  Double_t Y_cluster_3d_hit_x;
  Double_t Y_cluster_3d_hit_z;

  Double_t Y_cluster_3d_hit_y;


  Double_t X_reco_best=0.0;
  //Double_t Y_reco_best=0.0;
  Double_t Z_reco_best=0.0;

  //Double_t X_reco_best_nu=0.0;
  //Double_t Y_reco_best_nu=0.0;
  //Double_t Z_reco_best_nu=0.0;

  //Double_t track_point_length=0.0;
  //Double_t track_point_length_smallest=0.0;

  Double_t X_reco_best3d=0.0;
  Double_t Y_reco_best3d=0.0;
  Double_t Z_reco_best3d=0.0;





  Double_t pointdistance=0;
  //Double_t pointdistance_nu=0;

  //Double_t pointdistance3d=0;


  //Double_t V_cluster_z;
  //Double_t U_cluster_z;
  Double_t cluster_y_new;
  Double_t deltaY;



  Int_t plane;

  Double_t X_reco_smallest=0;
  Double_t Z_reco_smallest=0;

  //Double_t X_reco_smallest_nu=0;
  //Double_t Z_reco_smallest_nu=0;

  //Double_t X_reco_smallest3dV=0;
  //Double_t Z_reco_smallest3dV=0;
  //Double_t Y_reco_smallest3dV=0;

  //Double_t X_reco_smallest3dU=0;
  //Double_t Z_reco_smallest3dU=0;
  //Double_t Y_reco_smallest3dU=0;

  //Double_t X_reco_smallest3d=0;
  //Double_t Z_reco_smallest3d=0;
  //Double_t Y_reco_smallest3d=0;

  Double_t pointdistance_smallest;
  //Double_t pointdistance_smallest_nu;

  //Double_t pointdistance_smallestV;
  //Double_t pointdistance_smallestU;

  Double_t distance_smallest3d;
  //Double_t distance_smallestV;
  //Double_t distance_smallestU;
  Double_t distance_smallest;
  Double_t distance_smallest_nu;
  Double_t Y_charge;
  Double_t V_charge;
  Double_t U_charge;



  Double_t Y_cluster_charge;
  // Double_t Y_cluster_energy;

  Double_t V_cluster_charge;
  // Double_t V_cluster_energy;

  Double_t U_cluster_charge;
  // Double_t U_cluster_energy;


  Double_t deltaY_smallest;

  Double_t wire2cm,time2cm;


  Double_t start_tick_cluster2,start_tick_cluster1,start_tick_cluster0,end_tick_cluster2,end_tick_cluster1,end_tick_cluster0;

  std::vector<Double_t> Start_Cluster0;
  std::vector<Double_t> End_Cluster0;
  std::vector<Double_t> Start_Cluster1;
  std::vector<Double_t> End_Cluster1;
  std::vector<Double_t> Start_Cluster2;
  std::vector<Double_t> End_Cluster2;
  std::vector<Double_t> Y_clus_hitsize;
  std::vector<Double_t> V_clus_hitsize;
  std::vector<Double_t> U_clus_hitsize;
  std::vector<Int_t> Y_index_vector;
  std::vector<Int_t> V_index_vector;
  std::vector<Int_t> U_index_vector;
  std::vector<Int_t> Y_U_index;
  std::vector<Int_t> Y_V_index;






  std::vector<Double_t> IoU_U;
  std::vector<Double_t> IoU_V;

  std::vector<Double_t> y_U;
  std::vector<Double_t> y_V;

  Double_t V_t_min_abs; // smallest start point of the 2
  Double_t V_t_max_abs;    // largest start point of the three

  Double_t V_t_min_common; //  start point of the overlap
  Double_t V_t_max_common; // end point of the overlap

  Double_t U_t_min_abs; // smallest start point of the 2
  Double_t U_t_max_abs;    // largest start point of the three

  Double_t U_t_min_common; //  start point of the overlap
  Double_t U_t_max_common; // end point of the overlap


  Double_t  V_iou;
  Double_t V_biggest_iou;
  Int_t V_match_multiplicity;

  Double_t  U_iou;
  Double_t U_biggest_iou;
  Int_t U_match_multiplicity;

  Double_t Y_clus_lifetime;
  Double_t Y_clus_hitSize;
  Double_t V_clus_hitSize;
  Double_t U_clus_hitSize;


  TTree *Clustertree;

  TTree *Matchingtree;




};


BlipReco3D::BlipReco3D(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::SpacePoint > >();
  produces< art::Assns <recob::Cluster, recob::SpacePoint> >();
  // produces< art::Assns <recob::Spacepoint, recob::Hit> >();


  f2DcutY  = p.get<float>("TwoDcutY");
  f2DcutUV = p.get<float>("TwoDcutUV");
  fdeltaY  = p.get<float>("deltaYcut");

  fMC_track_tag   = p.get<std::string>("mctrack_tag"    );
  fReco_track_tag = p.get<std::string>("recotrack_tag"  );
  fCluster_tag = p.get<std::string>("cluster_tag"  );
  fHit_tag = p.get<std::string>("hit_tag"  );
  fpfparticle_tag=p.get<std::string>("pfparticle_tag");
  fVertex_tag=p.get<std::string>("vertex_tag");

}

void BlipReco3D::produce(art::Event & e)//START EVENT LOOP
{

  std::unique_ptr< std::vector< recob::SpacePoint> > SpacePoint_v(new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assns <recob::Cluster, recob::SpacePoint> > sps_clus_assn_v(new art::Assns<recob::Cluster,recob::SpacePoint> );
  // std::unique_ptr< art::Assns <recob::Hit, recob::SpacePoint> > sps_hit_assn_v(new art::Assns<recob::Hit,recob::SpacePoint> );

  // cluster pointer maker which we will use to retrieve a cluster art::Ptr given its absolute event index
  //art::PtrMaker<recob::Cluster> ClusPtrMaker(e);

  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() <<"**********************"<< std::endl;

  // Implementation of required member function here.
  art::Handle<std::vector<recob::Cluster> > cluster_handle;
  e.getByLabel(fCluster_tag,cluster_handle);

  art::Handle<std::vector<recob::Track> > recotrack_handle;
  e.getByLabel(fReco_track_tag,recotrack_handle);

  art::Handle<std::vector<recob::Hit> > hit_handle;
  e.getByLabel(fHit_tag,hit_handle);

  art::FindMany<recob::Hit> clus_hit_assn_v(cluster_handle, e, fCluster_tag);


  // art::FindMany<recob::Hit> sps_hit_assn_v(cluster_handle, e, fCluster_tag);

  art::Handle<std::vector<recob::PFParticle> > pfparticle_handle;
  e.getByLabel(fpfparticle_tag,pfparticle_handle);

  art::Handle<std::vector<recob::Vertex> > vertex_handle;
  e.getByLabel(fVertex_tag,vertex_handle);

  art::FindMany<recob::Vertex> pfp_vertex_assn_v(pfparticle_handle, e, fVertex_tag);

  auto const& cluster(*cluster_handle);

  // load neutrino vertex
  //   loop through tracks and identify those that are associated to the vertex (i.e. come very close to it)
  // if so, flag that track index as "neutrino track"


  // Nu_track_index.clear();




  recob::Vertex nuvtx;
  size_t neutrinos = 0;

  for (size_t p=0; p < pfparticle_handle->size(); p++) {
    auto pfp = pfparticle_handle->at(p);

    if (pfp.IsPrimary() == false) continue;

    auto PDG = fabs(pfp.PdgCode());
    if ( (PDG == 12) || (PDG == 14) ) {

      neutrinos += 1;

      // auto pfpkey = p;
      auto ass_vtx_v  =pfp_vertex_assn_v.at( p );
      if (ass_vtx_v.size() != 1)
      std::cout << "ERROR. Neutrino not associated with a single vertex..." << std::endl;
      nuvtx = *(ass_vtx_v.at(0));
      // V_xyz = nuvtx.position();

      // std::cout<<V_xyz.x()<<std::endl;

    }

  }





  // std::cout<<"size_cluster: "<<cluster_handle->size()<<std::endl;
  // std::cout<<"size_track: "<<recotrack_handle->size()<<std::endl;

  for (size_t i_c = 0, size_cluster = cluster_handle->size(); i_c != size_cluster; ++i_c) { //start cluster FOR loop for calculating 2-D distance

    // std::cout<<"Cluster # "<<i_c<<std::endl;
    //
    // COMMENT OUT FROM HERE TO THE END OF IF CONTINUE STATEMENT BEFORE THE RECO TRACK FOR LOOP IF RUNNING OVER DATA OR MONTE CARLO WITH COSMICS
    //
    // UNCOMMENT IF RUNNING OVER DATA OR MONTE CARLO WITH NO COSMICS
    //
    //
    //
    // auto hits = clus_hit_assn_v.at(i_c);
    //
    //
    //
    // for (auto const& hit : hits) {//START CLUSTER HIT LOOP
    //
    //   cluster_hit_z = hit->WireID().Wire * wire2cm;//Also equal to Cluster_hit_wire_cm
    //   cluster_hit_x = (hit->PeakTime() * time2cm)-44.575 ;//Also equal to Cluster_hit_time_cm
    //
    // }
    //
    // // std::cout<<"Number of reco tracks: "<<recotrack_handle->size()<<std::endl;
    //
    //
    //
    // if(cluster[i_c].View()==2){
    //   // std::cout<<"HERE 2 !"<<std::endl;
    //   plane=2;
    //   Start_Cluster2.push_back((cluster[i_c].StartTick ())-3.0);//added +- 3.0 time tick tolerances
    //   End_Cluster2.push_back((cluster[i_c].EndTick ())+3.0);
    //   // Start_Cluster2.push_back((cluster[i_c].StartTick ()));//added +- 3.0 time tick tolerances
    //   // End_Cluster2.push_back((cluster[i_c].EndTick ()));
    //   Y_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
    //   Y_index_vector.push_back(i_c); //Y Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
    // }
    //
    // if(cluster[i_c].View()==1){//IF LOOP TO CHECK WHAT PLANE A CLUSTER BELONGS TO
    //   // std::cout<<"HERE 1 !"<<std::endl;
    //   plane=1;
    //   Start_Cluster1.push_back((cluster[i_c].StartTick ())-3.0);
    //   End_Cluster1.push_back((cluster[i_c].EndTick ())+3.0);
    //   V_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
    //   // Start_Cluster1.push_back((cluster[i_c].StartTick ()));
    //   // End_Cluster1.push_back((cluster[i_c].EndTick ()));
    //   V_index_vector.push_back(i_c);//V Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
    // }
    //
    // if(cluster[i_c].View()==0){//IF LOOP TO CHECK WHAT PLANE A CLUSTER BELONGS TO
    //   // std::cout<<"HERE 0 !"<<std::endl;
    //   plane=0;
    //   //
    //   Start_Cluster0.push_back((cluster[i_c].StartTick ())-3.0);
    //   End_Cluster0.push_back((cluster[i_c].EndTick ())+3.0);
    //   U_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
    //   // Start_Cluster0.push_back((cluster[i_c].StartTick ()));
    //   // End_Cluster0.push_back((cluster[i_c].EndTick ()));
    //   U_index_vector.push_back(i_c);//U Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
    // }
    //
    //
    //
    //
    // distance_smallest=1e10; //Variable for 2D distance between a cluster and nearest reco track, initialized to a large number for comparison
    // distance_smallest_nu=1e10;
    //
    //
    //
    //
    //
    // Clustertree->Fill();
    //
    // if (recotrack_handle->size()==0)
    //
    // continue;



    distance_smallest=1e10; //Variable for 2D distance between a cluster and nearest reco track, initialized to a large number for comparison
    distance_smallest_nu=1e10;




    for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP

      auto hits = clus_hit_assn_v.at(i_c);
      auto const& track = recotrack_handle->at(i_t);

      //Neutrino Correlated Tracks are skipped. So all activity close to a neutrino track will be reconstructed.

      if ( (sqrt((pow(nuvtx.position().x()-track.Start().X(),2))+(pow(nuvtx.position().y()-track.Start().Y(),2))+ (pow(nuvtx.position().z()-track.Start().Z(),2))) <5.0) ||  (sqrt((pow(nuvtx.position().x()-track.End().X(),2))+(pow(nuvtx.position().y()-track.End().Y(),2))+ (pow(nuvtx.position().z()-track.End().Z(),2)))<5.0))
      continue;


      pointdistance_smallest=1e10;////Variable for distance between a cluster hit and a reco track point, initialized to a large number for comparison

      for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP


        X_reco=track.LocationAtPoint(m).X();
        Y_reco=track.LocationAtPoint(m).Y();
        Z_reco=track.LocationAtPoint(m).Z();


        auto const* geom = ::lar::providerFrom<geo::Geometry>();
        auto V_wire_cm = geom->WireCoordinate(Y_reco,Z_reco,geo::PlaneID(0,0,1)) * wire2cm;

        auto V_time_cm = X_reco;

        auto U_wire_cm = geom->WireCoordinate(Y_reco,Z_reco,geo::PlaneID(0,0,0)) * wire2cm;

        auto U_time_cm = X_reco;


        for (auto const& hit : hits) {//START CLUSTER HIT LOOP

          cluster_hit_z = hit->WireID().Wire * wire2cm;//Also equal to Cluster_hit_wire_cm
          cluster_hit_x = (hit->PeakTime() * time2cm)-44.575 ;//Also equal to Cluster_hit_time_cm

          if(cluster[i_c].View()==2){
            plane=2;
            pointdistance=sqrt((pow(cluster_hit_z-Z_reco,2))+(pow(cluster_hit_x-X_reco,2)));

          }
          if (cluster[i_c].View()==1){
            plane=1;

            pointdistance=sqrt((pow(cluster_hit_z-V_wire_cm,2))+(pow(cluster_hit_x-V_time_cm,2)));
          }

          if (cluster[i_c].View()==0){
            plane=0;
            pointdistance=sqrt((pow(cluster_hit_z-U_wire_cm,2))+(pow(cluster_hit_x-U_time_cm,2)));
          }

          if(pointdistance<pointdistance_smallest){//comparison IF loop

            pointdistance_smallest=pointdistance;

            X_reco_smallest=X_reco;
            Z_reco_smallest=Z_reco;

          }

        }//END CLUSTER HIT LOOP



      }//END RECO POINT LOOP

      // std::cout<<"distance_nu_smallest: "<<distance_nu_smallest<<std::endl;
      if(pointdistance_smallest<distance_smallest){
        distance_smallest=pointdistance_smallest;
        X_reco_best=X_reco_smallest; //variables for the coordinates of the nearest reco track
        Z_reco_best=Z_reco_smallest; //variables for the coordinates of the nearest reco track
      }


      //}//END RECO TRACK FOR LOOP







      // if(cluster[i_c].View()==2 && (distance_smallest > f2DcutY || (distance_smallest_nu > 2.0 && track_point_length_smallest>4.0 ) || (distance_smallest_nu > (2.0*track_point_length_smallest/4.0) && track_point_length_smallest<4.0 )  )) {// IF LOOP TO CHECK WHAT PLANE A CLUSTER BELONGS TO
      if(cluster[i_c].View()==2 && (distance_smallest > f2DcutY)){


        Start_Cluster2.push_back((cluster[i_c].StartTick ())-3.0);//added +- 3.0 time tick tolerances
        End_Cluster2.push_back((cluster[i_c].EndTick ())+3.0);
        // Start_Cluster2.push_back((cluster[i_c].StartTick ()));//added +- 3.0 time tick tolerances
        // End_Cluster2.push_back((cluster[i_c].EndTick ()));
        Y_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
        Y_index_vector.push_back(i_c); //Y Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
      }

      if(cluster[i_c].View()==1 && (distance_smallest > f2DcutUV )){//IF LOOP TO CHECK WHAT PLANE A CLUSTER BELONGS TO


        Start_Cluster1.push_back((cluster[i_c].StartTick ())-3.0);
        End_Cluster1.push_back((cluster[i_c].EndTick ())+3.0);
        V_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
        // Start_Cluster1.push_back((cluster[i_c].StartTick ()));
        // End_Cluster1.push_back((cluster[i_c].EndTick ()));
        V_index_vector.push_back(i_c);//V Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
      }

      if(cluster[i_c].View()==0 && (distance_smallest > f2DcutUV )){//IF LOOP TO CHECK WHAT PLANE A CLUSTER BELONGS TO

        //
        Start_Cluster0.push_back((cluster[i_c].StartTick ())-3.0);
        End_Cluster0.push_back((cluster[i_c].EndTick ())+3.0);
        U_clus_hitsize.push_back(clus_hit_assn_v.at(i_c).size());
        // Start_Cluster0.push_back((cluster[i_c].StartTick ()));
        // End_Cluster0.push_back((cluster[i_c].EndTick ()));
        U_index_vector.push_back(i_c);//U Index vector to store the event index for a given cluster. Very important variable for getting cluster-hit associaton
      }


    }//END RECO TRACK FOR LOOP (Thanks to Ivan Lepetic for pointing this out)

    Clustertree->Fill();
  }//end cluster FOR loop for calculating 2-D distance



  V_t_min_common = -9999.; //  start point of the overlap
  V_t_max_common = -9999.; // end point of the overlap

  U_t_min_common = -9999.; //  start point of the overlap
  U_t_max_common = -9999.; // end point of the overlap

  auto const* geom = ::lar::providerFrom<geo::Geometry>();

  for (size_t i = 0; i < Y_index_vector.size(); i++) {//START PLANE MATCHING FOR LOOP (Y-COLLECTION PLANE LOOP). A Match is to be found for this cluster always if it exists

    //Match variables (vectors) cleared or set to crazy values at the beginning of each Y Cluster Loop
    IoU_U.clear();
    IoU_V.clear();
    y_U.clear();
    y_V.clear();

    Y_V_index.clear();
    Y_U_index.clear();

    Int_t Y_Match_Ev_Index=-45;//Resetting the indexes to random -ve integers
    Int_t U_Match_Ev_Index=-46;
    Int_t V_Match_Ev_Index=-47;

    Y_clus_hitSize=Y_clus_hitsize[i];

    distance_smallest3d=1e10;

    deltaY=-9999.0;

    auto hiti = clus_hit_assn_v.at(Y_index_vector[i]);

    Y_Wire=hiti[0]->WireID().Wire;



    Y_cluster_charge=0;
    U_cluster_charge=0;
    V_cluster_charge=0;






    Y_clus_lifetime=End_Cluster2[i]-Start_Cluster2[i];




    V_iou=-9999.0;



    for (size_t j = 0; j < V_index_vector.size(); j++) {//START V-INDUCTION PLANE FOR LOOP)


      // V_clus_hitSize=V_clus_hitsize[j];

      V_t_min_abs = 9600.0;//variable for smallest time of the two clusters
      V_t_max_abs = 0.0;//variable for largest time of the two clusters
      if ( Start_Cluster2[i] < V_t_min_abs )
      V_t_min_abs = Start_Cluster2[i];



      if ( End_Cluster2[i]   > V_t_max_abs )
      V_t_max_abs = End_Cluster2[i];




      if ( Start_Cluster1[j] < V_t_min_abs )
      V_t_min_abs = Start_Cluster1[j];



      if ( End_Cluster1[j]   > V_t_max_abs )
      V_t_max_abs = End_Cluster1[j];





      if (Start_Cluster2[i]<Start_Cluster1[j])
      V_t_min_common = Start_Cluster1[j];
      if (Start_Cluster2[i]>Start_Cluster1[j])
      V_t_min_common = Start_Cluster2[i];

      if (End_Cluster2[i]<End_Cluster1[j])
      V_t_max_common = End_Cluster2[i];
      if (End_Cluster2[i]>End_Cluster1[j])
      V_t_max_common = End_Cluster1[j];





      auto hitj = clus_hit_assn_v.at(V_index_vector[j]);





      geo::WireIDIntersection clusterintersectionV;
      auto Vintersect = geom->WireIDsIntersect(hiti[0]->WireID(),hitj[0]->WireID(),clusterintersectionV);

      if (Vintersect==0)
      continue;


      if ((V_t_min_common > V_t_max_common)){//NO OVERLAP
        V_iou = -1;
      }
      if ((V_t_min_common < V_t_max_common)){//OVERLAP
        V_iou = (V_t_max_common - V_t_min_common) / (V_t_max_abs - V_t_min_abs);


      }

      if (V_iou>0){

        Y_V_index.push_back(j);
        IoU_V.push_back(V_iou);
        y_V.push_back(clusterintersectionV.y);
        // std::cout<<"Y Cluster hit wire: "<<hiti[0]->WireID()<<std::endl;
        // std::cout<<"V Cluster hit wire: "<<hitj[0]->WireID()<<std::endl;
        // std::cout<<"V Intersect: "<<Vintersect<<std::endl;

      }






    }//END V-INDUCTION PLANE FOR LOOP




    U_iou=-9999.0;


    for (size_t k = 0; k < U_index_vector.size(); k++) {//START U-INDUCTION PLANE FOR LOOP

      // U_clus_hitSize=U_clus_hitsize[k];

      U_t_min_abs = 9600.0;//variable for smallest time of the two clusters
      U_t_max_abs = 0.0;//variable for largest time of the two clusters


      if ( Start_Cluster2[i] < U_t_min_abs )
      U_t_min_abs = Start_Cluster2[i];
      // start_tick_cluster20=Start_Cluster2[i];

      if ( End_Cluster2[i]   > U_t_max_abs )
      U_t_max_abs = End_Cluster2[i];
      // end_tick_cluster20=End_Cluster2[i];


      if ( Start_Cluster0[k] < U_t_min_abs )
      U_t_min_abs = Start_Cluster0[k];
      // start_tick_cluster0=Start_Cluster0[k];


      if ( End_Cluster0[k]   > U_t_max_abs )
      U_t_max_abs = End_Cluster0[k];
      // end_tick_cluster0=End_Cluster0[k];


      if (Start_Cluster2[i]<Start_Cluster0[k])
      U_t_min_common = Start_Cluster0[k];
      if (Start_Cluster2[i]>Start_Cluster0[k])
      U_t_min_common = Start_Cluster2[i];

      if (End_Cluster2[i]<End_Cluster0[k])
      U_t_max_common = End_Cluster2[i];
      if (End_Cluster2[i]>End_Cluster0[k])
      U_t_max_common = End_Cluster0[k];


      auto hitk = clus_hit_assn_v.at(U_index_vector[k]);



      geo::WireIDIntersection clusterintersectionU;

      auto Uintersect = geom->WireIDsIntersect(hiti[0]->WireID(),hitk[0]->WireID(),clusterintersectionU);//

      if (Uintersect==0)
      continue;



      if (U_t_min_common > U_t_max_common)//NO OVERLAP
      U_iou = -1;
      if ((U_t_min_common < U_t_max_common) )//OVERLAP
      {
        U_iou = (U_t_max_common - U_t_min_common) / (U_t_max_abs - U_t_min_abs);


      }

      if (U_iou>0){
        Y_U_index.push_back(k);
        IoU_U.push_back(U_iou);
        y_U.push_back(clusterintersectionU.y);
        // std::cout<<"Y Cluster hit wire: "<<hiti[0]->WireID()<<std::endl;
        // std::cout<<"U Cluster hit wire: "<<hitk[0]->WireID()<<std::endl;
        // std::cout<<"U Intersect: "<<Uintersect<<std::endl;
      }

    }//END U-INDUCTION PLANE FOR LOOP


    U_match_multiplicity=IoU_U.size();
    V_match_multiplicity=IoU_V.size();


    cluster_y_new=-9999.0;//variable for storing the value of the y coordinate initialized to a crazy value

    //********************V=0,U=0 bin Start********************//

    if (IoU_U.size()==0 && IoU_V.size()==0 ) { //V=0,U=0 bin
      cluster_y_new=-9999.0;
      V_biggest_iou=-1.0;
      U_biggest_iou=-1.0;
      YV_iou=-1.0;
      YU_iou=-1.0;
    }//V=0,U=0 bin

    //********************V=0,U=0 bin End********************//




    //********************V=0,U=1 bin Start********************//

    if (IoU_U.size()==1 && IoU_V.size()==0 ) { //V=0,U=1 bin


      cluster_y_new=y_U.at(0);
      YU_iou=IoU_U.at(0);
      U_biggest_iou=IoU_U.at(0);
      YV_iou=-1.0;
      V_biggest_iou=-1.0;

      U_Match_Ev_Index=U_index_vector[Y_U_index[0]];//index of the cluster from the main event
      auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[0]]);

      for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP
        //  charge=0;
        U_charge = hitnewU->Integral();
        U_cluster_charge += U_charge;
        // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;
      }//END CLUSTER HIT LOOP


    }//V=0,U=1 bin


    //********************V=0,U=1 bin End********************//





    //********************V=1,U=0 bin Start********************//


    if (IoU_U.size()==0 && IoU_V.size()==1 ) { //V=1,U=0 bin


      cluster_y_new=y_V.at(0);
      YV_iou=IoU_V.at(0);
      V_biggest_iou=IoU_V.at(0);
      YU_iou=-1.0;
      U_biggest_iou=-1.0;

      V_Match_Ev_Index=V_index_vector[Y_V_index[0]];//index of the cluster from the main event
      auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[0]]);

      for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

        V_charge = hitnewV->Integral();
        V_cluster_charge += V_charge;
        // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;


      }//END CLUSTER HIT LOOP


    }//V=1,U=0 bin


    //********************V=0,U=1 bin End********************//





    //********************V=1,U=1 bin Start********************//


    if (IoU_U.size()==1 && IoU_V.size()==1 ) { //V=1,U=1 bin

      YU_iou=IoU_U.at(0);
      YV_iou=IoU_V.at(0);

      deltaY=y_V.at(0)-y_U.at(0);
      // std::cout<<"Y-Coordinate from V-Plane"<<y_V.at(0)<<std::endl;
      // std::cout<<"Y-Coordinate from U-Plane"<<y_U.at(0)<<std::endl;
      // std::cout<<"Delta Y: "<<deltaY<<std::endl;
      deltaY_smallest=deltaY;

      if (abs(deltaY)< fdeltaY){
        auto avg=0.5*(y_V.at(0)+y_U.at(0));
        cluster_y_new=avg;



        U_biggest_iou=IoU_U.at(0);
        V_biggest_iou=IoU_V.at(0);

        V_Match_Ev_Index=V_index_vector[Y_V_index[0]];
        auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[0]]);

        for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

          V_charge = hitnewV->Integral();
          V_cluster_charge += V_charge;
          // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;


        }//END CLUSTER HIT LOOP

        U_Match_Ev_Index=U_index_vector[Y_U_index[0]];
        auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[0]]);

        for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP

          U_charge = hitnewU->Integral();
          U_cluster_charge += U_charge;
          // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;


        }//END CLUSTER HIT LOOP

      }
      else{
        cluster_y_new=-9999.0;
      }

      start_tick_cluster2=Start_Cluster2[i];
      end_tick_cluster2=End_Cluster2[i];
      start_tick_cluster1=Start_Cluster1[Y_V_index[0]];
      end_tick_cluster1=End_Cluster1[Y_V_index[0]];
      start_tick_cluster0=Start_Cluster0[Y_U_index[0]];
      end_tick_cluster0=End_Cluster0[Y_U_index[0]];
      V_clus_hitSize=V_clus_hitsize[Y_V_index[0]];
      U_clus_hitSize=U_clus_hitsize[Y_U_index[0]];




    }//V=1,U=1 bin



    //********************V=1,U=1 bin End********************//



    //********************V>1,U=0 bin Start********************//


    if (IoU_U.size()==0 && IoU_V.size()>1 ) { //V>1,U=0 bin

      U_biggest_iou=-1.0;
      YU_iou=-1.0;
      Double_t V_biggest_ious=-1.0;


      for (size_t Vious=0; Vious<IoU_V.size(); Vious++) {

        YV_iou=IoU_V.at(Vious);

        if ((YV_iou > V_biggest_ious)){
          V_biggest_ious=YV_iou;
          cluster_y_new=y_V.at(Vious);
          V_biggest_iou=V_biggest_ious;
          V_Match_Ev_Index=V_index_vector[Y_V_index[Vious]];
          auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[Vious]]);

          for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

            V_charge = hitnewV->Integral();
            V_cluster_charge += V_charge;
            // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;


          }//END CLUSTER HIT LOOP



        }


      }

    }//V>1,U=0 bin


    //********************V>1,U=0 bin End********************//



    //********************V=0,U=>1 bin Start********************//

    if (IoU_V.size()==0 && IoU_U.size()>1 ) { //V=0,U>1 bin

      YV_iou=-1.0;
      V_biggest_iou=-1.0;
      Double_t U_biggest_ious=-1.0;


      for (size_t Uious=0; Uious<IoU_U.size(); Uious++) {

        YU_iou=IoU_U.at(Uious);

        if ((YU_iou > U_biggest_ious)){
          U_biggest_ious=YU_iou;
          cluster_y_new=y_U.at(Uious);
          U_biggest_iou=U_biggest_ious;

          U_Match_Ev_Index=U_index_vector[Y_U_index[Uious]];


          auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[Uious]]);

          for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP

            U_charge = hitnewU->Integral();
            U_cluster_charge += U_charge;
            // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;


          }//END CLUSTER HIT LOOP

        }
      }

    }//V=0,U>1 bin


    //********************V=0,U=>1 bin End********************//



    //********************V=1,U>1 bin Start********************//

    if (IoU_V.size()==1 && IoU_U.size()>1 ) { //V=1,U>1 bin



      YV_iou=IoU_V.at(0);
      V_biggest_iou=IoU_V.at(0);
      Double_t U_biggest_ious=-1.0;
      deltaY_smallest=1e5;




      for (size_t Uious=0; Uious<IoU_U.size(); Uious++) {

        YU_iou=IoU_U.at(Uious);


        if ((YU_iou > U_biggest_ious)){
          U_biggest_ious=YU_iou;

          U_biggest_iou=U_biggest_ious;
        }

        for (size_t Vious=0; Vious<IoU_V.size(); Vious++){


          deltaY=((y_V.at(Vious))-(y_U.at(Uious)));

          if(abs(deltaY)<deltaY_smallest){
            deltaY_smallest=abs(deltaY);

            if(abs(deltaY_smallest)<fdeltaY){
              cluster_y_new=0.5*(y_V.at(Vious)+y_U.at(Uious));

              U_Match_Ev_Index=U_index_vector[Y_U_index[Uious]];

              auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[Uious]]);

              for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP

                U_charge = hitnewU->Integral();
                U_cluster_charge += U_charge;
                // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;


              }//END CLUSTER HIT LOOP

              auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[0]]);

              V_Match_Ev_Index=V_index_vector[Y_V_index[0]];

              for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

                V_charge = hitnewV->Integral();
                V_cluster_charge += V_charge;
                // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;


              }//END CLUSTER HIT LOOP

            }

            else{
              cluster_y_new=-9999.0;
            }

          }

        }

      }

    }//V=1,U>1 bin


    //********************V=1,U>1 bin End********************//



    //********************V>1,U=1 bin Start********************//
    if (IoU_U.size()==1 && IoU_V.size()>1 ) { //V>1,U=1 bin

      deltaY_smallest=10e5;
      YU_iou=IoU_U.at(0);
      U_biggest_iou=IoU_U.at(0);
      Double_t V_biggest_ious=-1.0;




      for (size_t Vious=0; Vious<IoU_V.size(); Vious++) {

        YV_iou=IoU_V.at(Vious);

        if ((YV_iou > V_biggest_ious)){
          V_biggest_ious=YV_iou;

          V_biggest_iou=V_biggest_ious;
        }

        for (size_t Uious=0; Uious<IoU_U.size(); Uious++) {


          deltaY=((y_V.at(Vious))-(y_U.at(Uious)));



          if(abs(deltaY)<deltaY_smallest){
            deltaY_smallest=abs(deltaY);

            if(abs(deltaY_smallest)<fdeltaY){
              cluster_y_new=0.5*(y_V.at(Vious)+y_U.at(Uious));


              V_Match_Ev_Index=V_index_vector[Y_V_index[Vious]];
              auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[Vious]]);

              for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

                V_charge = hitnewV->Integral();
                V_cluster_charge += V_charge;
                // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;

              }//END CLUSTER HIT LOOP



              U_Match_Ev_Index=U_index_vector[Y_U_index[0]];
              auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[0]]);

              for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP

                U_charge = hitnewU->Integral();
                U_cluster_charge += U_charge;
                // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;


              }//END CLUSTER HIT LOOP

            }
          }


        }

      }

    }//V>1,U=1 bin

    //********************V>1,U=1 bin End********************//





    //********************V>1,U>1 bin Start********************//


    if (IoU_V.size()>1 && IoU_U.size()>1 ) { //V>1,U>1 bin

      deltaY_smallest=10e5;
      for (size_t Vious=0; Vious<IoU_V.size(); Vious++) {

        Double_t U_biggest_ious=-1.0;
        Double_t V_biggest_ious=-1.0;

        YV_iou=IoU_V.at(Vious);
        if ((YV_iou > V_biggest_ious)){
          V_biggest_ious=YV_iou;

          V_biggest_iou=V_biggest_ious;

        }


        for (size_t Uious=0; Uious<IoU_U.size(); Uious++) {


          YU_iou=IoU_U.at(Uious);

          if ((YU_iou > U_biggest_ious)){
            U_biggest_ious=YU_iou;

            U_biggest_iou=U_biggest_ious;

          }


          deltaY=((y_V.at(Vious))-(y_U.at(Uious)));

          if(abs(deltaY)<deltaY_smallest){
            deltaY_smallest=abs(deltaY);


            if (deltaY_smallest<fdeltaY){
              cluster_y_new=0.5*(y_V.at(Vious)+y_U.at(Uious));

              U_Match_Ev_Index=U_index_vector[Y_U_index[Uious]];
              auto hitknew = clus_hit_assn_v.at(U_index_vector[Y_U_index[Uious]]);

              for (auto const& hitnewU : hitknew) {//START CLUSTER HIT LOOP

                U_charge = hitnewU->Integral();
                U_cluster_charge += U_charge;
                // U_cluster_energy += U_charge*240*23.6*1e-6/0.5;


              }//END CLUSTER HIT LOOP


              V_Match_Ev_Index=V_index_vector[Y_V_index[Vious]];
              auto hitjnew = clus_hit_assn_v.at(V_index_vector[Y_V_index[Vious]]);

              for (auto const& hitnewV : hitjnew) {//START CLUSTER HIT LOOP

                V_charge = hitnewV->Integral();
                V_cluster_charge += V_charge;
                // V_cluster_energy += V_charge*240*23.6*1e-6/0.5;

              }//END CLUSTER HIT LOOP

            }

            else{
              cluster_y_new=-9999.0;

            }

          }

        }

      }

    }//V>1,U>1 bin






    for (auto const& hitnewY : hiti) {//START Y CLUSTER HIT LOOP

      Y_charge = hitnewY->Integral();
      Y_cluster_charge += Y_charge;
      // Y_cluster_energy += Y_charge*240*23.6*1e-6/0.5;

      Y_cluster_3d_hit_z = hitnewY->WireID().Wire * wire2cm;//Also equal to Cluster_hit_wire_cm
      Y_cluster_3d_hit_x = (hitnewY->PeakTime() * time2cm)-44.575 ;//Also equal to Cluster_hit_time_cm

      Y_cluster_3d_hit_y=cluster_y_new;

    }//END Y CLUSTER HIT LOOP




    Matchingtree->Fill();


    start_tick_cluster2=0;
    start_tick_cluster1=0;
    start_tick_cluster0=0;
    end_tick_cluster2=0;
    end_tick_cluster1=0;
    end_tick_cluster0=0;


    Y_Match_Ev_Index=Y_index_vector[i];





    if (cluster_y_new ==-9999.0)
    continue;



    Double_t weighted_numX = 0;
    Double_t weighted_numZ = 0;
    Double_t Y_cluster_charge_w = 0;
    Double_t Y_charge_w = 0;





    for (auto const& hitnewY : hiti) {//START Y CLUSTER HIT LOOP

      Y_charge_w = hitnewY->Integral();
      Y_cluster_charge_w += Y_charge_w;

      weighted_numX +=Y_charge_w*((hitnewY->PeakTime() * time2cm)-44.575) ;
      weighted_numZ +=Y_charge_w*(hitnewY->WireID().Wire * wire2cm);


    }//END Y CLUSTER HIT LOOP






    double x,y,z,errXYZ,chisq;
    Double32_t xyz[3];
    Double32_t ErrXYZ[6];
    Double32_t  Chisq;

    x= weighted_numX/Y_cluster_charge_w;
    y= cluster_y_new;
    z= weighted_numZ/Y_cluster_charge_w;
    errXYZ=0;
    chisq=0;



    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    ErrXYZ[0]=errXYZ;
    ErrXYZ[1]=errXYZ;
    ErrXYZ[2]=errXYZ;
    ErrXYZ[3]=errXYZ;
    ErrXYZ[4]=errXYZ;
    ErrXYZ[5]=errXYZ;
    Chisq=chisq;



    recob::SpacePoint newpt(xyz,ErrXYZ,Chisq);
    SpacePoint_v->emplace_back(newpt);





    if((U_match_multiplicity==1 && V_match_multiplicity==0) || (U_match_multiplicity>1 && V_match_multiplicity==0)){
      art::Ptr<recob::Cluster> ClusPtrY(cluster_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Cluster> ClusPtrU(cluster_handle,U_Match_Ev_Index);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrY, *sps_clus_assn_v);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrU, *sps_clus_assn_v);

      art::Ptr<recob::Hit> HitPtrY(hit_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Hit> HitPtrU(hit_handle,U_Match_Ev_Index);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrY, *sps_hit_assn_v);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrU, *sps_hit_assn_v);

    }

    if((U_match_multiplicity==0 && V_match_multiplicity==1) || (U_match_multiplicity==0 && V_match_multiplicity>1)){
      art::Ptr<recob::Cluster> ClusPtrY(cluster_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Cluster> ClusPtrV(cluster_handle,V_Match_Ev_Index);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrY, *sps_clus_assn_v);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrV, *sps_clus_assn_v);

      art::Ptr<recob::Hit> HitPtrY(hit_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Hit> HitPtrV(hit_handle,V_Match_Ev_Index);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrY, *sps_hit_assn_v);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrV, *sps_hit_assn_v);
    }


    if((U_match_multiplicity==1 && V_match_multiplicity>1) || (U_match_multiplicity>1 && V_match_multiplicity==1) || (U_match_multiplicity>1 && V_match_multiplicity>1) || (U_match_multiplicity==1 && V_match_multiplicity==1) ){
      art::Ptr<recob::Cluster> ClusPtrY(cluster_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Cluster> ClusPtrU(cluster_handle,U_Match_Ev_Index);
      art::Ptr<recob::Cluster> ClusPtrV(cluster_handle,V_Match_Ev_Index);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrY, *sps_clus_assn_v);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrV, *sps_clus_assn_v);
      util::CreateAssn(*this, e, *SpacePoint_v,ClusPtrU, *sps_clus_assn_v);

      art::Ptr<recob::Hit> HitPtrY(hit_handle,Y_Match_Ev_Index);
      art::Ptr<recob::Hit> HitPtrU(hit_handle,U_Match_Ev_Index);
      art::Ptr<recob::Hit> HitPtrV(hit_handle,V_Match_Ev_Index);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrY, *sps_hit_assn_v);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrV, *sps_hit_assn_v);
      // util::CreateAssn(*this, e, *SpacePoint_v,HitPtrU, *sps_hit_assn_v);


    }


  }//END PLANE MATCHING FOR LOOP (Y-COLLECTION PLANE LOOP)

  //Lets clear these variables at the end of every event

  Start_Cluster0.clear();
  End_Cluster0.clear();
  Start_Cluster1.clear();
  End_Cluster1.clear();
  Start_Cluster2.clear();
  End_Cluster2.clear();
  Y_clus_hitsize.clear();
  V_clus_hitsize.clear();
  U_clus_hitsize.clear();
  Y_index_vector.clear();
  V_index_vector.clear();
  U_index_vector.clear();




  e.put(std::move(SpacePoint_v));
  e.put(std::move(sps_clus_assn_v));
  // e.put(std::move(sps_hit_assn_v));
}//END EVENT LOOP

void BlipReco3D::beginJob()
{
  // Implementation of optional member function here.

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();


  wire2cm = geom->WirePitch(0,0,0);
  time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );



  art::ServiceHandle<art::TFileService> tfs;

  Clustertree = tfs->make<TTree>("Clustertree",    "Clustertree");
  Matchingtree = tfs->make<TTree>("Matchingtree",    "Matchingtree");

  Clustertree->Branch("cluster_hit_z",&cluster_hit_z,"cluster_hit_z/D");
  Clustertree->Branch("cluster_hit_x",&cluster_hit_x,"cluster_hit_x/D");
  Clustertree->Branch("Z_reco_best",&Z_reco_best,"Z_reco_best/D");
  Clustertree->Branch("X_reco_best",&X_reco_best,"X_reco_best/D");
  Clustertree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
  Clustertree->Branch("plane",&plane,"plane/I");


  Matchingtree->Branch("V_biggest_iou",&V_biggest_iou,"V_biggest_iou/D");
  Matchingtree->Branch("V_match_multiplicity",&V_match_multiplicity,"V_match_multiplicity/I");
  Matchingtree->Branch("U_biggest_iou",&U_biggest_iou,"U_biggest_iou/D");
  Matchingtree->Branch("U_match_multiplicity",&U_match_multiplicity,"U_match_multiplicity/I");
  Matchingtree->Branch("Y_clus_lifetime",&Y_clus_lifetime,"Y_clus_lifetime/D");
  Matchingtree->Branch("Y_clus_hitSize",&Y_clus_hitSize,"Y_clus_hitSize/D");
  Matchingtree->Branch("V_clus_hitSize",&V_clus_hitSize,"V_clus_hitSize/D");
  Matchingtree->Branch("U_clus_hitSize",&U_clus_hitSize,"U_clus_hitSize/D");
  Matchingtree->Branch("Y_cluster_3d_hit_z",&Y_cluster_3d_hit_z,"Y_cluster_3d_hit_z/D");
  Matchingtree->Branch("Y_cluster_3d_hit_x",&Y_cluster_3d_hit_x,"Y_cluster_3d_hit_x/D");
  Matchingtree->Branch("Y_cluster_3d_hit_y",&Y_cluster_3d_hit_y,"Y_cluster_3d_hit_y/D");
  Matchingtree->Branch("Z_reco_best3d",&Z_reco_best3d,"Z_reco_best3d/D");
  Matchingtree->Branch("X_reco_best3d",&X_reco_best3d,"X_reco_best3d/D");
  Matchingtree->Branch("Y_reco_best3d",&Y_reco_best3d,"Y_reco_best3d/D");
  Matchingtree->Branch("distance_smallest3d",&distance_smallest3d,"distance_smallest3d/D");
  Matchingtree->Branch("Y_cluster_charge",&Y_cluster_charge,"Y_cluster_charge/D");
  Matchingtree->Branch("V_cluster_charge",&V_cluster_charge,"V_cluster_charge/D");
  Matchingtree->Branch("U_cluster_charge",&U_cluster_charge,"U_cluster_charge/D");
  Matchingtree->Branch("YV_iou",&YV_iou,"YV_iou/D");
  Matchingtree->Branch("YU_iou",&YU_iou,"YU_iou/D");
  Matchingtree->Branch("deltaY",&deltaY,"deltaY/D");
  Matchingtree->Branch("deltaY_smallest",&deltaY_smallest,"deltaY_smallest/D");
  Matchingtree->Branch("Y_Wire",&Y_Wire,"Y_Wire/i");
  Matchingtree->Branch("start_tick_cluster0",&start_tick_cluster0,"start_tick_cluster0/D");
  Matchingtree->Branch("start_tick_cluster1",&start_tick_cluster1,"start_tick_cluster1/D");
  Matchingtree->Branch("start_tick_cluster2",&start_tick_cluster2,"start_tick_cluster2/D");
  Matchingtree->Branch("end_tick_cluster0",&end_tick_cluster0,"end_tick_cluster0/D");
  Matchingtree->Branch("end_tick_cluster1",&end_tick_cluster1,"end_tick_cluster1/D");
  Matchingtree->Branch("end_tick_cluster2",&end_tick_cluster2,"end_tick_cluster2/D");

}

void BlipReco3D::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(BlipReco3D)
