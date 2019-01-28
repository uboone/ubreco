////////////////////////////////////////////////////////////////////////
// Class:       ClusterTrackDistance
// Plugin Type: analyzer (art v2_05_01)
// File:        ClusterTrackDistance_module.cc
//
// Generated at Mon Jul 16 11:12:21 2018 by Avinay Bhat using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
//#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <memory>
#include <map>
// Services
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
//#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"

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


class ClusterTrackDistance;


class ClusterTrackDistance : public art::EDAnalyzer {
public:
    explicit ClusterTrackDistance(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ClusterTrackDistance(ClusterTrackDistance const &) = delete;
    ClusterTrackDistance(ClusterTrackDistance &&) = delete;
    ClusterTrackDistance & operator = (ClusterTrackDistance const &) = delete;
    ClusterTrackDistance & operator = (ClusterTrackDistance &&) = delete;
    
    // Required functions.
    void analyze(art::Event const & e) override;
    
    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    
private:
    
    // Declare member data here.
    
    
    std::string fMC_track_tag,fReco_track_tag,fCluster_tag,fHit_tag;
    
    
    
    Int_t ClusterStartWire=0;
    Int_t ClusterEndWire=0;
    Double_t ClusterStartTick=0;
    Double_t ClusterEndTick=0;
    
    Double_t X_reco=0.0;
    Double_t Y_reco=0.0;
    Double_t Z_reco=0.0;
    
    Double_t cluster_hit_z=0.0;
    Double_t cluster_hit_x=0.0;
    
    Double_t X_reco_best=0.0;
    Double_t Y_reco_best=0.0;
    Double_t Z_reco_best=0.0;
    
    
    
    //   Double_t Wire_reco=0;
    //    Double_t Tick_reco=0.0;
    
    Double_t pointdistance=0;
    Double_t       X_reco_smallest=0;
    Double_t       Z_reco_smallest=0;
    
    
    Double_t pointdistance_smallest;
    
    //    Double_t distance;
    Double_t distance_smallest;
    Double_t charge;
    //   Double_t charge1;
    Double_t cluster_charge;
    Double_t cluster_energy;
    
    Double_t wire2cm,time2cm;
    
    TTree *Clustertree;
    
    
    
    
    
    
};





ClusterTrackDistance::ClusterTrackDistance(fhicl::ParameterSet const & p)
: EDAnalyzer(p)  // ,
// More initializers here.
{
    
    //All tags get filled in the fcl file
    fMC_track_tag   = p.get<std::string>("mctrack_tag"    );
    fReco_track_tag = p.get<std::string>("recotrack_tag"  );
    fCluster_tag = p.get<std::string>("cluster_tag"  );
    fHit_tag = p.get<std::string>("hit_tag"  );
    
    
}

void ClusterTrackDistance::analyze(art::Event const & e)
{
    
    //    std::cout << "on event " << e.event() << " at run " << e.run() << std::endl;
    
    // Implementation of required member function here.
    art::Handle<std::vector<recob::Cluster> > cluster_handle;
    e.getByLabel(fCluster_tag,cluster_handle);
    
    art::Handle<std::vector<recob::Track> > recotrack_handle;
    e.getByLabel(fReco_track_tag,recotrack_handle);
    
    art::Handle<std::vector<recob::Hit> > hit_handle;
    e.getByLabel(fHit_tag,hit_handle);
    
    art::FindMany<recob::Hit> clus_hit_assn_v(cluster_handle, e, fCluster_tag);
    
    Int_t clustercounter=0; //Variable to keep track of cluster index
    
    
    for (size_t i_c = 0, size_cluster = cluster_handle->size(); i_c != size_cluster; ++i_c) { //START CLUSTER FOR LOOP
        clustercounter++;
        //    cout<<"***************************************************Cluster Counter: "<<clustercounter<<endl;
        
        
        distance_smallest=1e10; //Variable for distance between a cluster and nearest reco track, initialized to a large number for comparison
        cluster_charge=0.0;//Cluster variable initialized to zero
        cluster_energy=0.0;//Cluster variable initialized to zero
        
        
        auto hits = clus_hit_assn_v.at(i_c);
        
        Int_t trackcounter=0;//Variable to keep track of cluster index
        //Int_t trackcounter_smallest=0;
        
        
        for (size_t i_t = 0, size_track = recotrack_handle->size(); i_t != size_track; ++i_t) {//START RECO TRACK FOR LOOP
            trackcounter++;
            //       cout<<"***************************************Track Counter: "<<trackcounter<<endl;
            
            
            auto const& track = recotrack_handle->at(i_t);
            
            pointdistance_smallest=1e10;////Variable for distance between a cluster hit and a reco track point, initialized to a large number for comparison
            for(size_t m=0;m<(track.NumberTrajectoryPoints());m++){//START RECO POINT LOOP
                
                
                X_reco=track.LocationAtPoint(m).X();
                Y_reco=track.LocationAtPoint(m).Y();
                Z_reco=track.LocationAtPoint(m).Z();
                //   cout<<"X_reco: "<<X_reco<<endl;
                
                
                for (auto const& hit : hits) {//START CLUSTER HIT LOOP
                    
                    cluster_hit_z = hit->WireID().Wire * wire2cm;//What is this?
                    cluster_hit_x = (hit->PeakTime() * time2cm)-44.575 ;//What is this?
                    pointdistance=sqrt((pow(cluster_hit_z-Z_reco,2))+(pow(cluster_hit_x-X_reco,2)));
                    
                    //    cout<<"pointdistance: "<<pointdistance<<endl;
                    
                    if(pointdistance<pointdistance_smallest){//comparison IF loop
                        
                        pointdistance_smallest=pointdistance;
                        
                        X_reco_smallest=X_reco;
                        Z_reco_smallest=Z_reco;
                        
                    }
                    
                }//END CLUSTER HIT LOOP
                
                
                //       cout<<"Pointdistance_smallest1: "<<pointdistance_smallest<<endl;
                
            }//END RECO POINT LOOP
            
            
            if(pointdistance_smallest<distance_smallest){
                distance_smallest=pointdistance_smallest;
                //   trackcounter_smallest=trackcounter;
                X_reco_best=X_reco_smallest;
                Z_reco_best=Z_reco_smallest;
            }
            //      cout<<"distance_smallest: "<<distance_smallest<<endl;
            
        }//END RECO TRACK FOR LOOP
        
        
        
        
        for (auto const& hit : hits) {//START CLUSTER HIT LOOP
            //  charge=0;
            charge = hit->Integral();
            cluster_charge += charge;
            cluster_energy += charge*240*23.6*1e-6/0.5;
            
        }//END CLUSTER HIT LOOP
        
        
        
        //      cout<<"Smallest Distance: "<<distance_smallest<<endl;
        //      cout<<"X_reco_best: "<<X_reco_best<<endl;
        //     cout<<"cluster_energy: "<<cluster_energy<<" MeV"<<endl;
        
        Clustertree->Fill();
    }//END CLUSTER FOR LOOP
    
    
}


void ClusterTrackDistance::beginJob()
{ // get detector specific properties
    
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    wire2cm = geom->WirePitch(0,0,0);
    time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    
    
    
    art::ServiceHandle<art::TFileService> tfs;
    
    Clustertree = tfs->make<TTree>("Clustertree",    "Clustertree");
    
    Clustertree->Branch("cluster_hit_z",&cluster_hit_z,"cluster_hit_z/D");
    Clustertree->Branch("cluster_hit_x",&cluster_hit_x,"cluster_hit_x/D");
    Clustertree->Branch("Z_reco_best",&Z_reco_best,"Z_reco_best/D");
    Clustertree->Branch("X_reco_best",&X_reco_best,"X_reco_best/D");
    Clustertree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
    Clustertree->Branch("cluster_charge",&cluster_charge,"cluster_charge/D");
    Clustertree->Branch("cluster_energy",&cluster_energy,"cluster_energy/D");
}

void ClusterTrackDistance::endJob()
{
    // Implementation of optional member function here.
}









DEFINE_ART_MODULE(ClusterTrackDistance)
