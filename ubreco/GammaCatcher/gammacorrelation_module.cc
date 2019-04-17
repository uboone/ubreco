////////////////////////////////////////////////////////////////////////
// Class:       gammacorrelation
// Plugin Type: analyzer (art v3_01_02)
// File:        gammacorrelation_module.cc
//
// Generated at Wed Mar 20 10:42:52 2019 by Avinay Bhat using cetskelgen
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


#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
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
  bool endSubRun(art::SubRun &subrun);

private:

  // Declare member data here.

  std::string fHit_tag,fpfparticle_tag,fvertex_tag,fsps_tag,fcluster_tag;

  Double32_t sps_x,sps_y,sps_z,sps_hit_charge,sps_cluster_charge,sps_cluster_charge10,sps_cluster_charge20,sps_cluster_charge50;
  Double_t distance, distance_smallest,Event_cluster_charge;
  Double_t Vertex_x,Vertex_y,Vertex_z;
  TRandom3 rand;
  Double_t _rand_vtx_x, _rand_vtx_y, _rand_vtx_z, distance_rand_vtx, distance_smallest_rand_vtx;
  Int_t neutrinos,N_sps,N_Event,N_Run,N_SubRun,N_sps10,N_sps20,N_sps50;


  TTree *Event_Correlationtree;
  TTree *Sps_Correlationtree;

  TTree* _subrun_tree;
    int _run_sr;                  // The run number
    int _sub_sr;                  // The subRun number
    float _pot;                   // The total amount of POT for the current sub run


  Double_t fidVolMinX =    0; //Fiducial Volume dimensions for MicroBooNE
  Double_t fidVolMaxX =  256;
  Double_t fidVolMinY = -116;
  Double_t fidVolMaxY =  116;
  Double_t fidVolMinZ =    0;
  Double_t fidVolMaxZ = 1030;

  bool fData;
  art::InputTag fMCTproducer;

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
  fData     = p.get< bool >             ("IsData");
  fMCTproducer = p.get< art::InputTag > ("MCTproducer");

}

void gammacorrelation::analyze(art::Event const& e)
{


  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() << " SubRun Number: " << e.subRun() <<"**********************"<< std::endl;

  N_Event=e.event();
  N_Run=e.run();
  N_SubRun=e.subRun();

  art::Handle<std::vector<recob::Hit> > hit_handle;
  e.getByLabel(fHit_tag,hit_handle);

  art::Handle<std::vector<recob::Cluster> > cluster_handle;
  e.getByLabel(fcluster_tag,cluster_handle);

  art::Handle<std::vector<recob::PFParticle> > pfparticle_handle;
  e.getByLabel(fpfparticle_tag,pfparticle_handle);

  art::Handle<std::vector<recob::Vertex> > vertex_handle;
  e.getByLabel(fvertex_tag,vertex_handle);

  art::FindMany<recob::Vertex> pfp_vertex_assn_v(pfparticle_handle, e, fvertex_tag);

  art::Handle<std::vector<recob::SpacePoint> > spacepoint_handle;
  e.getByLabel(fsps_tag,spacepoint_handle);

  art::FindMany<recob::Cluster> sps_clus_assn_v(spacepoint_handle, e, fsps_tag);

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


      auto ass_vtx_v  =pfp_vertex_assn_v.at( p );
      if (ass_vtx_v.size() != 1)
      std::cout << "ERROR. Neutrino not associated with a single vertex..." << std::endl;
      nuvtx = *(ass_vtx_v.at(0));


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


  for(size_t s=0;s<spacepoint_handle->size();s++){
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
  }


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

  _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
  _subrun_tree->Branch("run"   , &_run_sr   , "run/I");
  _subrun_tree->Branch("subRun", &_sub_sr   , "subRun/I");

  if (!fData)
        _subrun_tree->Branch("pot", &_pot, "pot/F");
}

void gammacorrelation::endJob(){}

bool gammacorrelation::endSubRun(art::SubRun &subrun)
{
    std::cout << "DAVIDC END SUBRUN" << std::endl;
    if (!fData)
    {
        std::cout << "DAVIDC MC!" << std::endl;
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(fMCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
        std::cout << "POT is " << _pot << std::endl;
        std::cout << "[LArPandoraExternalEventBuilding::endSubRun] Storing POT info!" << std::endl;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();
    return true;
}


DEFINE_ART_MODULE(gammacorrelation)
