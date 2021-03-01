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
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
//#include "larevt/SpaceADCServices/SpaceADCService.h"
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

  Int_t pdg_particle=-123456;
  Int_t pdg_particle_pr=-123456;
  Double_t MC_Particle_Energy=0.0;
  Double_t MC_Particle_Energy_pr=0.0;
  Double_t Event_Energy_total=0.0;
  Double_t Event_Energy_total_pr=0.0;
  // Double_t Event_Energy_5=0.0;
  // Double_t Event_Energy_10=0.0;
  // Double_t Event_Energy_20=0.0;
  // Double_t Event_Energy_30=0.0;
  // Double_t Event_Energy_40=0.0;
  // Double_t Event_Energy_50=0.0;
  // Double_t Event_Energy_100=0.0;
  Double_t position_X=-9999.0;
  Double_t position_Y=-9999.0;
  Double_t position_Z=-9999.0;
  Double_t position_X_pr=-9999.0;
  Double_t position_Y_pr=-9999.0;
  Double_t position_Z_pr=-9999.0;
  // Double_t distance_center=-9999.0;
  //Double_t MC_Particle_Energy_MeV=0.0;
  //Double_t Event_Energy_total_MeV=0.0;
  // Double_t Event_Energy_5_MeV=0.0;
  // Double_t Event_Energy_10_MeV=0.0;
  // Double_t Event_Energy_20_MeV=0.0;
  // Double_t Event_Energy_30_MeV=0.0;
  // Double_t Event_Energy_40_MeV=0.0;
  // Double_t Event_Energy_50_MeV=0.0;
  // Double_t Event_Energy_100_MeV=0.0;
  Int_t StatusCode=-9999;
  // Double_t distance_true_reco=999999.0;
  // Double_t distance_true_reco_smallest=999999.0;
  //Double_t MC_Particle_Energy_MeV_smallest=99999.0;
  //Double_t sps_ADC_smallest=99999.0;

  TTree *Eventtree;
  TTree *MCParticletree;
  TTree *Hittree;


  // s = std.string();

  // Declare member data here.

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  std::string fHit_tag,fpfparticle_tag,fvertex_tag,fsps_tag,fcluster_tag,fReco_track_tag,fMatch_tag,endprocess,process;

  Int_t evttime=0;

  //Double32_t sps_x,sps_y,sps_z,sps_hit_ADC,bktrkd_particle_energy,bktrkd_particle_energy_per_sps,sps_ADC_U,sps_ADC_V,sps_ADC_Y,sps_ADC_match,sps_ADC10,sps_ADC20,sps_ADC50,hit_ADC_U;
  Double32_t sps_x,sps_y,sps_z,bktrkd_particle_energy,bktrkd_particle_energy_per_sps,sps_ADC_U,sps_ADC_V,sps_ADC_Y,sps_ADC10,sps_ADC20,sps_ADC50,hit_ADC_U;

  Double_t Event_sps_ADC_U,Event_sps_ADC_V,Event_sps_ADC_Y,hit_ADC_U_sum,numElectrons,numElectrons_sum;

  TRandom3 rand;

  //Int_t neutrinos,N_Event,N_Run,N_SubRun,N_sps, N_total_e,N_total_pr_e, N_total_bktrk_e,N_true_pr_bktrk_e,hit_counter;
  Int_t neutrinos,N_Event,N_Run,N_SubRun,N_sps, N_total_e,N_total_pr_e, hit_counter;

  float _maxTrkLen;
  int   _neutrinoshowers;
  int   _neutrinotracks;
  float _muon_px, _muon_py, _muon_pz;
  //Double_t tracklength = 0.0,X_reco=0.0,Y_reco=0.0,Z_reco=0.0;

  Double_t matched_sps_true_energy, matched_sps_reco_ADC;

  //TTree *Event_Correlationtree;
  TTree *Sps_Correlationtree;

  bool isData,fData;

  /** Setup root trees  */
  TTree *potTree;
  double sr_pot = 0;
  int sr_run = 0;
  int sr_sub_run = 0;


  //Double_t fidVolMinX =    0; //Fiducial Volume dimensions for MicroBooNE
  //Double_t fidVolMaxX =  256;
  //Double_t fidVolMinY = -116;
  //Double_t fidVolMaxY =  116;
  //Double_t fidVolMinZ =    0;
  //Double_t fidVolMaxZ = 1030;


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
  fMatch_tag = p.get<std::string>("match_tag");
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

  art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_handle, e, fcluster_tag);

  // art::FindManyP<recob::Hit> sps_hit_assn_v(spacepoint_handle, e, fsps_tag);

  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_handle(hit_handle,e,fMatch_tag);


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
  Event_sps_ADC_U=0,Event_sps_ADC_V=0,Event_sps_ADC_Y=0,sps_ADC50=0,sps_ADC20=0,sps_ADC10=0;





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


  auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCProducer_tag);

  N_sps=0;

  for(size_t s=0;s<spacepoint_handle->size();s++){ //START SPACEPOINT LOOP
    N_sps++;



    matched_sps_reco_ADC= -9999.0;
    matched_sps_true_energy= -9999.0;

    std::cout<<"SpacePoint_Number: "<<s<<std::endl;


    auto sps = spacepoint_handle->at(s);
    auto cluster_v = sps_clus_assn_v.at(s);


    sps_x=sps.XYZ()[0];
    sps_y=sps.XYZ()[1];
    sps_z=sps.XYZ()[2];


    // std::cout<<"SpacePoint_Position X: "<<sps_x<<", Y:"<<sps_y<<", Z:"<<sps_z<<std::endl;
    //
    //
    // std::cout<<"SpacePoint Vector: "<<sps<<std::endl;

    for (auto const& cluster: cluster_v){
      auto plane = cluster->View();

      // std::cout<<"Cluster Plane: "<<cluster->View()<<std::endl;
      // std::cout<<"Cluster Vector Hit Size: "<<cluster->NHits()<<std::endl;
      // std::cout<<"Cluster Vector Summed ADC: "<<cluster->SummedADC()<<std::endl;






      if (plane==2){

        sps_ADC_Y = cluster->Integral();
        std::cout<<"sps_ADC Plane 2: "<<sps_ADC_Y<<std::endl;

      }

      if (plane==1){

        sps_ADC_V = cluster->Integral();
        std::cout<<"sps_ADC Plane 1: "<<sps_ADC_V<<std::endl;

      }

      if (plane==0){

        sps_ADC_U = cluster->Integral();
        std::cout<<"sps_ADC Plane 0: "<<sps_ADC_U<<std::endl;

      }

    }




    const std::vector<art::Ptr<recob::Hit>> hit_v = clus_hit_assn_v.at(s);

    hit_counter=0;
    hit_ADC_U_sum=0.0;
    hit_ADC_U=0.0;
    numElectrons_sum=0.0;

    for (art::Ptr<recob::Hit> hit : hit_v){//START HIT FOR LOOP

      hit_counter++;
      std::cout << "hit_counter: " <<hit_counter<< '\n';
      numElectrons=0.0;
      if (hit->View()==2) {
        std::cout << "Hit ADC Plane 2: " <<hit->Integral()<< '\n';

      }
      //
      if (hit->View()==1) {
        std::cout << "Hit ADC Plane 1: " <<hit->Integral()<< '\n';
      }
      //
      if (hit->View()==0) {
        std::cout << "Hit ADC Plane 0: " <<hit->Integral()<< '\n';
      }
      //
      // // if (hit->View() !=2)
      // //
      // // continue;



      hit_ADC_U= hit->Integral();
      hit_ADC_U_sum += hit->Integral();


      std::cout << "Hit Integral: " <<hit_ADC_U<< '\n';
      std::cout << "Hit ADC Sum: " <<hit_ADC_U_sum<< '\n';

      auto hitidx = hit.key();


      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      backtrack_handle.get(hitidx, particle_vec, match_vec);


      std::cout<<"MC Particle Size: "<<particle_vec.size()<<std::endl;

      bktrkd_particle_energy_per_sps=0.0;
      bktrkd_particle_energy=0.0;


      for (size_t p = 0; p < particle_vec.size(); p++){//START MC PARTICLE FOR LOOP FOR BACKTRACKER

        std::cout<<"MC Particle# "<<p<<std::endl;






        if (particle_vec.at(p)->PdgCode() != 11 || match_vec[p]->isMaxIDE!=1 || particle_vec.at(p)->Process() != "primary")
        continue;

        // if (particle_vec.at(p)->PdgCode() != 11)
        // continue;

        // if (match_vec[p]->isMaxIDE!=1 )
        // continue;


        std::cout<<"MC Particle# Backtracked "<<p<<std::endl;

        // std::cout<<"isMaxIDE for MC Particle "<<p<<" :"<<match_vec[p]->isMaxIDE<<std::endl;

        bktrkd_particle_energy=particle_vec.at(p)->E()* 1000.0;
        numElectrons=match_vec[p]->numElectrons;
        numElectrons_sum += match_vec[p]->numElectrons;
        std::cout << "bktrkd_particle_energy in the particle loop: " <<bktrkd_particle_energy<< '\n';
        std::cout<<"Number of electrons bktrkd: "<<numElectrons<<std::endl;
        std::cout<<"Number of electrons bktrkd (numElectrons_sum): "<<numElectrons_sum<<std::endl;

        bktrkd_particle_energy_per_sps +=bktrkd_particle_energy;
      }//END MC PARTICLE LOOP FOR BACKTRACKER


      Hittree->Fill();
    }//END HIT FOR LOOP





    std::cout<<"sps_ADC_Y: "<<sps_ADC_Y<<std::endl;
    // std::cout<<"matched_sps_reco_ADC: "<<sps_ADC<<std::endl;

    std::cout << "bktrkd_particle_energy in the sps loop: " <<bktrkd_particle_energy_per_sps<< '\n';

    Sps_Correlationtree->Fill();

    Event_sps_ADC_U +=sps_ADC_U;

    Event_sps_ADC_V +=sps_ADC_V;
    Event_sps_ADC_Y +=sps_ADC_Y;
  }//END SPACEPOINT LOOP







  Event_Energy_total=0.0;
  Event_Energy_total_pr=0.0;

  N_total_pr_e=0;
  N_total_e=0;


  for (size_t i_p = 0; i_p < mcp_h->size(); i_p++){//START MC PARTICLE FOR LOOP
    auto mcp = mcp_h->at(i_p);

    // std::cout << "MCParticle #: " <<i_p<< '\n';


    N_total_e++;


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

    MC_Particle_Energy= mcp.E()*1000;

    process= mcp.Process();



    std::cout<<"MC Particle Energy: "<<MC_Particle_Energy<<std::endl;

    Event_Energy_total += MC_Particle_Energy;
    std::cout<<"Process: "<<process<<std::endl;
    std::cout<<"Pdg: "<<pdg_particle<<std::endl;




    if(mcp.PdgCode() == 11 && mcp.Process() == "primary" ){
      N_total_pr_e++;


      StatusCode=mcp.StatusCode();
      position_X_pr=-9999.0;
      position_Y_pr=-9999.0;
      position_Z_pr=-9999.0;
      pdg_particle_pr=-9999;
      MC_Particle_Energy_pr= -9999.0;

      pdg_particle_pr=mcp.PdgCode();


      position_X_pr =mcp.Vx();
      position_Y_pr =mcp.Vy();
      position_Z_pr =mcp.Vz();

      MC_Particle_Energy_pr= mcp.E()*1000;

      endprocess= mcp.EndProcess();

      std::cout<<"Primary MC Particle Energy: "<<MC_Particle_Energy_pr<<std::endl;
      std::cout<<"Primary MC endprocess: "<<endprocess<<std::endl;

      Event_Energy_total_pr += MC_Particle_Energy_pr;

    }

    MCParticletree->Fill();


  }//END MC PARTICLE LOOP






  std::cout<<"SPS ADC per Event# "<<Event_sps_ADC_Y<<std::endl;
  std::cout<<"MC Particle Energy MeV per Event: "<<Event_Energy_total<<std::endl;
  std::cout<<"Primary MC Particle Energy MeV per Event: "<<Event_Energy_total_pr<<std::endl;
  std::cout<<"Number of primary electrons per Event: "<<N_total_pr_e<<std::endl;

  Eventtree->Fill();



}//END OF EVENT LOOP

void TruthStudies::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;

  Eventtree = tfs->make<TTree>("Eventtree",    "Eventtree");
  MCParticletree = tfs->make<TTree>("MCParticletree",    "MCParticletree");
  Sps_Correlationtree = tfs->make<TTree>("Sps_Correlationtree",    "Sps_Correlationtree");
  potTree = tfs->make<TTree>("potTree","potTree");
  Hittree = tfs->make<TTree>("Hittree",    "Hittree");



  Eventtree->Branch("Event_Energy_total_MeV",&Event_Energy_total,"Event_Energy_total_MeV/D");
  Eventtree->Branch("Event_Energy_total_pr",&Event_Energy_total_pr,"Event_Energy_total_pr/D");
  // Eventtree->Branch("Event_Energy_5_MeV",&Event_Energy_5_MeV,"Event_Energy_5_MeV/D");
  // Eventtree->Branch("Event_Energy_10_MeV",&Event_Energy_10_MeV,"Event_Energy_10_MeV/D");
  // Eventtree->Branch("Event_Energy_20_MeV",&Event_Energy_20_MeV,"Event_Energy_20_MeV/D");
  // Eventtree->Branch("Event_Energy_30_MeV",&Event_Energy_30_MeV,"Event_Energy_30_MeV/D");
  // Eventtree->Branch("Event_Energy_40_MeV",&Event_Energy_40_MeV,"Event_Energy_40_MeV/D");
  // Eventtree->Branch("Event_Energy_50_MeV",&Event_Energy_50_MeV,"Event_Energy_50_MeV/D");
  // Eventtree->Branch("Event_Energy_100_MeV",&Event_Energy_100_MeV,"Event_Energy_100_MeV/D");
  Eventtree->Branch("evttime",&evttime,"evttime/I");
  // Eventtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  // Eventtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  // Eventtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  // Eventtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  // Eventtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  // Eventtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  // Eventtree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
  Eventtree->Branch("N_sps",&N_sps,"N_sps/I");
  Eventtree->Branch("Event_sps_ADC_U",&Event_sps_ADC_U,"Event_sps_ADC_U/D");
  Eventtree->Branch("Event_sps_ADC_V",&Event_sps_ADC_V,"Event_sps_ADC_V/D");
  Eventtree->Branch("Event_sps_ADC_Y",&Event_sps_ADC_Y,"Event_sps_ADC_Y/D");
  // Eventtree->Branch("distance_smallest_rand_vtx",&distance_smallest_rand_vtx,"distance_smallest_rand_vtx/D");
  // Eventtree->Branch("N_sps10",&N_sps10,"N_sps10/I");
  // Eventtree->Branch("N_sps20",&N_sps20,"N_sps20/I");
  // Eventtree->Branch("N_sps50",&N_sps50,"N_sps50/I");
  // Eventtree->Branch("sps_ADC10",&sps_ADC10,"sps_ADC10/D");
  // Eventtree->Branch("sps_ADC20",&sps_ADC20,"sps_ADC20/D");
  // Eventtree->Branch("sps_ADC50",&sps_ADC50,"sps_ADC50/D");
  Eventtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Eventtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Eventtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Eventtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Eventtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Eventtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Eventtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  // Eventtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  // Eventtree->Branch("distance_nu_smallest" ,&distance_nu_smallest ,"distance_nu_smallest/D" );
  // Eventtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  Eventtree->Branch("N_Run" ,&N_Run ,"N_Run/I" );
  Eventtree->Branch("N_SubRun" ,&N_SubRun ,"N_SubRun/I" );
  Eventtree->Branch("N_total_e" ,&N_total_e ,"N_total_e/I" );
  Eventtree->Branch("N_total_pr_e" ,&N_total_pr_e ,"N_total_pr_e/I" );
  Eventtree->Branch("endprocess" ,&endprocess );


  MCParticletree->Branch("MC_Particle_Energy",&MC_Particle_Energy,"MC_Particle_Energy/D");
  MCParticletree->Branch("position_X",&position_X,"position_X/D");
  MCParticletree->Branch("position_Y",&position_Y,"position_Y/D");
  MCParticletree->Branch("position_Z",&position_Z,"position_Z/D");
  MCParticletree->Branch("pdg_particle",&pdg_particle,"pdg_particle/I");
  MCParticletree->Branch("endprocess" ,&endprocess );
  MCParticletree->Branch("process" ,&process );



  Sps_Correlationtree->Branch("evttime",&evttime,"evttime/I");
  Sps_Correlationtree->Branch("sps_x",&sps_x,"sps_x/D");
  Sps_Correlationtree->Branch("sps_y",&sps_y,"sps_y/D");
  Sps_Correlationtree->Branch("sps_z",&sps_z,"sps_z/D");
  // Sps_Correlationtree->Branch("distance",&distance,"distance/D");
  Sps_Correlationtree->Branch("sps_ADC_U",&sps_ADC_U,"sps_ADC_U/D");
  Sps_Correlationtree->Branch("sps_ADC_V",&sps_ADC_V,"sps_ADC_V/D");
  Sps_Correlationtree->Branch("sps_ADC_Y",&sps_ADC_Y,"sps_ADC_Y/D");
  Sps_Correlationtree->Branch("N_Event",&N_Event,"N_Event/I");
  Sps_Correlationtree->Branch("N_Run",&N_Run,"N_Run/I");
  Sps_Correlationtree->Branch("N_SubRun",&N_SubRun,"N_SubRun/I");
  // Sps_Correlationtree->Branch("Vertex_x",&Vertex_x,"Vertex_x/D");
  // Sps_Correlationtree->Branch("Vertex_y",&Vertex_y,"Vertex_y/D");
  // Sps_Correlationtree->Branch("Vertex_z",&Vertex_z,"Vertex_z/D");
  // Sps_Correlationtree->Branch("distance_rand_vtx",&distance_rand_vtx,"distance_rand_vtx/D");
  // Sps_Correlationtree->Branch("_rand_vtx_x",&_rand_vtx_x,"_rand_vtx_x/D");
  // Sps_Correlationtree->Branch("_rand_vtx_y",&_rand_vtx_y,"_rand_vtx_y/D");
  // Sps_Correlationtree->Branch("_rand_vtx_z",&_rand_vtx_z,"_rand_vtx_z/D");
  Sps_Correlationtree->Branch("neutrinos",&neutrinos,"neutrinos/I");
  Sps_Correlationtree->Branch("neutrinoshowers",&_neutrinoshowers,"neutrinoshowers/I");
  Sps_Correlationtree->Branch("neutrinotracks" ,&_neutrinotracks ,"neutrinotracks/I" );
  Sps_Correlationtree->Branch("muon_px" ,&_muon_px ,"muon_px/F" );
  Sps_Correlationtree->Branch("muon_py" ,&_muon_py ,"muon_py/F" );
  Sps_Correlationtree->Branch("muon_pz" ,&_muon_pz ,"muon_pz/F" );
  Sps_Correlationtree->Branch("maxTrkLen" ,&_maxTrkLen ,"maxTrkLen/F" );
  // Sps_Correlationtree->Branch("tracklength" ,&tracklength ,"tracklength/D" );
  // Sps_Correlationtree->Branch("distance_nu_smallest" ,&distance_nu_smallest ,"distance_nu_smallest/D" );
  // Sps_Correlationtree->Branch("cosmic_trk_50" ,&cosmic_trk_50 ,"cosmic_trk_50/I" );
  // Sps_Correlationtree->Branch("pointdistance_trk_smallest" ,&pointdistance_trk_smallest ,"pointdistance_trk_smallest/D" );
  // Sps_Correlationtree->Branch("distance_trk_smallest" ,&distance_trk_smallest ,"distance_trk_smallest/D" );
  Sps_Correlationtree->Branch("bktrkd_particle_energy_per_sps",&bktrkd_particle_energy_per_sps,"bktrkd_particle_energy_per_sps/D");


  Hittree->Branch("numElectrons" ,&numElectrons ,"numElectrons/D" );
  Hittree->Branch("numElectrons_sum" ,&numElectrons_sum ,"numElectrons_sum/D" );
  Hittree->Branch("hit_ADC_U" ,&hit_ADC_U ,"hit_ADC_U/D" );
  Hittree->Branch("hit_ADC_U_sum" ,&hit_ADC_U_sum ,"hit_ADC_U_sum/D" );
  Hittree->Branch("hit_counter" ,&hit_counter ,"hit_counter/I" );


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
