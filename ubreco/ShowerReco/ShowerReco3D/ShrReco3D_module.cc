///////////////////////////////////////////////////////////////////////
// Class:       Showerreco
// Plugin Type: producer (art v2_09_06)
// File:        ShrReco3D_module.cc
//
// Generated at Fri Feb  9 16:38:52 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// shower-reco classes and utilities
#include "ProtoShower/ProtoShowerAlgBase.h"
#include "Base/ShrRecoManager.h"
// include specific protoshower and recomanager instances

#include "art_root_io/TFileService.h"

#include "art/Utilities/make_tool.h"

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()

#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include <Eigen/Dense>

#include "art/Persistency/Common/PtrMaker.h"

#include "TTree.h"

#include "Base/ShowerRecoException.h"

#include <memory>

class ShrReco3D;


class ShrReco3D : public art::EDProducer {
public:
  explicit ShrReco3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  ShrReco3D(ShrReco3D const &) = delete;
  ShrReco3D(ShrReco3D &&) = delete;
  ShrReco3D & operator = (ShrReco3D const &) = delete;
  ShrReco3D & operator = (ShrReco3D &&) = delete;
  
  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  
private:
  
  /// Input producer name
  std::string fPFPproducer;
  std::string fClusproducer;
  std::string fVtxproducer;
  std::string fBacktrackTag;
  // is this event neutrino or single particle?
  bool fNeutrinoEvent;
  // fill a ttree?
  bool fFillTree;
  
  /// map for backtracking which stores mcshower index to vector of track ids for the mcshower
  std::map<size_t, std::vector<unsigned int> > _MCShowerInfo;  
  
  /// Shower reco core class instance
  ::showerreco::ShrRecoManager* _manager;
  
  // ProtoShowerAlgBase to make protoshowers
  std::unique_ptr<::protoshower::ProtoShowerAlgBase> _psalg;
  //::protoshower::ProtoShowerAlgBase* _psalg;

  trkf::SpacePointAlg* spsalg;

  // reconstruction ttree
  TTree* _rcshr_tree;
  double _shr_px, _shr_py, _shr_pz;
  double _shr_x, _shr_y, _shr_z;
  std::vector<double> _shr_dedx_v, _shr_dedx_pl0_v, _shr_dedx_pl1_v, _shr_dedx_pl2_v;
  std::vector<double> _shr_e_v;
  double _mc_shr_e;
  double _mc_shr_px, _mc_shr_py, _mc_shr_pz;
  double _mc_shr_x, _mc_shr_y, _mc_shr_z;
  double _xtimeoffset;
  double _xsceoffset;
  double _completeness, _purity;
  int _mc_shr_pdg;

  /**
     @brief Save output showers produced by reconstruction algorithms
     #input art::Event for backtracking
     @input s : index in output shower vector, to associate back to PFP index
     @input shower : reconstructed shower to be translated to recob:: object
     @input Shower_v : vector of recob::Showers. This will be placed in the art::Event
     @input pfp_h : handle to input pfparticles
     @input pfp_clus_assn_v : input pfp -> clus assn vector
     @input pfp_hit_assn_v  : input pfp -> hit assn vector
     @input Shower_PFP_assn_v     : output shower -> pfp assn vector
     @input Shower_Cluster_assn_v : output shower -> cluster assn vector
     @input Shower_Hit_assn_v     : output shower -> hit assn vector
  */
  void SaveShower(detinfo::DetectorClocksData const& detClocks,
                  detinfo::DetectorPropertiesData const& detProperties,
                  art::Event & e,
		  const size_t idx,
		  const showerreco::Shower_t& shower,
		  std::unique_ptr< std::vector<recob::Shower> >& Shower_v,
		  const art::PtrMaker<recob::Shower> ShowerPtrMaker,
 		  std::unique_ptr< std::vector<recob::SpacePoint> >& SpacePoint_v,
 		  const art::PtrMaker<recob::SpacePoint> SpacePointPtrMaker,
 		  std::unique_ptr< std::vector<recob::PCAxis> >& PCAxis_v,
 		  const art::PtrMaker<recob::PCAxis> PCAxisPtrMaker,
 		  std::unique_ptr< std::vector<recob::Vertex> >& Vertex_v,
 		  const art::PtrMaker<recob::Vertex> VertexPtrMaker,
		  const art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
		  const art::FindManyP<recob::Cluster> pfp_clus_assn_v,
		  const std::vector<std::vector<art::Ptr<recob::Hit> > > pfp_hit_assn_v,
		  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> >& Shower_PFP_assn_v,
		  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    >& Shower_Cluster_assn_v,
 		  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        >& Shower_Hit_assn_v,
 		  std::unique_ptr< art::Assns <recob::Shower, recob::SpacePoint> >& Shower_Sps_assn_v,
 		  std::unique_ptr< art::Assns <recob::Shower, recob::PCAxis>     >& Shower_PCAxis_assn_v,
 		  std::unique_ptr< art::Assns <recob::Shower, recob::Vertex>     >& Shower_Vertex_assn_v,
 		  std::unique_ptr< art::Assns <recob::SpacePoint, recob::Hit>    >& Sps_Hit_assn_v );

  /**
     @brief Backtracking function to find best MCShower matched with hits of a reco shower
     This function loops thorugh mcshowers, and for each computes the purity and completeness based on the backtarcked energy associated to each hit in the shower. The entry with the best purity is determined as the best match and returned.
     @input art::Event with which to load backtracking metadata
     @input hit_idx_v : vector of hit indices (within event handle) associated to the shower (can be from a single plane)
     @input completeness : value passed by reference matched to the most pure MCShower match
     @input purity : see above
     @return index of associated mcshower
   */
  size_t BackTrack(detinfo::DetectorClocksData const& detClocks,
                   detinfo::DetectorPropertiesData const& detProperties,
                   art::Event & e, const std::vector<unsigned int>& hit_idx_v, float& completeness, float& purity);

  /**
     @brief find showers associated to neutrino interaction
     @input mct_h mctruth
     @input mcs_h mcshowers
     @ return map connecting mcshower index to vector of track IDs associated to that shower
   */
  std::map<size_t, std::vector<unsigned int> > GetMCShowerInfo(const art::ValidHandle<std::vector<simb::MCTruth> > mct_h, const art::Handle<std::vector<sim::MCShower> > mcs_h);

  /**
     @brief return list of mcshowers and their associated track IDs
     @input mcs_h mcshowers
     @ return map connecting mcshower index to vector of track IDs associated to that shower
   */
  std::map<size_t, std::vector<unsigned int> > GetMCShowerInfo(const art::Handle<std::vector<sim::MCShower> > mcs_h);
  
  recob::PCAxis ComputePCA(std::vector<recob::SpacePoint> spv);
  // Declare member data here.

  /**
     @brief Set TTree to be used to save reconstructed shower variables
   */
  void SetTTree();
  
};


ShrReco3D::ShrReco3D(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{

  fPFPproducer   = p.get<std::string>("PFPproducer"  );
  fClusproducer  = p.get<std::string>("Clusproducer" );
  fVtxproducer   = p.get<std::string>("Vtxproducer"  );
  fBacktrackTag  = p.get<std::string>("BacktrackTag","" );
  fNeutrinoEvent = p.get<bool>       ("NeutrinoEvent");
  fFillTree      = p.get<bool>       ("FillTree",false);
  
  const fhicl::ParameterSet& protoshower_pset = p.get<fhicl::ParameterSet>("ProtoShowerTool");  

  // grab algorithms for merging
  _manager = new showerreco::ShrRecoManager();
  _manager->Clear();
  const fhicl::ParameterSet& showerrecoTools = p.get<fhicl::ParameterSet>("ShowerRecoTools");
  for (const std::string& showerrecoTool : showerrecoTools.get_pset_names()) {
    const fhicl::ParameterSet& showerreco_pset = showerrecoTools.get<fhicl::ParameterSet>(showerrecoTool);
    _manager->AddAlgo(art::make_tool<showerreco::ShowerRecoModuleBase>(showerreco_pset));
  }// for all algorithms to be added

  _manager->SetDebug(false);

  //_manager = new ::showerreco::Pi0RecoAlgorithm();
  _psalg = art::make_tool<::protoshower::ProtoShowerAlgBase>(protoshower_pset);

  spsalg = new trkf::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));

  produces<std::vector<recob::Shower> >();
  produces<std::vector<recob::SpacePoint> >();
  produces<std::vector<recob::PCAxis> >();
  produces<std::vector<recob::Vertex> >();
  produces<art::Assns <recob::Shower, recob::PFParticle> >();
  produces<art::Assns <recob::Shower, recob::Cluster>    >();
  produces<art::Assns <recob::Shower, recob::Hit>        >();
  produces<art::Assns <recob::Shower, recob::SpacePoint> >();
  produces<art::Assns <recob::Shower, recob::PCAxis> >();
  produces<art::Assns <recob::Shower, recob::Vertex> >();
  produces<art::Assns <recob::SpacePoint, recob::Hit> >();

  _manager->Initialize();

  auto recomb = p.get<double>("recombination");
  auto adctoe = p.get<std::vector<double> >("ADCtoE");
  if (adctoe.size() != 3) {
    std::cout << "ERROR provided !3 planes for calorimetry" << std::endl;
  }

  std::vector<double> calib = {adctoe[0] * 0.0000236 / recomb,
			       adctoe[1] * 0.0000236 / recomb,
			       adctoe[2] * 0.0000236 / recomb };

  _psalg->setCalorimetry(calib);

  SetTTree();
  
}

void ShrReco3D::produce(art::Event & e)
{

  // produce recob::Showers
  std::unique_ptr< std::vector<recob::Shower> > Shower_v(new std::vector<recob::Shower> );
  std::unique_ptr< std::vector<recob::SpacePoint> > SpacePoint_v(new std::vector<recob::SpacePoint> );
  std::unique_ptr< std::vector<recob::PCAxis> > PCAxis_v(new std::vector<recob::PCAxis> );
  std::unique_ptr< std::vector<recob::Vertex> > Vertex_v(new std::vector<recob::Vertex> );
  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> > Shower_PFP_assn_v    ( new art::Assns<recob::Shower, recob::PFParticle>);
  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    > Shower_Cluster_assn_v( new art::Assns<recob::Shower, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        > Shower_Hit_assn_v    ( new art::Assns<recob::Shower, recob::Hit>       );
  std::unique_ptr< art::Assns <recob::Shower, recob::SpacePoint> > Shower_Sps_assn_v    ( new art::Assns<recob::Shower, recob::SpacePoint>);
  std::unique_ptr< art::Assns <recob::Shower, recob::PCAxis>     > Shower_PCAxis_assn_v ( new art::Assns<recob::Shower, recob::PCAxis>    );
  std::unique_ptr< art::Assns <recob::Shower, recob::Vertex>     > Shower_Vertex_assn_v ( new art::Assns<recob::Shower, recob::Vertex>    );
  std::unique_ptr< art::Assns <recob::SpacePoint, recob::Hit>    > Sps_Hit_assn_v       ( new art::Assns<recob::SpacePoint, recob::Hit>   );

  // shower pointer maker for later to create associations
  art::PtrMaker<recob::Shower> ShowerPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> SpacePointPtrMaker(e);
  art::PtrMaker<recob::PCAxis> PCAxisPtrMaker(e);
  art::PtrMaker<recob::Vertex> VertexPtrMaker(e);

  // pass event to ProtoShowerAlgBase to create ProtoShower objects
  // which will then be fed to shower reco algorithm chain
  std::vector<protoshower::ProtoShower> event_protoshower_v;
  _psalg->GenerateProtoShowers(e, event_protoshower_v);

  // set protoshowers for algorithms
  _manager->SetProtoShowers(event_protoshower_v);

  // output showers to be saved to event
  std::vector< ::showerreco::Shower_t> output_shower_v;
  auto const detClocks = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProperties = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, detClocks);
  util::GeometryUtilities const gser{*lar::providerFrom<geo::Geometry>(),
                                     art::ServiceHandle<geo::WireReadout const>()->Get(),
                                     detClocks,
                                     detProperties};
  _manager->Reconstruct(gser, output_shower_v);

  // if using truth, backtrack and load mcshowers
  if (fBacktrackTag != ""){
    // load mcshowers & mctruth
    art::Handle< std::vector<sim::MCShower> > mcs_h;
    e.getByLabel("mcreco",mcs_h);
    auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
    if (!mcs_h.isValid()) {
      _MCShowerInfo.clear();
    }
    else {
      if (fNeutrinoEvent)
	_MCShowerInfo = GetMCShowerInfo(mct_h,mcs_h);
      else
	_MCShowerInfo = GetMCShowerInfo(mcs_h);
    }
  }

  // load PFP, clus, hit so that associations to showers can be stored
  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  // grab clusters associated with PFParticles
  art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);
  // ADDITION FROM PETRILLO
  e.getValidHandle<std::vector<recob::Cluster>>(fClusproducer);
  // grab the hits associated to the PFParticles
  auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

  // save output showers
  for (size_t s=0; s < output_shower_v.size(); s++) {
    try {
      SaveShower(detClocks, detProperties,
                 e, s, output_shower_v.at(s), Shower_v, ShowerPtrMaker,
		 SpacePoint_v, SpacePointPtrMaker, PCAxis_v, PCAxisPtrMaker, Vertex_v, VertexPtrMaker,
		 pfp_h, pfp_clus_assn_v, pfp_hit_assn_v,
		 Shower_PFP_assn_v, Shower_Cluster_assn_v, Shower_Hit_assn_v,
		 Shower_Sps_assn_v, Shower_PCAxis_assn_v, Shower_Vertex_assn_v, Sps_Hit_assn_v);
    }
    catch (showerreco::ShowerRecoException const& e) {
      std::cout << e.what() << std::endl;
    }
  }// for all output reconstructed showers
  
  e.put(std::move(Shower_v));
  e.put(std::move(SpacePoint_v));
  e.put(std::move(PCAxis_v));
  e.put(std::move(Vertex_v));
  e.put(std::move(Shower_PFP_assn_v));
  e.put(std::move(Shower_Cluster_assn_v));
  e.put(std::move(Shower_Hit_assn_v));
  e.put(std::move(Shower_Sps_assn_v));
  e.put(std::move(Shower_PCAxis_assn_v));
  e.put(std::move(Shower_Vertex_assn_v));
  e.put(std::move(Sps_Hit_assn_v));

}

void ShrReco3D::beginJob()
{
  // Implementation of optional member function here.
}

void ShrReco3D::endJob()
{
  _manager->Finalize();
}


void ShrReco3D::SaveShower(detinfo::DetectorClocksData const& detClocks,
                           detinfo::DetectorPropertiesData const& detProperties,
                           art::Event & e,
			   const size_t idx,
			   const showerreco::Shower_t& shower,
			   std::unique_ptr< std::vector<recob::Shower> >& Shower_v,
			   const art::PtrMaker<recob::Shower> ShowerPtrMaker,
                           std::unique_ptr< std::vector<recob::SpacePoint> >& SpacePoint_v,
                           const art::PtrMaker<recob::SpacePoint> SpacePointPtrMaker,
                           std::unique_ptr< std::vector<recob::PCAxis> >& PCAxis_v,
                           const art::PtrMaker<recob::PCAxis> PCAxisPtrMaker,
                           std::unique_ptr< std::vector<recob::Vertex> >& Vertex_v,
                           const art::PtrMaker<recob::Vertex> VertexPtrMaker,
			   const art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
			   const art::FindManyP<recob::Cluster> pfp_clus_assn_v,
			   const std::vector<std::vector<art::Ptr<recob::Hit> > > pfp_hit_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> >& Shower_PFP_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::Cluster> >& Shower_Cluster_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::Hit> >& Shower_Hit_assn_v,
                           std::unique_ptr< art::Assns <recob::Shower, recob::SpacePoint> >& Shower_Sps_assn_v,
                           std::unique_ptr< art::Assns <recob::Shower, recob::PCAxis> >& Shower_PCAxis_assn_v,
                           std::unique_ptr< art::Assns <recob::Shower, recob::Vertex> >& Shower_Vertex_assn_v,
                           std::unique_ptr< art::Assns <recob::SpacePoint, recob::Hit> >& Sps_Hit_assn_v)
  
{

  if (shower.fPassedReconstruction == false) {
    std::stringstream ss;
    ss << "Shower " << idx << " did not pass reconstruction.";
    throw showerreco::ShowerRecoException(ss.str());
  }
  
  // filter out showers with garbage values
  if (shower.fXYZStart.Mag2()  == 0) {
    std::stringstream ss;
    ss << "Shower " << idx << " does not have a valid XYZ start point";
    throw showerreco::ShowerRecoException(ss.str());
  }
  if (shower.fDCosStart.Mag2() == 0) {
    std::stringstream ss;
    ss << "Shower " << idx << " does not have a valid 3D direction";
    throw showerreco::ShowerRecoException(ss.str());
  }

  // make sure associations exist for PFP, clusters, and hits
  if (pfp_h->size() <= shower.fIndex) {
    std::stringstream ss;
    ss << "Shower " << idx << " does not have PFP associated to index " << shower.fIndex;
    throw showerreco::ShowerRecoException(ss.str());
  }
  if (pfp_clus_assn_v.size() <= shower.fIndex) {
    std::stringstream ss;
    ss << "Shower " << idx << " does not have CLUSTERS associated to index " << shower.fIndex;
    throw showerreco::ShowerRecoException(ss.str());
  }
  if (pfp_hit_assn_v.size() <= shower.fIndex) {
    std::stringstream ss;
    ss << "Shower " << idx << " does not have HITS associated to index " << shower.fIndex;
    throw showerreco::ShowerRecoException(ss.str());
  }

  // save shower variables to TTree
  _shr_dedx_pl0_v = shower.fdEdx_v_v[0];
  _shr_dedx_pl1_v = shower.fdEdx_v_v[1];
  _shr_dedx_pl2_v = shower.fdEdx_v_v[2];
  _shr_dedx_v     = shower.fdEdx_v;
  _shr_px     = shower.fDCosStart[0];
  _shr_py     = shower.fDCosStart[1];
  _shr_pz     = shower.fDCosStart[2];
  _shr_x      = shower.fXYZStart[0];
  _shr_y      = shower.fXYZStart[1];
  _shr_z      = shower.fXYZStart[2];
  _shr_e_v    = shower.fTotalEnergy_v;


  recob::Shower s;
  s.set_id ( Shower_v->size() );
  s.set_total_energy          ( shower.fTotalEnergy_v         );
  s.set_total_energy_err      ( shower.fSigmaTotalEnergy_v    );
  s.set_total_best_plane      ( shower.fBestPlane.Plane       );
  s.set_direction             ( shower.fDCosStart             );
  s.set_direction_err         ( shower.fSigmaDCosStart        );
  s.set_start_point           ( shower.fXYZStart              );
  s.set_start_point_err       ( shower.fSigmaXYZStart         );
  s.set_dedx                  ( shower.fdEdx_v                );
  s.set_dedx_err              ( shower.fSigmadEdx_v           );
  s.set_length                ( shower.fLength                );
  s.set_open_angle            ( shower.fOpeningAngle          );
  
  Shower_v->emplace_back(s);

  art::Ptr<recob::Shower> const ShrPtr = ShowerPtrMaker(Shower_v->size()-1);

  // now take care of associations
  
  // step 1 : pfp
  const art::Ptr<recob::PFParticle> PFPPtr(pfp_h, shower.fIndex);
  Shower_PFP_assn_v->addSingle( ShrPtr, PFPPtr );

  // step 2 : clusters
  std::vector<art::Ptr<recob::Cluster> > clus_v = pfp_clus_assn_v.at(shower.fIndex);
  for (size_t c=0; c < clus_v.size(); c++)
    Shower_Cluster_assn_v->addSingle( ShrPtr, clus_v.at(c) );

  // step 3 : hits
  std::vector<art::Ptr<recob::Hit> > hit_v = pfp_hit_assn_v.at(shower.fIndex);
  // for backtracking purposes save the vector of collection-plane hit indices associated to the shower
  std::vector<unsigned int> hit_idx_v;
  for (size_t h=0; h < hit_v.size(); h++) {
    Shower_Hit_assn_v->addSingle( ShrPtr, hit_v.at(h) );
    if (hit_v.at(h)->WireID().Plane == 2) { hit_idx_v.push_back( hit_v.at(h).key() ); }
  }// for all hits

  // step 4 : sps
  std::vector<recob::SpacePoint> spts;
  art::PtrVector<recob::Hit> hit_ptr_v;
  for (auto h : hit_v) hit_ptr_v.push_back(h);
  spsalg->makeSpacePoints(detClocks, detProperties, hit_ptr_v, spts);
  for (auto sp : spts) {
    art::PtrVector<recob::Hit> hits = spsalg->getAssociatedHits(sp);
    for (auto hp : hits) {
      //make sure we have one sps for each hit, as that is what Pandora does
      SpacePoint_v->push_back(sp);
      art::Ptr<recob::SpacePoint> const SpsPtr = SpacePointPtrMaker(SpacePoint_v->size()-1);
      Shower_Sps_assn_v->addSingle( ShrPtr, SpsPtr );
      //
      Sps_Hit_assn_v->addSingle( SpsPtr, hp );
    }
  }

  // step 5 : PCAxis
  PCAxis_v->push_back(ComputePCA(spts));
  art::Ptr<recob::PCAxis> const PcaPtr = PCAxisPtrMaker(PCAxis_v->size()-1);
  Shower_PCAxis_assn_v->addSingle( ShrPtr, PcaPtr );

  // step 6 : Vertex
  double vtx[3] = {s.ShowerStart().x(),s.ShowerStart().y(),s.ShowerStart().z()};
  Vertex_v->push_back(recob::Vertex(vtx,int(Shower_v->size()-1)));
  art::Ptr<recob::Vertex> const VtxPtr = VertexPtrMaker(Vertex_v->size()-1);
  Shower_Vertex_assn_v->addSingle( ShrPtr, VtxPtr );
  
  if (fBacktrackTag != ""){
    float completeness_max, purity_max;
    BackTrack(detClocks, detProperties, e, hit_idx_v, completeness_max, purity_max);
  }

  if (fFillTree)
    _rcshr_tree->Fill();
  
  return;
}// end of SaveShower function

size_t ShrReco3D::BackTrack(detinfo::DetectorClocksData const& detClocks,
                            detinfo::DetectorPropertiesData const& detProperties,
                            art::Event & e, const std::vector<unsigned int>& hit_idx_v, float& completeness_max, float& purity_max)
{
  
  art::InputTag BacktrackTag { fBacktrackTag };
  auto const& gaushit_h = e.getValidHandle<std::vector<recob::Hit> > ("gaushit");
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_h(gaushit_h,e,BacktrackTag);

  // score the match for each MCShower
  purity_max = 0.;
  completeness_max = 0.;
  size_t mcs_idx_match = 0;

  art::Handle< std::vector<sim::MCShower> > mcs_h;
  e.getByLabel("mcreco",mcs_h);
  if (!mcs_h.isValid()) {
    return mcs_idx_match;
  }

  for (auto const& mcshower : _MCShowerInfo) {

    size_t s = mcshower.first; // index of mcshower in mcs_h
    
    float purity = 0;
    float completeness = 0;

    // float BackTrackEnergy       = 0; // unused
    float BackTrackShowerEnergy = 0;
    float BackTrackCharge       = 0;
    float BackTrackShowerCharge = 0;

    auto const& mcs = mcs_h->at(s);
    auto shrtrackIDs = mcshower.second;
    
    //std::cout << "\t ANCESTOR comparing with MCShower of energy " << mcs.Start().E() << std::endl;
    //std::cout << "\t ANCESTOR start is [" << mcs.Start().X() << ", " << mcs.Start().Y() << ", " << mcs.Start().Z() << "]" << std::endl;
    //std::cout << "\t ANCESTOR HAS " << shrtrackIDs.size() << " particles" << std::endl;
    
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;  
    
    for (auto const& hit_idx : hit_idx_v) {
      
      particle_vec.clear(); match_vec.clear();

      backtrack_h.get(hit_idx, particle_vec, match_vec);
      
      // does this hit match to the mcshower?
      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){            
	
	auto mctrkid = particle_vec.at(i_p)->TrackId();
	auto charge  = match_vec[i_p]->numElectrons;
	auto energy  = match_vec[i_p]->energy;

	BackTrackCharge += charge;
	// BackTrackEnergy += energy; // unused
	// does this trackID match that of the MCShower?
	for (auto const& shrtrkid : shrtrackIDs) {
	  if ( shrtrkid == (unsigned int)mctrkid ){
	    BackTrackShowerCharge += charge;
	    BackTrackShowerEnergy += energy;
	    break;
	  }
	}
      }// for all particles associated to this hit
    }// for all hits
    
    purity       = BackTrackShowerCharge / BackTrackCharge;
    completeness = BackTrackShowerEnergy / mcs.Start().E();

    //std::cout << "ANCESTOR Purity : " << BackTrackShowerCharge << " / " << BackTrackCharge << " = " << purity << std::endl;
    
    if (purity > purity_max) {
      purity_max = purity;
      completeness_max = completeness;
      mcs_idx_match = s;
    }
    
  }// end of MCShower loop

  //std::cout << "ANCESTOR max purity : " << purity_max << std::endl;

  _completeness = completeness_max;
  _purity       = purity_max;
  auto matched_mcs = mcs_h->at(mcs_idx_match);
  _mc_shr_e = matched_mcs.Start().E();
  _mc_shr_pdg = matched_mcs.PdgCode();
  _mc_shr_x = matched_mcs.Start().X();
  _mc_shr_y = matched_mcs.Start().Y();
  _mc_shr_z = matched_mcs.Start().Z();

  // get X offset due to time w.r. trigger time
  if (fNeutrinoEvent) {
    auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
    auto gen = mct_h->at(0);
    double g4Ticks = detClocks.TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties.GetXTicksOffset(0, 0, 0) - trigger_offset(detClocks);
    //std::cout << "nu vtx @ [" << gen.GetNeutrino().Nu().Vx() << ", " << gen.GetNeutrino().Nu().Vy() << ", " << gen.GetNeutrino().Nu().Vz() << " ]" << std::endl;
    _xtimeoffset = detProperties.ConvertTicksToX(g4Ticks, 0, 0, 0);
  }
  else { _xtimeoffset = 0.; }

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(geo::Point_t(_mc_shr_x,_mc_shr_y,_mc_shr_z));
  //_mc_shr_x += offset.X() + xtrueoffset;
  _xsceoffset = offset.X();
  _mc_shr_y += offset.Y();
  _mc_shr_z += offset.Z();

  auto mom = matched_mcs.Start().Momentum().Vect().Mag();
  _mc_shr_px = matched_mcs.Start().Px() / mom;
  _mc_shr_py = matched_mcs.Start().Py() / mom;
  _mc_shr_pz = matched_mcs.Start().Pz() / mom;

  // fill TTree variables
  
  return mcs_idx_match;
}// function end


std::map<size_t, std::vector<unsigned int> > ShrReco3D::GetMCShowerInfo(const art::Handle<std::vector<sim::MCShower> > mcs_h) {

  std::map<size_t, std::vector<unsigned int> > event_shower_map;
  // map connecting e+/e- trackID in mcshower to mcshower index
  //std::map<unsigned int, size_t> event_mcpart_map;

  for (size_t i=0; i < mcs_h->size(); i++) {

    auto const& mcs = mcs_h->at(i);

    std::vector<unsigned int> shrtrackIDs = mcs.DaughterTrackID();
    shrtrackIDs.push_back( mcs.TrackID() );
    //std::cout << "\t\t shower track ID : " << mcs.TrackID() << std::endl; 
    // get daughter track IDs:
    auto daughterIDs = mcs.DaughterTrackID();
    for (auto const& id : daughterIDs)
      if (id != mcs.TrackID()) { shrtrackIDs.push_back(id); } //std::cout << "\t\t shower track ID : " << id << std::endl; }
    
    event_shower_map[ i ] = shrtrackIDs;
    
    //std::cout << "\t ANCESTOR mother PDG is " << mcs.MotherPdgCode() << std::endl;
    //std::cout << "\t ANCESTOR start is [" << mcs.Start().X() << ", " << mcs.Start().Y() << ", " << mcs.Start().Z() << "]" << std::endl;
    ///std::cout << "\t ANCESTOR Process is " << mcs.Process() << std::endl;
    //std::cout << "\t ANCESTOR energy is " << mcs.Start().E() << std::endl; 
    //std::cout << "\t ANCESTOR number of daughters is " << shrtrackIDs.size() << std::endl; 
    
  }// for all mcshowers
  
  return event_shower_map;
}

std::map<size_t, std::vector<unsigned int> > ShrReco3D::GetMCShowerInfo(const art::ValidHandle<std::vector<simb::MCTruth> > mct_h, const art::Handle<std::vector<sim::MCShower> > mcs_h) {
  
  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino().Nu();
  auto vtx = neutrino.Position(0).Vect();

  Double_t xyz[3] = {};
  xyz[0] = vtx.X();
  xyz[1] = vtx.Y();
  xyz[2] = vtx.Z();
  
  //std::cout << "ANCESTOR neutrino vertex @ [ " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;

  // loop through MCShowers and identify those originating from the pi0
  // map connecting mcshower index to track ID vector for all e+/e- in MCShower
  std::map<size_t, std::vector<unsigned int> > event_shower_map;
  // map connecting e+/e- trackID in mcshower to mcshower index
  //std::map<unsigned int, size_t> event_mcpart_map;

  for (size_t i=0; i < mcs_h->size(); i++) {
    auto const& mcs = mcs_h->at(i);

    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (xyz[0] - x) * (xyz[0] - x) ) +
		     ( (xyz[1] - y) * (xyz[1] - y) ) +
		     ( (xyz[2] - z) * (xyz[2] - z) ) );

    if ( (d < 0.01) ){// || ( (mcs.Process() == "primary") && (mcs.MotherPdgCode() == 22) ) ) {
      std::vector<unsigned int> shrtrackIDs = mcs.DaughterTrackID();
      shrtrackIDs.push_back( mcs.TrackID() );
      //std::cout << "\t\t shower track ID : " << mcs.TrackID() << std::endl; 
      // get daughter track IDs:
      auto daughterIDs = mcs.DaughterTrackID();
      for (auto const& id : daughterIDs)
	if (id != mcs.TrackID()) { shrtrackIDs.push_back(id); } //std::cout << "\t\t shower track ID : " << id << std::endl; }
      
      event_shower_map[ i ] = shrtrackIDs;

      //std::cout << "\t ANCESTOR mother PDG is " << mcs.MotherPdgCode() << std::endl;
      //std::cout << "\t ANCESTOR start is [" << mcs.Start().X() << ", " << mcs.Start().Y() << ", " << mcs.Start().Z() << "]" << std::endl;
      ///std::cout << "\t ANCESTOR Process is " << mcs.Process() << std::endl;
      //std::cout << "\t ANCESTOR energy is " << mcs.Start().E() << std::endl; 
      //std::cout << "\t ANCESTOR number of daughters is " << shrtrackIDs.size() << std::endl; 

    }// if mcshower matched to pi0
  }// for all mcshowers

  return event_shower_map;
}

recob::PCAxis ShrReco3D::ComputePCA(std::vector<recob::SpacePoint> sps){

  //float TotalCharge = 0;
  float sumWeights = 0;
  float xx = 0;
  float yy = 0;
  float zz = 0;
  float xy = 0;
  float xz = 0;
  float yz = 0;

  float ShowerCentre_c[3] = {0,0,0};
  for (auto const& sp : sps) {
    for (int i=0;i<3;i++) {
      ShowerCentre_c[i] += sp.XYZ()[i];
    }
  }
  for (int i=0;i<3;i++) ShowerCentre_c[i] *= (1. / sps.size());
  geo::Point_t ShowerCentre(ShowerCentre_c[0],ShowerCentre_c[1],ShowerCentre_c[2]);

  // //Get the Shower Centre
  // if (fChargeWeighted) {
  //   ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, sps, fmh, TotalCharge);
  // }
  // else {
  //   ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(sps);
  // }

  //Normalise the spacepoints, charge weight and add to the PCA.
  for (auto& sp : sps) {

    float wht = 1;

    //Normalise the spacepoint position.
    auto const sp_position = geo::Point_t(sp.XYZ()[0], sp.XYZ()[1], sp.XYZ()[2]) - ShowerCentre;

    // if (fChargeWeighted) {

    //   //Get the charge.
    //   float Charge = IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(sp, fmh);

    //   //Get the time of the spacepoint
    //   float Time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(sp, fmh);

    //   //Correct for the lifetime at the moment.
    //   Charge *= std::exp((sampling_rate(clockData) * Time) / (detProp.ElectronLifetime() * 1e3));

    //   //Charge Weight
    //   wht *= std::sqrt(Charge / TotalCharge);
    // }

    xx += sp_position.X() * sp_position.X() * wht;
    yy += sp_position.Y() * sp_position.Y() * wht;
    zz += sp_position.Z() * sp_position.Z() * wht;
    xy += sp_position.X() * sp_position.Y() * wht;
    xz += sp_position.X() * sp_position.Z() * wht;
    yz += sp_position.Y() * sp_position.Z() * wht;
    sumWeights += wht;
  }

  // Using Eigen package
  Eigen::Matrix3f matrix;

  // Construct covariance matrix
  matrix << xx, xy, xz, xy, yy, yz, xz, yz, zz;

  // Normalise from the sum of weights
  matrix /= sumWeights;

  // Run the PCA
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMatrix(matrix);

  Eigen::Vector3f eigenValuesVector = eigenMatrix.eigenvalues();
  Eigen::Matrix3f eigenVectorsMatrix = eigenMatrix.eigenvectors();

  // Put in the required form for a recob::PCAxis
  const bool svdOk = true; //TODO: Should probably think about this a bit more
  const int nHits = sps.size();
  // For some reason eigen sorts the eigenvalues from smallest to largest, reverse it
  const double eigenValues[3] = {
    eigenValuesVector(2), eigenValuesVector(1), eigenValuesVector(0)};
  std::vector<std::vector<double>> eigenVectors = {
    {eigenVectorsMatrix(0, 2), eigenVectorsMatrix(1, 2), eigenVectorsMatrix(2, 2)},
    {eigenVectorsMatrix(0, 1), eigenVectorsMatrix(1, 1), eigenVectorsMatrix(2, 1)},
    {eigenVectorsMatrix(0, 0), eigenVectorsMatrix(1, 0), eigenVectorsMatrix(2, 0)}};
  const double avePos[3] = {ShowerCentre.X(), ShowerCentre.Y(), ShowerCentre.Z()};

  return recob::PCAxis(svdOk, nHits, eigenValues, eigenVectors, avePos);
}

void ShrReco3D::SetTTree() {

  art::ServiceHandle<art::TFileService> tfs;

  // MC shower-by-shower TTree
  _rcshr_tree = tfs->make<TTree>("_rcshr_tree","ShrReco3D Shower TTree");
  _rcshr_tree->Branch("_shr_x",&_shr_x,"shr_x/D");
  _rcshr_tree->Branch("_xtimeoffset",&_xtimeoffset,"xtimeoffset/D");
  _rcshr_tree->Branch("_xsceoffset",&_xsceoffset,"xsceoffset/D");
  _rcshr_tree->Branch("_shr_y",&_shr_y,"shr_y/D");
  _rcshr_tree->Branch("_shr_z",&_shr_z,"shr_z/D");
  _rcshr_tree->Branch("_shr_dedx_pl0_v","std::vector<double>",&_shr_dedx_pl0_v);
  _rcshr_tree->Branch("_shr_dedx_pl1_v","std::vector<double>",&_shr_dedx_pl1_v);
  _rcshr_tree->Branch("_shr_dedx_pl2_v","std::vector<double>",&_shr_dedx_pl2_v);
  _rcshr_tree->Branch("_shr_e_v","std::vector<double>",&_shr_e_v);
  _rcshr_tree->Branch("_shr_dedx_v","std::vector<double>",&_shr_dedx_v);
  _rcshr_tree->Branch("_shr_px",&_shr_px,"shr_px/D");
  _rcshr_tree->Branch("_shr_py",&_shr_py,"shr_py/D");
  _rcshr_tree->Branch("_shr_pz",&_shr_pz,"shr_pz/D");
  // mc variables
  _rcshr_tree->Branch("_completeness",&_completeness,"completeness/D");
  _rcshr_tree->Branch("_purity",&_purity,"purity/D");
  _rcshr_tree->Branch("_mc_shr_pdg",&_mc_shr_pdg,"mc_shr_pdg/I");
  _rcshr_tree->Branch("_mc_shr_e",&_mc_shr_e,"mc_shr_e/D");
  _rcshr_tree->Branch("_mc_shr_x",&_mc_shr_x,"mc_shr_x/D");
  _rcshr_tree->Branch("_mc_shr_y",&_mc_shr_y,"mc_shr_y/D");
  _rcshr_tree->Branch("_mc_shr_z",&_mc_shr_z,"mc_shr_z/D");
  _rcshr_tree->Branch("_mc_shr_px",&_mc_shr_px,"mc_shr_px/D");
  _rcshr_tree->Branch("_mc_shr_py",&_mc_shr_py,"mc_shr_py/D");
  _rcshr_tree->Branch("_mc_shr_pz",&_mc_shr_pz,"mc_shr_pz/D");

}

DEFINE_ART_MODULE(ShrReco3D)
