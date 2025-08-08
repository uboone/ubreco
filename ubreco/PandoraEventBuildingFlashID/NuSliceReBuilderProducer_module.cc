////////////////////////////////////////////////////////////////////////
// Class:       NuSliceReBuilderProducer
// Plugin Type: producer (art v3_06_03)
// File:        NuSliceReBuilderProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"

class NuSliceReBuilderProducer;

class NuSliceReBuilderProducer : public art::EDProducer {
public:
  explicit NuSliceReBuilderProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuSliceReBuilderProducer(NuSliceReBuilderProducer const&) = delete;
  NuSliceReBuilderProducer(NuSliceReBuilderProducer&&) = delete;
  NuSliceReBuilderProducer& operator=(NuSliceReBuilderProducer const&) = delete;
  NuSliceReBuilderProducer& operator=(NuSliceReBuilderProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  // Declare member data here.
  // std::string fPfpLabel;
  // std::string fSliceLabel;
  // std::string fHitLabel;
  // std::string fHitTruthLabel;
};

NuSliceReBuilderProducer::NuSliceReBuilderProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  // , fPfpLabel(p.get<std::string>("PfpLabel", "pandora"))
  // , fSliceLabel(p.get<std::string>("SliceLabel", "pandora"))
  // , fHitLabel(p.get<std::string>("HitLabel", "gaushit"))
  // , fHitTruthLabel(p.get<std::string>("HitTruthLabel", "gaushitTruthMatch"))
// More initializers here.A
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::PFParticle>>();
  produces<std::vector<recob::Slice>>();
  produces<std::vector<recob::PCAxis>>();
  produces<std::vector<recob::Cluster>>();
  produces<std::vector<recob::Vertex>>();
  produces<std::vector<recob::Shower>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<std::vector<recob::Track>>();
  produces<std::vector<larpandoraobj::PFParticleMetadata>>();
  //
  produces<art::Assns<recob::Slice,recob::Hit,void>>();
  produces<art::Assns<recob::PFParticle,recob::PCAxis,void>>();
  produces<art::Assns<recob::PFParticle,recob::Vertex,void>>();
  produces<art::Assns<recob::PFParticle,recob::Slice,void>>();
  produces<art::Assns<recob::PFParticle,recob::SpacePoint,void>>();
  produces<art::Assns<recob::Shower,recob::Hit,void>>();
  produces<art::Assns<recob::PFParticle,recob::Track,void>>();
  produces<art::Assns<recob::Shower,recob::PCAxis,void>>();
  produces<art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void>>();
  produces<art::Assns<recob::Cluster,recob::Hit,void>>();
  produces<art::Assns<recob::Track,recob::Hit,recob::TrackHitMeta>>();
  produces<art::Assns<recob::SpacePoint,recob::Hit,void>>();
  produces<art::Assns<recob::PFParticle,recob::Cluster,void>>();
  produces<art::Assns<recob::PFParticle,recob::Shower,void>>();
}

void NuSliceReBuilderProducer::produce(art::Event& e)
{
  //
  auto outputPFP        = std::make_unique<std::vector<recob::PFParticle>>();
  auto outputSlice      = std::make_unique<std::vector<recob::Slice>>();
  auto outputPCAxis     = std::make_unique<std::vector<recob::PCAxis>>();
  auto outputCluster    = std::make_unique<std::vector<recob::Cluster>>();
  auto outputVertex     = std::make_unique<std::vector<recob::Vertex>>();
  auto outputShower     = std::make_unique<std::vector<recob::Shower>>();
  auto outputSpacePoint = std::make_unique<std::vector<recob::SpacePoint>>();
  auto outputTrack      = std::make_unique<std::vector<recob::Track>>();
  auto outputPFMeta     = std::make_unique<std::vector<larpandoraobj::PFParticleMetadata>>();
  //
  auto outPFPPCAAssns = std::make_unique<art::Assns<recob::PFParticle,recob::PCAxis,void>>();
  auto outPFPVtxAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Vertex,void>>();
  auto outPFPSlcAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Slice,void>>();
  auto outPFPSpsAssns = std::make_unique<art::Assns<recob::PFParticle,recob::SpacePoint,void>>();
  auto outPFPTrkAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Track,void>>();
  auto outPFPMetAssns = std::make_unique<art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void>>();
  auto outPFPCluAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Cluster,void>>();
  auto outPFPShrAssns = std::make_unique<art::Assns<recob::PFParticle,recob::Shower,void>>();
  auto outSlcHitAssns = std::make_unique<art::Assns<recob::Slice,recob::Hit,void>>();
  auto outCluHitAssns = std::make_unique<art::Assns<recob::Cluster,recob::Hit,void>>();
  auto outSpsHitAssns = std::make_unique<art::Assns<recob::SpacePoint,recob::Hit,void>>();
  auto outTrkHitAssns = std::make_unique<art::Assns<recob::Track,recob::Hit,recob::TrackHitMeta>>();
  auto outShrHitAssns = std::make_unique<art::Assns<recob::Shower,recob::Hit,void>>();
  auto outShrPCAAssns = std::make_unique<art::Assns<recob::Shower,recob::PCAxis,void>>();

  art::PtrMaker<recob::Track> trkPtrMaker(e);
  art::PtrMaker<recob::Shower> shrPtrMaker(e);
  art::PtrMaker<recob::Slice> slcPtrMaker(e);
  art::PtrMaker<recob::Cluster> cluPtrMaker(e);
  art::PtrMaker<recob::PFParticle> pfpPtrMaker(e);
  art::PtrMaker<recob::Vertex> vtxPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spsPtrMaker(e);
  art::PtrMaker<recob::PCAxis> pcaPtrMaker(e);
  art::PtrMaker<larpandoraobj::PFParticleMetadata> pfmetaPtrMaker(e);

  // original pandora slice products
  art::ValidHandle<std::vector<recob::PFParticle>> inputPndrPFParticle = e.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  art::ValidHandle<std::vector<recob::Slice>> inputPndrSlice = e.getValidHandle<std::vector<recob::Slice>>("pandora");
  //art::ValidHandle<std::vector<recob::PCAxis>> inputPndrPCAxis = e.getValidHandle<std::vector<recob::PCAxis>>("pandora");
  art::ValidHandle<std::vector<recob::Cluster>> inputPndrCluster = e.getValidHandle<std::vector<recob::Cluster>>("pandora");
  //art::ValidHandle<std::vector<recob::Vertex>> inputPndrVertex = e.getValidHandle<std::vector<recob::Vertex>>("pandora");
  art::ValidHandle<std::vector<recob::SpacePoint>> inputPndrSpacePoint = e.getValidHandle<std::vector<recob::SpacePoint>>("pandora");
  art::ValidHandle<std::vector<recob::Shower>> inputPndrShower = e.getValidHandle<std::vector<recob::Shower>>("pandora");
  art::ValidHandle<std::vector<recob::Track>> inputPndrTrack = e.getValidHandle<std::vector<recob::Track>>("pandora");
  //art::ValidHandle<std::vector<larpandoraobj::PFParticleMetadata>> inputPndrPFParticleMetadata = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>("pandora");

  auto inputPndrPFPPCAAssns = std::unique_ptr<art::FindOneP<recob::PCAxis>>    (new art::FindOneP<recob::PCAxis>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPVtxAssns = std::unique_ptr<art::FindOneP<recob::Vertex>>    (new art::FindOneP<recob::Vertex>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPSlcAssns = std::unique_ptr<art::FindOneP<recob::Slice>>     (new art::FindOneP<recob::Slice>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPTrkAssns = std::unique_ptr<art::FindOneP<recob::Track>>     (new art::FindOneP<recob::Track>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPShrAssns = std::unique_ptr<art::FindOneP<recob::Shower>>    (new art::FindOneP<recob::Shower>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPCluAssns = std::unique_ptr<art::FindManyP<recob::Cluster>>   (new art::FindManyP<recob::Cluster>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPSpsAssns = std::unique_ptr<art::FindManyP<recob::SpacePoint>>(new art::FindManyP<recob::SpacePoint>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrPFPMetAssns = std::unique_ptr<art::FindManyP<larpandoraobj::PFParticleMetadata>>(new art::FindManyP<larpandoraobj::PFParticleMetadata>(inputPndrPFParticle, e, "pandora"));
  auto inputPndrSlcHitAssns = std::unique_ptr<art::FindManyP<recob::Hit>>       (new art::FindManyP<recob::Hit>(inputPndrSlice, e, "pandora"));
  auto inputPndrCluHitAssns = std::unique_ptr<art::FindManyP<recob::Hit>>       (new art::FindManyP<recob::Hit>(inputPndrCluster, e, "pandora"));
  auto inputPndrSpsHitAssns = std::unique_ptr<art::FindManyP<recob::Hit>>       (new art::FindManyP<recob::Hit>(inputPndrSpacePoint, e, "pandora"));
  auto inputPndrShrHitAssns = std::unique_ptr<art::FindManyP<recob::Hit>>       (new art::FindManyP<recob::Hit>(inputPndrShower, e, "pandora"));
  auto inputPndrTrkHitAssns = std::unique_ptr<art::FindManyP<recob::Hit,recob::TrackHitMeta>>(new art::FindManyP<recob::Hit,recob::TrackHitMeta>(inputPndrTrack, e, "pandora"));

  // keep track of what pfps were removed
  art::ValidHandle<std::vector<recob::PFParticle>> removedPFParticle = e.getValidHandle<std::vector<recob::PFParticle>>("nugraphshowerhits");

  // shower pfps from NG2 that needs to be added
  art::ValidHandle<std::vector<recob::Shower>> inputShowers = e.getValidHandle<std::vector<recob::Shower>>("showerreco3d");
  art::ValidHandle<std::vector<recob::SpacePoint>> inputShowerSpacePoints = e.getValidHandle<std::vector<recob::SpacePoint>>("showerreco3d");
  art::ValidHandle<std::vector<recob::Cluster>> inputShowerClusters = e.getValidHandle<std::vector<recob::Cluster>>("cmerger");
  //
  auto assocShrPfp = std::unique_ptr<art::FindOneP<recob::PFParticle>>(new art::FindOneP<recob::PFParticle>(inputShowers, e, "showerreco3d"));
  auto assocShrClu = std::unique_ptr<art::FindManyP<recob::Cluster>>(new art::FindManyP<recob::Cluster>(inputShowers, e, "showerreco3d"));
  auto assocShrSps = std::unique_ptr<art::FindManyP<recob::SpacePoint>>(new art::FindManyP<recob::SpacePoint>(inputShowers, e, "showerreco3d"));
  auto assocShrPca = std::unique_ptr<art::FindOneP<recob::PCAxis>>(new art::FindOneP<recob::PCAxis>(inputShowers, e, "showerreco3d"));
  auto assocShrHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputShowers, e, "showerreco3d"));
  auto assocShrVtx = std::unique_ptr<art::FindOneP<recob::Vertex>>(new art::FindOneP<recob::Vertex>(inputShowers, e, "showerreco3d"));
  //
  auto assocSpsHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputShowerSpacePoints, e, "showerreco3d"));
  auto assocCluHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputShowerClusters, e, "cmerger"));

  std::cout << "input sizes=" << inputPndrPFParticle->size() << " " << inputPndrPFPTrkAssns->size() << std::endl;

  struct pfphierarchy {
    int origid;
    int newid;
    int parentid;
  };
  std::vector<pfphierarchy> pfphvec;

  size_t d_pfp_idx = 1;//reserve 0 for the nu pfp
  for (unsigned int ip=0; ip<inputPndrPFParticle->size(); ip++) {
    art::Ptr<recob::PFParticle> pfp(inputPndrPFParticle, ip);
    bool removed = false;
    for (unsigned int ir=0; ir<removedPFParticle->size(); ir++) {
      art::Ptr<recob::PFParticle> rpfp(removedPFParticle, ir);
      if (pfp->Self()==rpfp->Self()) {
	removed = true;
	break;
      }      
    }
    /////////
    {
      int nhits = 0;
      auto clusters = inputPndrPFPCluAssns->at(pfp.key());
      for (auto clu : clusters) {
	auto hits = inputPndrCluHitAssns->at(clu.key());
	nhits += hits.size();
      }
      std::cout << "pfp->PdgCode()=" << pfp->PdgCode() << " primary=" << pfp->IsPrimary() << " key=" << pfp.key() << " self=" << pfp->Self() << " parent=" << pfp->Parent() << " ndaughters=" << pfp->NumDaughters() << " nhits=" << nhits << " removed=" << removed << std::endl;
    }
    /////////
    if (removed) continue;
    //
    if (pfp->IsPrimary()) {
      //if this is the neutrino we need a different strategy
      recob::PFParticle pfp_copy = recob::PFParticle(pfp->PdgCode(),0,std::numeric_limits<size_t>::max(),std::vector< size_t >());
      outputPFP->push_back(pfp_copy);
      pfphvec.push_back( {int(pfp->Self()), 0, -1} );
      //
      art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
      auto pfmetas = inputPndrPFPMetAssns->at(pfp.key());
      for (auto pfm : pfmetas) {
	// for (auto pfmm : pfm->GetPropertiesMap()) std::cout << " nu slice " << pfmm.first << " " << pfmm.second << std::endl;
	outputPFMeta->push_back(*pfm);
	auto pfm_ptr = pfmetaPtrMaker(outputPFMeta->size()-1);
	outPFPMetAssns->addSingle(pfp_ptr,pfm_ptr);
      }
      //
    } else {
      //skip pfparticles with no space points, they will have invalid track pointers
      if (inputPndrPFPSpsAssns->at(pfp.key()).size()==0 || inputPndrPFPTrkAssns->at(pfp.key()).isAvailable()==0) continue;
      //NG2 told us this is not a shower, so let's make it a track
      recob::PFParticle pfp_copy(13,d_pfp_idx,0,std::vector< size_t >());
      outputPFP->push_back(pfp_copy);
      pfphvec.push_back( {int(pfp->Self()), int(d_pfp_idx), int(pfp->Parent())} );
      d_pfp_idx++;
    }
    art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
    if (pfp->PdgCode()==13) {
      //track
      //std::cout << "key=" << pfp.key() << " size=" << inputPndrPFPTrkAssns->size() << std::endl;
      //auto tp = inputPndrPFPTrkAssns->at(pfp.key());
      //std::cout << "tp valid=" << tp.isAvailable() << " key=" << tp.key() << " sps=" << inputPndrPFPSpsAssns->at(pfp.key()).size() << " cls=" << inputPndrPFPCluAssns->at(pfp.key()).size() << std::endl;
      recob::Track track = *inputPndrPFPTrkAssns->at(pfp.key());
      outputTrack->push_back(track);
      auto trk_ptr = trkPtrMaker(outputTrack->size()-1);
      outPFPTrkAssns->addSingle(pfp_ptr,trk_ptr);
      auto hits = inputPndrTrkHitAssns->at(inputPndrPFPTrkAssns->at(pfp.key()).key());
      auto meta = inputPndrTrkHitAssns->data(inputPndrPFPTrkAssns->at(pfp.key()).key());
      for (size_t ih=0; ih<hits.size(); ih++) {
	outTrkHitAssns->addSingle(trk_ptr,hits[ih],*meta[ih]);
      }
      //
      auto pfmetas = inputPndrPFPMetAssns->at(pfp.key());
      for (auto pfm : pfmetas) {
	//for (auto pfmm : pfm->GetPropertiesMap()) std::cout << pfmm.first << " " << pfmm.second << std::endl;
	outputPFMeta->push_back(*pfm);
	auto pfm_ptr = pfmetaPtrMaker(outputPFMeta->size()-1);
	outPFPMetAssns->addSingle(pfp_ptr,pfm_ptr);
      }
    } else if (pfp->PdgCode()==11) {
      /*
      //shower -- DO WE WANT TO DO THIS??? NG2 told us this is not a shower
      recob::Shower shower = *inputPndrPFPShrAssns->at(pfp.key());
      outputShower->push_back(shower);
      auto shr_ptr = shrPtrMaker(outputShower->size()-1);
      outPFPShrAssns->addSingle(pfp_ptr,shr_ptr);
      //
      auto pca = inputPndrPFPPCAAssns->at(pfp.key());
      outputPCAxis->push_back(*pca);
      auto pca_ptr = pcaPtrMaker(outputPCAxis->size()-1);
      outPFPPCAAssns->addSingle(pfp_ptr,pca_ptr);
      outShrPCAAssns->addSingle(shr_ptr,pca_ptr);
      auto hits = inputPndrShrHitAssns->at(inputPndrPFPShrAssns->at(pfp.key()).key());
      for (auto h : hits) {
	outShrHitAssns->addSingle(shr_ptr,h);	
      }
      */
      // force it to have track-like score
      std::map<std::string,float> mmap;
      mmap.emplace("TrackScore",1.);
      larpandoraobj::PFParticleMetadata pfm(mmap);
      outputPFMeta->push_back(pfm);
      auto pfm_ptr = pfmetaPtrMaker(outputPFMeta->size()-1);
      outPFPMetAssns->addSingle(pfp_ptr,pfm_ptr);
    }
    //
    auto clusters = inputPndrPFPCluAssns->at(pfp.key());
    for (auto clu : clusters) {
      outputCluster->push_back(*clu);
      auto clu_ptr = cluPtrMaker(outputCluster->size()-1);
      outPFPCluAssns->addSingle(pfp_ptr,clu_ptr);
      auto hits = inputPndrCluHitAssns->at(clu.key());
      for (auto h : hits) {
	outCluHitAssns->addSingle(clu_ptr,h);	
      }
    }
    //
    auto spacepoints = inputPndrPFPSpsAssns->at(pfp.key());
    for (auto sps : spacepoints) {
      outputSpacePoint->push_back(*sps);
      auto sps_ptr = spsPtrMaker(outputSpacePoint->size()-1);
      outPFPSpsAssns->addSingle(pfp_ptr,sps_ptr);
      auto hits = inputPndrSpsHitAssns->at(sps.key());
      for (auto h : hits) {
	outSpsHitAssns->addSingle(sps_ptr,h);	
      }
    }
    //
    auto vtx = inputPndrPFPVtxAssns->at(pfp.key());
    outputVertex->push_back(*vtx);
    auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
    outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
    //
    auto slice = inputPndrPFPSlcAssns->at(pfp.key());
    if (outputSlice->size()==0) {//it should be guaranteed that we don't have >1 slice at this point
      outputSlice->push_back(*slice);
      auto slc_ptr = slcPtrMaker(outputSlice->size()-1);
      auto hits = inputPndrSlcHitAssns->at(slice.key());
      for (auto h : hits) {
	outSlcHitAssns->addSingle(slc_ptr,h);	
      }
    }
    auto slc_ptr = slcPtrMaker(outputSlice->size()-1);
    outPFPSlcAssns->addSingle(pfp_ptr,slc_ptr);
  }

  //
  for (unsigned int is=0; is<inputShowers->size(); is++) {
    art::Ptr<recob::Shower> shr(inputShowers, is);
    auto pfp = assocShrPfp->at(shr.key());
    // we need to redefine the pfp id, etc.
    {
      std::cout << "pfp->PdgCode()=" << pfp->PdgCode() << " self=" << pfp->Self() << " parent=" << pfp->Parent() << " ndaughters=" << pfp->NumDaughters() << " nhits=" << assocShrHit->at(shr.key()).size() << " added" << std::endl;
    }
    recob::PFParticle pfp_copy(pfp->PdgCode(),d_pfp_idx,0,std::vector< size_t >());
    outputPFP->push_back(pfp_copy);
    pfphvec.push_back( {-1, int(d_pfp_idx), -2} );
    d_pfp_idx++;
    art::Ptr<recob::PFParticle> pfp_ptr = pfpPtrMaker(outputPFP->size()-1);
    //
    // larpandoraobj::PFParticleMetadata
    std::map<std::string,float> mmap;
    mmap.emplace("TrackScore",0);
    larpandoraobj::PFParticleMetadata pfm(mmap);
    outputPFMeta->push_back(pfm);
    auto pfm_ptr = pfmetaPtrMaker(outputPFMeta->size()-1);
    outPFPMetAssns->addSingle(pfp_ptr,pfm_ptr);
    //shower
    outputShower->push_back(*shr);
    auto shr_ptr = shrPtrMaker(outputShower->size()-1);
    outPFPShrAssns->addSingle(pfp_ptr,shr_ptr);
    //
    auto pca = assocShrPca->at(shr.key());
    outputPCAxis->push_back(*pca);
    auto pca_ptr = pcaPtrMaker(outputPCAxis->size()-1);
    outPFPPCAAssns->addSingle(pfp_ptr,pca_ptr);
    outShrPCAAssns->addSingle(shr_ptr,pca_ptr);
    auto hits = assocShrHit->at(shr.key());
    for (auto h : hits) {
      outShrHitAssns->addSingle(shr_ptr,h);	
    }
    //
    auto clusters = assocShrClu->at(shr.key());
    for (auto clu : clusters) {
      outputCluster->push_back(*clu);
      auto clu_ptr = cluPtrMaker(outputCluster->size()-1);
      outPFPCluAssns->addSingle(pfp_ptr,clu_ptr);
      auto hits = assocCluHit->at(clu.key());
      for (auto h : hits) {
	outCluHitAssns->addSingle(clu_ptr,h);	
      }
    }
    //
    auto spacepoints = assocShrSps->at(shr.key());
    for (auto sps : spacepoints) {
      outputSpacePoint->push_back(*sps);
      auto sps_ptr = spsPtrMaker(outputSpacePoint->size()-1);
      outPFPSpsAssns->addSingle(pfp_ptr,sps_ptr);
      auto hits = assocSpsHit->at(sps.key());
      for (auto h : hits) {
	outSpsHitAssns->addSingle(sps_ptr,h);	
      }
    }
    //
    auto vtx = assocShrVtx->at(shr.key());
    outputVertex->push_back(*vtx);
    auto vtx_ptr = vtxPtrMaker(outputVertex->size()-1);
    outPFPVtxAssns->addSingle(pfp_ptr,vtx_ptr);
    //
    auto slc_ptr = slcPtrMaker(outputSlice->size()-1);
    outPFPSlcAssns->addSingle(pfp_ptr,slc_ptr);
  }

  //restore the PFP hierarchy, to the extent that we can
  for (size_t idx=0;idx<outputPFP->size();idx++) {
    pfphierarchy h = pfphvec.at(idx);
    std::vector<size_t> daughters;
    size_t parent = 0;
    for (size_t idx2=0;idx2<outputPFP->size();idx2++) {
      if (idx2==idx) continue;
      pfphierarchy h2 = pfphvec.at(idx2);
      if (h.origid>=0 && h2.parentid==h.origid) daughters.push_back(h2.newid); // add daughters from old hierarchy
      if (h.newid==0 && h2.origid<0) daughters.push_back(h2.newid); // add newly created pfps as daughters of the nu pfp
      if (h2.origid>=0 && h.parentid==h2.origid) parent = h2.newid; // set the parent from the old hierarchy, otherwise use the nu pfp as default
    }
    if (h.origid>=0 && h.parentid==-1) parent = std::numeric_limits<size_t>::max();//restore the parent for the neutrino
    outputPFP->at(idx) = recob::PFParticle(outputPFP->at(idx).PdgCode(),h.newid,parent,daughters);
    std::cout << "output pfp idx=" << idx << " pdg=" << outputPFP->at(idx).PdgCode() << " newid=" << h.newid << " parent=" << parent << " ndaughters=" << daughters.size();
    std::cout << " daughters=";
    if (daughters.size()) {
      for (auto d : daughters) std::cout << d << " ";
      std::cout << std::endl;
    }
    else std::cout << "none" << std::endl;
  }

  //review pfps
  for (auto pfp : *outputPFP) std::cout << "final pfp->PdgCode()=" << pfp.PdgCode() << " self=" << pfp.Self() << " parent=" << pfp.Parent() << " ndaughters=" << pfp.NumDaughters() << std::endl;

  e.put(std::move(outputPFP));
  e.put(std::move(outputSlice));  
  e.put(std::move(outputPCAxis));
  e.put(std::move(outputCluster));
  e.put(std::move(outputVertex));
  e.put(std::move(outputShower));
  e.put(std::move(outputSpacePoint));
  e.put(std::move(outputTrack));
  e.put(std::move(outputPFMeta));
  e.put(std::move(outSlcHitAssns));
  e.put(std::move(outPFPPCAAssns));
  e.put(std::move(outPFPVtxAssns));
  e.put(std::move(outPFPSlcAssns));
  e.put(std::move(outPFPSpsAssns));
  e.put(std::move(outShrHitAssns));
  e.put(std::move(outPFPTrkAssns));
  e.put(std::move(outShrPCAAssns));
  e.put(std::move(outPFPMetAssns));
  e.put(std::move(outCluHitAssns));
  e.put(std::move(outTrkHitAssns));
  e.put(std::move(outSpsHitAssns));
  e.put(std::move(outPFPCluAssns));
  e.put(std::move(outPFPShrAssns));
}

DEFINE_ART_MODULE(NuSliceReBuilderProducer)
