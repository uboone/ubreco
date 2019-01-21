////////////////////////////////////////////////////////////////////////
// Class:       PandoraEventAnalyser
// Plugin Type: analyzer (art v2_11_03)
// File:        PandoraEventAnalyser_module.cc
//
// Generated at Sun Oct 28 14:38:04 2018 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Utilities/make_tool.h"
#include "FlashMatchingToolBase_tool.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "TTree.h"

// art TOOLS
#include "art/Utilities/ToolMacros.h"

class PandoraEventAnalyser;


class PandoraEventAnalyser : public art::EDAnalyzer {
public:
  explicit PandoraEventAnalyser(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraEventAnalyser(PandoraEventAnalyser const &) = delete;
  PandoraEventAnalyser(PandoraEventAnalyser &&) = delete;
  PandoraEventAnalyser & operator = (PandoraEventAnalyser const &) = delete;
  PandoraEventAnalyser & operator = (PandoraEventAnalyser &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:
  std::string fPFPproducer, fSpacePointproducer;

  std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

  TTree* _tree;
  float _score;

  // PFP map
  std::map<unsigned int, unsigned int> _pfpmap;

};

DEFINE_ART_MODULE(PandoraEventAnalyser)

PandoraEventAnalyser::PandoraEventAnalyser(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  fPFPproducer  = p.get<std::string>("PFPproducer");
  fSpacePointproducer  = p.get<std::string>("SpacePointproducer");

  const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
  
  std::cout << "Configuring flsh config" << std::endl;
  _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);
  //_flashmatchTool->configure(flashconfig);
  std::cout << "Done" << std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("flashmatch","flashmatching tree");
  _tree->Branch("_score",&_score,"score/F");

}

void PandoraEventAnalyser::analyze(art::Event const & e)
{
  std::cout << "Test: Hello World" << std::endl;
  // Implementation of required member function here.

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  
  // grab clusters associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
  
  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPFPproducer);    
  
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);
  
  /*
  // ADDITION FROM PETRILLO
  e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointproducer);
  
  // grab the hits associated to the PFParticles
  auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::SpacePoint>::find(pfp_h, e, fPFPproducer);
  */
  
  std::cout << "There are " << pfp_h->size() << " pfparticles in the event " << std::endl;
  
  // fill map: pfparticle Self() -> index/key
  _pfpmap.clear();
  for (unsigned int p=0; p < pfp_h->size(); p++)
    _pfpmap[pfp_h->at(p).Self()] = p;

  for (unsigned int p=0; p < pfp_h->size(); p++){
    
    auto const& pfp = pfp_h->at(p);
    
    /*
    // get metadata for this PFP
    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));
    
    if (!pfParticleMetadataList.empty()) {
    
    for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
    {
    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
    //const larpandoraobj::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
    if (!pfParticlePropertiesMap.empty())
    //std::cout << " Found PFParticle " << pfp.Self() << " with: " << std::endl;
    for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
    std::cout << it->first << " ... " << it->second << std::endl;
    }
    }
    }
    */
    
    // start from primary PFParticles
    if (pfp.IsPrimary() == false) continue;
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    
    // now build vectors of PFParticles, space-points, and hits for this slice
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint>>> spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit>>> hit_v_v;
    
    std::cout << "creating PFP hierarchy." << std::endl;
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);
    std::cout << "There are " << pfp_ptr_v.size() << " PFParticles in this hierarchy " << std::endl << std::endl;
    
    // go through these pfparticles and fill info needed for matching
    for (size_t i=0; i < pfp_ptr_v.size(); i++) {
      
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      
      //std::cout << "key is " << key << std::endl;
      
      pfp_v.push_back(pfp);
      
      
      //auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      
      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	auto const& spkey = spacepoint_ptr_v.at(sp).key();
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	  hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	}// for all hits associated to this spacepoint
      }// fpr all spacepoints
      
      
      //std::cout << "there are " << hit_ptr_v.size() << " hits" << std::endl;
      //std::cout << "there are " << spacepoint_ptr_v.size() << " spacepoints" << std::endl;
      
      spacepoint_v_v.push_back( spacepoint_ptr_v );
      hit_v_v.push_back( hit_ptr_v );
      
    }// for all pfp pointers
    
    // ready to call flash-matching
    //float score = _flashmatchTool->ClassifySlice();
    _score = _flashmatchTool->ClassifySlice(e, pfp_ptr_v, spacepoint_v_v, hit_v_v);
    //float score = _flashmatchTool->ClassifySlice(e);//, pfp_ptr_v, spacepoint_v_v, hit_v_v);
    
    std::cout << "score is " << _score << std::endl;

    _tree->Fill();
    
    
  }// for all PFParticles
  
  return;
}

void PandoraEventAnalyser::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
  auto daughters = pfp_ptr->Daughters();
  
  pfp_v.push_back(pfp_ptr);
  
  std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
  for(auto const& daughterid : daughters) {

    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
    AddDaughters(pfp_ptr, pfp_h, pfp_v);
    
  }// for all daughters
  
  return;
}
