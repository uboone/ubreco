////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeFlashMatchChi2
// Plugin Type: analyzer (art v3_01_02)
// File:        AnalyzeFlashMatchChi2_module.cc
//
// Generated at Tue Jun  4 13:36:23 2019 by David Caratelli using cetskelgen
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

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

class AnalyzeFlashMatchChi2;


class AnalyzeFlashMatchChi2 : public art::EDAnalyzer {
public:
  explicit AnalyzeFlashMatchChi2(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeFlashMatchChi2(AnalyzeFlashMatchChi2 const&) = delete;
  AnalyzeFlashMatchChi2(AnalyzeFlashMatchChi2&&) = delete;
  AnalyzeFlashMatchChi2& operator=(AnalyzeFlashMatchChi2 const&) = delete;
  AnalyzeFlashMatchChi2& operator=(AnalyzeFlashMatchChi2&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  
  TTree* _tree;
  int _evt, _run, _sub, _pdg;
  float _chi2, _time, _nuscore, _clearcosmic;

  art::InputTag fPFPproducer, fT0producer;
  

};


AnalyzeFlashMatchChi2::AnalyzeFlashMatchChi2(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fT0producer  = p.get<art::InputTag>("T0producer" );

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("slicet0tree","slice t0 tree");
  _tree->Branch("evt",&_evt,"evt/I");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("sub",&_sub,"sub/I");
  _tree->Branch("chi2",&_chi2,"chi2/F");
  _tree->Branch("pdg",&_pdg,"pdg/I");
  _tree->Branch("time",&_time,"time/F");
  _tree->Branch("nuscore",&_nuscore,"nuscore/F");
  _tree->Branch("clearcosmic",&_clearcosmic,"clearcosmic/F");

}

void AnalyzeFlashMatchChi2::analyze(art::Event const& e)
{

  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  
  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfp_meta_assn_v(pfp_h, e, fPFPproducer);

  // grab associated T0
  art::FindManyP< anab::T0 > pfp_t0_assn_v( pfp_h, e, fT0producer);
  
  // loop through PFParticles
  for (size_t p=0; p < pfp_h->size(); p++) {

    _clearcosmic = 0;
    _nuscore     = 0;
    _pdg         = 0;
    _time        = 0;
    _chi2        = 0;
    
    const recob::PFParticle pfp = pfp_h->at(p);

    // skip non-primary PFP
    if (pfp.IsPrimary() == false) continue;

    _pdg = pfp.PdgCode();

    // get t0 for this PFP
    const std::vector< art::Ptr<anab::T0> > &pfp_t0_v(pfp_t0_assn_v.at(p));

    if (pfp_t0_v.size() == 1) {
      _chi2 = pfp_t0_v[0]->TriggerConfidence();
      std::cout << "Chi2 is " << _chi2 << std::endl;
      _time = pfp_t0_v[0]->Time();
    }
    
    // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfp_meta_assn_v.at(p));
    
    if (!pfParticleMetadataList.empty()) {
      
      for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
	{
	  const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	    
	    if (!pfParticlePropertiesMap.empty())
	      
	      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {

		// if clear cosmic
		if ( (it->first == "IsClearCosmic") && (it->second == 1) ) 
		  _clearcosmic = 1;
		
		// if there is a nu score
		if (it->first == "NuScore") 
		  _nuscore = it->second;

	      }// if neutrino score too low
	}// for metadata map
    }// if metadata exits

    _tree->Fill();

  }// for all PFP
  
  return;
}

void AnalyzeFlashMatchChi2::beginJob()
{
  // Implementation of optional member function here.
}

void AnalyzeFlashMatchChi2::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(AnalyzeFlashMatchChi2)
