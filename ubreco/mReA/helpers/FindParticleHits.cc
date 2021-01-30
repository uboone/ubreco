//////////////////////////////////////////////////////////
// Implementation file for FindParticleHits class
// select hits associated with input PF or MC particle
// author: M. Rosenberg (2/12/2020)
//////////////////////////////////////////////////////////


#include "ubreco/mReAReco2/helpers/FindParticleHits.h"

#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
//#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/Cluster.h"


namespace mReA {


FindParticleHits::FindParticleHits(){
  m_pfp_producer = "pandora";
  m_cluster_producer = "pandora";
  m_hit_producer = "gaushit";
  m_pfpToCls_producer = "pandora";
  m_clsToHit_producer = "pandora";
  m_mcp_producer = "largeant";
  m_backTrack_producer = "gaushitTruthMatch";
  m_useDaughters = true;
  PrintConfig();
}

FindParticleHits::FindParticleHits(std::string pfp, std::string cluster, std::string hit,
 std::string pfpToCls, std::string clsToHit, std::string mcp, std::string backTrack,
 std::string hitProc, std::string mcpProc, std::string backTrackProc, bool useDs){
  m_pfp_producer = pfp;
  m_cluster_producer = cluster;
  m_hit_producer = hit;
  m_pfpToCls_producer = pfpToCls;
  m_clsToHit_producer = clsToHit;
  if(!mcp.empty()) m_mcp_producer = mcp;
  if(!backTrack.empty()) m_backTrack_producer = backTrack;
  if(!hitProc.empty()) m_hit_producer += "::"+hitProc;
  if(!mcpProc.empty()) m_mcp_producer += "::"+mcpProc;
  if(!backTrackProc.empty()) m_backTrack_producer += "::"+backTrackProc;
  m_useDaughters = useDs;
  PrintConfig();
}

void FindParticleHits::Initialize(std::string pfp, std::string cluster, std::string hit,
 std::string pfpToCls, std::string clsToHit, std::string mcp, std::string backTrack,
 std::string hitProc, std::string mcpProc, std::string backTrackProc, bool useDs){
  m_pfp_producer = pfp;
  m_cluster_producer = cluster;
  m_hit_producer = hit;
  m_pfpToCls_producer = pfpToCls;
  m_clsToHit_producer = clsToHit;
  m_mcp_producer = mcp;
  m_backTrack_producer = backTrack;
  if(!hitProc.empty()) m_hit_producer += "::"+hitProc;
  if(!mcpProc.empty()) m_mcp_producer += "::"+mcpProc;
  if(!backTrackProc.empty()) m_backTrack_producer += "::"+backTrackProc;
  m_useDaughters = useDs;
  PrintConfig();
}

void FindParticleHits::Initialize(std::string pndr, std::string hit, std::string mcp,
 std::string backTrk, bool useDs){
  m_pfp_producer = pndr;
  m_cluster_producer = pndr;
  m_hit_producer = hit;
  m_pfpToCls_producer = pndr;
  m_clsToHit_producer = pndr;
  m_mcp_producer = mcp;
  m_backTrack_producer = backTrk;
  m_useDaughters = useDs;
  PrintConfig();
}

void FindParticleHits::PrintConfig() const {
    std::cout << "FindParticleHits initialized/re-initialized with..." << std::endl;
    std::cout << "  m_pfp_producer = " << m_pfp_producer << std::endl;
    std::cout << "  m_cluster_producer = " << m_cluster_producer << std::endl;
    std::cout << "  m_hit_producer = " << m_hit_producer << std::endl;
    std::cout << "  m_pfpToCls_producer = " << m_pfpToCls_producer << std::endl;
    std::cout << "  m_clsToHit_producer = " << m_clsToHit_producer << std::endl;
    std::cout << "  m_mcp_producer = " << m_mcp_producer << std::endl;
    std::cout << "  m_backTrack_producer = " << m_backTrack_producer << std::endl;
    std::cout << "  m_useDaughters = " << m_useDaughters << std::endl;
}

void FindParticleHits::getPFPHits(const std::map<size_t, art::Ptr<recob::PFParticle> >& pfParticleIdMap,
 const art::Event& evt, const art::Ptr<recob::PFParticle>& inputPFPart,
 std::vector< art::Ptr<recob::Hit> >& inputPFPartHits){

  //Get inputPFPart daughters
  std::vector<size_t> inputParticleChain;
  inputParticleChain.push_back(inputPFPart->Self());
  if(m_useDaughters){
    const std::vector<size_t>& partDaughters = inputPFPart->Daughters();
    if(!partDaughters.empty()){
      inputParticleChain.insert(inputParticleChain.end(), partDaughters.begin(), partDaughters.end());
      std::vector<size_t> inputDaughters(partDaughters);
      std::vector<size_t> newDaughters;
      int failSafe = 0;
      while(true){
        getRecoDaughterChain(pfParticleIdMap, inputParticleChain, inputDaughters, newDaughters);
        if(newDaughters.empty()) break;
        inputDaughters.clear();
        inputDaughters.insert(inputDaughters.end(), newDaughters.begin(), newDaughters.end());
        newDaughters.clear();
        if(failSafe > 1000000) throw cet::exception("FindParticleHits") <<
          "got stuck on infinite loop when locating PFParticle daughters." << std::endl;
        ++failSafe;
      }
    }
  }

  //Get inputPFPart (and daughter particle) hits
  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  evt.getByLabel(m_pfp_producer, pfParticleHandle);
  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  evt.getByLabel(m_cluster_producer, clusterHandle);
  art::FindManyP<recob::Cluster> pfPartToClusterAssoc(pfParticleHandle, evt, m_pfpToCls_producer);
  art::FindManyP<recob::Hit> clusterToHitAssoc(clusterHandle, evt, m_clsToHit_producer);
  inputPFPartHits.clear();
  inputPFPartHits.reserve(10000);

  for(auto iP = inputParticleChain.begin(); iP != inputParticleChain.end(); ++iP){
    const auto associatedClusters = pfPartToClusterAssoc.at( pfParticleIdMap.at(*iP).key() );
    for(unsigned int i = 0; i < associatedClusters.size(); ++i){
      const auto associatedHits = clusterToHitAssoc.at( associatedClusters[i].key() );
      for(unsigned int j = 0; j < associatedHits.size(); ++j){
        inputPFPartHits.push_back(associatedHits[j]);
      }
    }
  }

}


void FindParticleHits::getPFPHits(const art::Event& evt, const art::Ptr<recob::PFParticle>& inputPFPart,
 std::vector< art::Ptr<recob::Hit> >& inputPFPartHits){

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  evt.getByLabel(m_pfp_producer, pfParticleHandle);

  //fill map between PFParticle->Self() and PFParticle
  std::map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;
  for(unsigned int iPF = 0; iPF < pfParticleHandle->size(); ++iPF){
    const art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, iPF);
    if(!pfParticleIdMap.emplace(pfParticle->Self(), pfParticle).second)
      throw cet::exception("FindParticleHits") << "Repeated PFParticles in input!" << std::endl;
  }

  getPFPHits(pfParticleIdMap, evt, inputPFPart, inputPFPartHits);

}


void FindParticleHits::getMCPHitMap(const std::map<int, art::Ptr<simb::MCParticle> >& mcParticleIdMap,
 const art::Event& evt, const art::Ptr<simb::MCParticle>& inputMCPart,
 std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> >& inputMCPartHitsMap){

  //Get inputMCPart daughters
  std::vector<int> inputMCParticleChain;
  if(m_useDaughters){
    inputMCParticleChain.reserve(2000);
    inputMCParticleChain.push_back(inputMCPart->TrackId());
    int nDaughters = inputMCPart->NumberDaughters();
    std::vector<int> inputMCDaughters;
    inputMCDaughters.reserve(nDaughters);
    for(int i = 0; i < nDaughters; ++i){
      int daughterId = inputMCPart->Daughter(i);
      if(mcParticleIdMap.count(daughterId)){
        inputMCParticleChain.push_back(daughterId);
        inputMCDaughters.push_back(daughterId);
      }
    }
    if(!inputMCDaughters.empty()){
      std::vector<int> newMCDaughters;
      newMCDaughters.reserve(1000);
      int failSafe = 0;
      while(true){
        getMCDaughterChain(mcParticleIdMap, inputMCParticleChain, inputMCDaughters, newMCDaughters);
        if(newMCDaughters.empty()) break;
        inputMCDaughters.clear();
        inputMCDaughters.reserve(newMCDaughters.size());
        inputMCDaughters.insert(inputMCDaughters.end(), newMCDaughters.begin(), newMCDaughters.end());
        newMCDaughters.clear();
        newMCDaughters.reserve(1000);
        if(failSafe > 1000000) throw cet::exception("FindParticleHits") <<
          "got stuck on infinite loop when locating MCParticle daughters." << std::endl;
        ++failSafe;
      }
    }
  }
  else{ inputMCParticleChain.push_back(inputMCPart->TrackId()); }

  //Get inputMCPart (and daughter particle) hits
  //art::InputTag hitTag(m_hit_producer, "", m_hit_process);
  //art::InputTag backTrackTag(m_backTrack_producer, "", m_backTrack_process);
  art::Handle< std::vector<recob::Hit> > hitHandle;
  evt.getByLabel(m_hit_producer, hitHandle);
  art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcPartsAssoc(hitHandle,evt,m_backTrack_producer);

  for(unsigned int iH = 0; iH < hitHandle->size(); ++iH){

    const art::Ptr<recob::Hit> thisHit(hitHandle, iH);
    std::vector< art::Ptr<simb::MCParticle> > matchedMCPart_vec;
    std::vector< const anab::BackTrackerHitMatchingData* > match_vec;
    std::vector< BackTrackerWrapper > backTracker_vec;
    mcPartsAssoc.get(thisHit.key(), matchedMCPart_vec, match_vec);
    bool inputMCPartContrib = false;

    for(unsigned int iP = 0; iP < matchedMCPart_vec.size(); ++iP){
      if(matchedMCPart_vec[iP].isNull()) continue;
      int matchedId = matchedMCPart_vec[iP]->TrackId();
      bool isInputMCPart = false;
      auto iT = std::find(inputMCParticleChain.begin(), inputMCParticleChain.end(), matchedId);
      if(iT != inputMCParticleChain.end()){
        isInputMCPart = true;
        inputMCPartContrib = true;
      }
      BackTrackerWrapper backTracker(matchedId, isInputMCPart, *match_vec[iP]);
      backTracker_vec.emplace_back(std::move(backTracker));
    }

    if(inputMCPartContrib) inputMCPartHitsMap[thisHit] = backTracker_vec;

  }

}


void FindParticleHits::getMCPHitMap(const art::Event& evt, const art::Ptr<simb::MCParticle>& inputMCPart,
 std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> >& inputMCPartHitsMap){

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
  //art::InputTag mcpTag(m_mcp_producer, "", m_mcp_process);
  evt.getByLabel(m_mcp_producer, mcParticleHandle);

  std::map<int, art::Ptr<simb::MCParticle> > mcParticleIdMap;
  for (unsigned int i = 0; i < mcParticleHandle->size(); ++i){
    const art::Ptr<simb::MCParticle> mcParticle(mcParticleHandle, i);
    if(!mcParticleIdMap.emplace(mcParticle->TrackId(), mcParticle).second)
      throw cet::exception("FindParticleHits") << "Repeated MCParticles in input!" << std::endl;
  }

  getMCPHitMap(mcParticleIdMap, evt, inputMCPart, inputMCPartHitsMap);

}


void FindParticleHits::getMCPHits(const std::map<int, art::Ptr<simb::MCParticle> >& mcParticleIdMap,
 const art::Event& evt, const art::Ptr<simb::MCParticle>& inputMCPart,
 std::vector< art::Ptr<recob::Hit> >& inputMCPartHits, bool bestMatchOnly){

  std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> > mcPartHitMap;
  getMCPHitMap(mcParticleIdMap, evt, inputMCPart, mcPartHitMap);
  FillMCHitVector(mcPartHitMap, inputMCPartHits, bestMatchOnly);

}


void FindParticleHits::getMCPHits(const art::Event& evt, const art::Ptr<simb::MCParticle>& inputMCPart,
 std::vector< art::Ptr<recob::Hit> >& inputMCPartHits, bool bestMatchOnly){

  std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> > mcPartHitMap;
  getMCPHitMap(evt, inputMCPart, mcPartHitMap);
  FillMCHitVector(mcPartHitMap, inputMCPartHits, bestMatchOnly);

}


void FindParticleHits::FillMCHitVector(const std::map<art::Ptr<recob::Hit>,std::vector<BackTrackerWrapper> >&
 mcPartHitMap, std::vector< art::Ptr<recob::Hit> >& inputMCPartHits, bool bestMatchOnly){

  inputMCPartHits.reserve(10000);
  for(auto iM = mcPartHitMap.begin(); iM != mcPartHitMap.end(); ++iM){
    bool inputPartHit = false;
    for(auto iBT = (iM->second).begin(); iBT != (iM->second).end(); ++iBT){
      if(iBT->fromInput){
        if(bestMatchOnly && iBT->bt_isMaxIDE){ inputPartHit = true; break; }
        else{ inputPartHit = true; break; }
      }
    }
    if(inputPartHit) inputMCPartHits.push_back(iM->first);
  }

}


void FindParticleHits::getRecoDaughterChain(const std::map<size_t, art::Ptr<recob::PFParticle> >&
 pfParticleIdMap, std::vector<size_t>& partChainVec, const std::vector<size_t>& inputDaughters,
 std::vector<size_t>& newDaughters){

  for(auto it = inputDaughters.begin(); it != inputDaughters.end(); ++it){
    const std::vector<size_t>& currentDaughters = pfParticleIdMap.at(*it)->Daughters();
    if(!currentDaughters.empty()){
      partChainVec.insert(partChainVec.end(), currentDaughters.begin(), currentDaughters.end());
      newDaughters.insert(newDaughters.end(), currentDaughters.begin(), currentDaughters.end());
    }
  }

}


void FindParticleHits::getMCDaughterChain(const std::map<int, art::Ptr<simb::MCParticle> >&
 mcParticleIdMap, std::vector<int>& partChainVec, const std::vector<int>& inputDaughters,
 std::vector<int>& newDaughters){

  for(auto it = inputDaughters.begin(); it != inputDaughters.end(); ++it){
    int nDaughters = mcParticleIdMap.at(*it)->NumberDaughters();
    for(int i = 0; i < nDaughters; ++i){
      int daughterTrackId = mcParticleIdMap.at(*it)->Daughter(i);
      if(mcParticleIdMap.count(daughterTrackId)){
        partChainVec.push_back(daughterTrackId);
        newDaughters.push_back(daughterTrackId);
      }
    }
  }

}


} //end namespace mReA
