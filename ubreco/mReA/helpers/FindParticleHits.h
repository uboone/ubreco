//////////////////////////////////////////////////////////
// Header for FindParticleHits class
// select hits associated with input PF or MC particle
// author: M. Rosenberg (2/12/2020)
//////////////////////////////////////////////////////////

#ifndef UBRECO_MREA_FINDPARTICLE_HITS_H
#define UBRECO_MREA_FINDPARTICLE_HITS_H

#include <string>
#include <vector>
#include <map>

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "nusimdata/SimulationBase/MCParticle.h"


namespace mReA {


struct BackTrackerWrapper {
  int mcPartTrackId; // MCParticle->TrackId() for this particle
  bool fromInput; // is this the getMCPHits input particle (or one if its daughters, or their daughters, etc.)?
  //BackTrackerHitMatchingData members
  float bt_ideFraction;  //fraction of energy in hit from this particle
  int bt_isMaxIDE;       //is this particle the max contributor to this hit?
  float bt_ideNFraction; // fraction of number of electrons on the wire in hit from this particle
  int bt_isMaxIDEN;      // is this particle the max contributor to this hit in terms of number of electrons?
  float bt_numElectrons; // Number of electrons collected at the readout wire
  float bt_energy;       // energy deposited by ionization by this track ID [MeV]
  BackTrackerWrapper(){
    mcPartTrackId = -1;
    fromInput = false;
    bt_ideFraction = -1.;
    bt_isMaxIDE = -1;
    bt_ideNFraction = -1.;
    bt_isMaxIDEN = -1;
    bt_numElectrons = -1.;
    bt_energy = -1.;
  }
  BackTrackerWrapper(int track, bool isInput, const anab::BackTrackerHitMatchingData& matchData){
    mcPartTrackId = track;
    fromInput = isInput;
    bt_ideFraction = matchData.ideFraction;
    bt_isMaxIDE = matchData.isMaxIDE;
    bt_ideNFraction = matchData.ideNFraction;
    bt_isMaxIDEN = matchData.isMaxIDEN;
    bt_numElectrons = matchData.numElectrons;
    bt_energy = matchData.energy;
  }
};


class FindParticleHits {

  public:

    //constructor, initialize functions set producer labels
    FindParticleHits();
    FindParticleHits(std::string pfp, std::string cluster, std::string hit,
     std::string pfpToCls, std::string clsToHit, std::string mcp, std::string backTrack,
     std::string hitProc="", std::string mcpProc="", std::string backTrackProc="", bool useDs=true);
    void Initialize(std::string pfp, std::string cluster, std::string hit,
     std::string pfpToCls, std::string clsToHit, std::string mcp, std::string backTrack,
     std::string hitProc="", std::string mcpProc="", std::string backTrackProc="", bool useDs=true);
    void Initialize(std::string pndr, std::string hit, std::string mcp, std::string backTrk, bool useDs=true);

    //Fill inputPFPartHits vec with hits associated with inputPFPart and daughter particles
     //from map of PFParticle->Self() to PFParticle
    void getPFPHits(const std::map<size_t, art::Ptr<recob::PFParticle> >& pfParticleIdMap,
                    const art::Event& evt,
                    const art::Ptr<recob::PFParticle>& inputPFPart,
                    std::vector< art::Ptr<recob::Hit> >& inputPFPartHits);
     //from art event
    void getPFPHits(const art::Event& evt,
                    const art::Ptr<recob::PFParticle>& inputPFPart,
                    std::vector< art::Ptr<recob::Hit> >& inputPFPartHits);

    //Fill inputMCPartHitsMap map(hit -> vec<BackTrackerWrapper>) with all hits that have a
    //contribution from inputMCPart and daughter particles (if requested)
     //from map of MCParticle->TrackID() to MCParticle
    void getMCPHitMap(const std::map<int, art::Ptr<simb::MCParticle> >& mcParticleIdMap,
                    const art::Event& evt,
                    const art::Ptr<simb::MCParticle>& inputMCPart,
                    std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> >& inputMCPartHitsMap);
     //from art event
    void getMCPHitMap(const art::Event& evt,
                    const art::Ptr<simb::MCParticle>& inputMCPart,
                    std::map< art::Ptr<recob::Hit>, std::vector<BackTrackerWrapper> >& inputMCPartHitsMap);

    //Fill inputMCPartHits vec with all hits for which the backTracker maximum contributing MCParticle
    //is inputMCPart or daughter particles (if requested)
     //from map of MCParticle->TrackID() to MCParticle
    void getMCPHits(const std::map<int, art::Ptr<simb::MCParticle> >& mcParticleIdMap,
                    const art::Event& evt,
                    const art::Ptr<simb::MCParticle>& inputMCPart,
                    std::vector< art::Ptr<recob::Hit> >& inputMCPartHits,
                    bool bestMatchOnly=true);
     //from art event
    void getMCPHits(const art::Event& evt,
                    const art::Ptr<simb::MCParticle>& inputMCPart,
                    std::vector< art::Ptr<recob::Hit> >& inputMCPartHits,
                    bool bestMatchOnly=true);

  private:

    std::string m_pfp_producer;
    std::string m_cluster_producer;
    std::string m_hit_producer;
    std::string m_pfpToCls_producer;
    std::string m_clsToHit_producer;
    std::string m_mcp_producer;
    std::string m_backTrack_producer;
    //std::string m_hit_process;
    //std::string m_mcp_process;
    //std::string m_backTrack_process;
    bool m_useDaughters;

    void PrintConfig() const;

    void getRecoDaughterChain(const std::map<size_t, art::Ptr<recob::PFParticle> >& pfParticleIdMap,
                              std::vector<size_t>& partChainVec,
                              const std::vector<size_t>& inputDaughters,
                              std::vector<size_t>& newDaughters);
    void getMCDaughterChain(const std::map<int, art::Ptr<simb::MCParticle> >& mcParticleIdMap,
                            std::vector<int>& partChainVec,
                            const std::vector<int>& inputDaughters,
                            std::vector<int>& newDaughters);

    void FillMCHitVector(const std::map<art::Ptr<recob::Hit>,std::vector<BackTrackerWrapper> >& mcPartHitMap,
                         std::vector< art::Ptr<recob::Hit> >& inputMCPartHits, bool bestMatchOnly);

};

} // end namespace mReA

#endif // UBRECO_MREA_FINDPARTICLE_HITS_H
