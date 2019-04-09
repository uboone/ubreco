////////////////////////////////////////////////////////////////////////
// Class:       StopMuMCFilter
// Plugin Type: filter (art v3_01_02)
// File:        StopMuMCFilter_module.cc
//
// Generated at Mon Apr  8 13:33:45 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>

class StopMuMCFilter;


class StopMuMCFilter : public art::EDFilter {
public:
  explicit StopMuMCFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StopMuMCFilter(StopMuMCFilter const&) = delete;
  StopMuMCFilter(StopMuMCFilter&&) = delete;
  StopMuMCFilter& operator=(StopMuMCFilter const&) = delete;
  StopMuMCFilter& operator=(StopMuMCFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::InputTag fMCPproducer;
  float fTrkLen;

};


StopMuMCFilter::StopMuMCFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fMCPproducer = p.get< art::InputTag > ("MCPproducer");
  fTrkLen      = p.get< float         > ("TrkLen");
}

bool StopMuMCFilter::filter(art::Event& e)
{
  
  // load  MCParticle from largeant
  auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);

  // identified number of stopping muons
  int nstopmuon = 0;

  // loop through MCParticles and cut on whether end position is in TPC boundary
  // adn on in-detector track-length

  for (size_t p=0; p < mcp_h->size(); p++) {

    auto const& mcp = mcp_h->at(p);

    if ( fabs(mcp.PdgCode()) != 13) continue;

    std::cout << "Found a muon!" << std::endl;

    double sx, sy, sz, ex, ey, ez;
    sx = 1e4;
    sy = 1e4;
    sz = 1e4;
    bool entered = false;

    // get particle time for drift velocity offset
    auto time = mcp.T(); // in ns
    auto xoffset = 1.1098 * (1e-4) * time; // drift velocity in cm / ns

    ex = mcp.EndX() + xoffset;
    ey = mcp.EndY();
    ez = mcp.EndZ();

    std::cout << "\t end coordinates : [ " << ex << ", " << ey << ", " << ez << " ]" << std::endl;

    if ( (ex < 0)    || (ex > 256. ) ) continue; 
    if ( (ey < -116) || (ey > 116. ) ) continue; 
    if ( (ez < 0)    || (ez > 1036.) ) continue; 

    for (size_t t=0; t < mcp.NumberTrajectoryPoints(); t++) {

      auto x = mcp.Vx(t) + xoffset;
      auto y = mcp.Vy(t);
      auto z = mcp.Vz(t);

      if ( (x < 0)    || (x > 256. ) ) continue; 
      if ( (y < -116) || (y > 116. ) ) continue; 
      if ( (z < 0)    || (z > 1036.) ) continue; 

      // made it this far -> the point is contained
      sx = x;
      sy = y;
      sz = z;

      entered = true;
      break;
    }// for all trajectory points

    if (entered == false) continue;

    double trklen = sqrt( pow((sx-ex),2) + pow((sy-ey),2) + pow((sz-ez),2) );

    std::cout << "\t track length : " << trklen << std::endl;

    if (trklen > fTrkLen) nstopmuon += 1;
    
  }// for all MCParticles

  if (nstopmuon >= 1) {
    std::cout << "Selected! "<< std::endl << std::endl;
    return true;
  }

    std::cout << "Rejected "<< std::endl << std::endl;
  return false;
}

void StopMuMCFilter::beginJob()
{
  // Implementation of optional member function here.
}

void StopMuMCFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(StopMuMCFilter)
