////////////////////////////////////////////////////////////////////////
// Class:       ACPTtrigMCFilter
// Plugin Type: filter (art v3_01_02)
// File:        ACPTtrigMCFilter_module.cc
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

class ACPTtrigMCFilter;


class ACPTtrigMCFilter : public art::EDFilter {
public:
  explicit ACPTtrigMCFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ACPTtrigMCFilter(ACPTtrigMCFilter const&) = delete;
  ACPTtrigMCFilter(ACPTtrigMCFilter&&) = delete;
  ACPTtrigMCFilter& operator=(ACPTtrigMCFilter const&) = delete;
  ACPTtrigMCFilter& operator=(ACPTtrigMCFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::InputTag fMCPproducer;
  //float fTrkLen;

  float _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  bool fCheckAnode;
  bool fCheckCathode;

  bool Contained(const float& x, const float& y, const float& z);

};


ACPTtrigMCFilter::ACPTtrigMCFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fMCPproducer = p.get< art::InputTag > ("MCPproducer");
  //  fTrkLen      = p.get< float         > ("TrkLen");

  fCheckAnode   = p.get< bool         > ("CheckAnode");
  fCheckCathode = p.get< bool         > ("CheckCathode");

  _xmin = 0.;
  _xmax = 256.;
  _ymin = -116.;
  _ymax = 116.;
  _zmin = 0.;
  _zmax = 1036.;

}

bool ACPTtrigMCFilter::filter(art::Event& e)
{
  
  // load  MCParticle from largeant
  auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);

  // identified number of stopping muons
  int nACPT = 0;

  // loop through MCParticles and cut on whether end position is in TPC boundary
  // adn on in-detector track-length

  for (size_t p=0; p < mcp_h->size(); p++) {

    auto const& mcp = mcp_h->at(p);

    if ( fabs(mcp.PdgCode()) != 13) continue;

    std::cout << "Found a muon!" << std::endl;

    // does the track cross the cathode?
    //for each trajectory step, check if cathode is crossed

    // get particle time for drift velocity offset
    //auto time = mcp.T(); // in ns
    //auto xoffset = 1.1098 * (1e-4) * time; // drift velocity in cm / ns

    for (size_t t=0; t < mcp.NumberTrajectoryPoints()-1; t++) {

      auto x1 = mcp.Vx(t);
      auto y1 = mcp.Vy(t);
      auto z1 = mcp.Vz(t);

      auto x2 = mcp.Vx(t+1); // + xoffset;
      auto y2 = mcp.Vy(t+1);
      auto z2 = mcp.Vz(t+1);

      // if track has entered into the TPC
      if ( ( Contained(x1,y1,z1) == false ) && ( Contained(x2,y2,z2) == true ) ) {

	// does the point cross the cathode boundary?
	if ( (x1 > _xmax) && (x2 < _xmax) && fCheckCathode )
	  nACPT += 1;
	
	// does the point cross the anode boundary?
	if ( (x1 < _xmin) && (x2 > _xmin) && fCheckAnode)
	  nACPT += 1;
	
      }// if track has entered into the TPC

      // if track has exited the TPC
      if ( ( Contained(x1,y1,z1) == true ) && ( Contained(x2,y2,z2) == false ) ) {

	// does the point cross the x boundary?
	if ( (x1 < _xmax) && (x2 > _xmax) && fCheckCathode)
	  nACPT += 1;

	// does the point cross the anode boundary?
	if ( (x1 > _xmin) && (x2 < _xmin) && fCheckAnode)
	  nACPT += 1;
	
      }// if track has exited the TPC
      
    }// for all trajectory points
    
  }// for all MCParticles
  
  if (nACPT >= 1) {
    std::cout << "Selected! "<< std::endl << std::endl;
    return true;
  }
  
  std::cout << "Rejected "<< std::endl << std::endl;
  return false;
}

bool ACPTtrigMCFilter::Contained(const float& x, const float& y, const float& z) {

  if (x < _xmin) return false;
  if (x > _xmax) return false;
  if (y < _ymin) return false;
  if (y > _ymax) return false;
  if (z < _zmin) return false;
  if (z > _zmax) return false;

  return true;

}

void ACPTtrigMCFilter::beginJob()
{
  // Implementation of optional member function here.
}

void ACPTtrigMCFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ACPTtrigMCFilter)
