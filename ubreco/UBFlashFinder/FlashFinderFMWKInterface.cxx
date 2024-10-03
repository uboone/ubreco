#ifndef __FLASHFINDERFMWKINTERFACE_CXX__
#define __FLASHFINDERFMWKINTERFACE_CXX__

//#include "cetlib/exception.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "FlashFinderFMWKInterface.h"
namespace pmtana {

  size_t NOpDets() {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->NOpDets();
  }

  size_t OpDetFromOpChannel(size_t opch) {
    auto const& channelMapAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
    return channelMapAlg.OpDetFromOpChannel(opch);
  }

  void OpDetCenterFromOpChannel(size_t opch, double *xyz) {
    auto const& channelMapAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const tmp = channelMapAlg.OpDetGeoFromOpChannel(opch).GetCenter();
    xyz[0] = tmp.X();
    xyz[1] = tmp.Y();
    xyz[2] = tmp.Z();
  }

}
#endif
