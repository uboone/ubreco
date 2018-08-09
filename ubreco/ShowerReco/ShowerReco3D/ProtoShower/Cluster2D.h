#ifndef CLUSTER2D_CLUSTER2D_H
#define CLUSTER2D_CLUSTER2D_H

#include <TVector3.h>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/PxUtils.h"

namespace cluster2d {


class Cluster2D
{

public:
  Cluster2D(){};
  ~Cluster2D(){};

  /**
     @brief Reset containers for cluster
   */
  void Reset() {_hits.clear();};

  /**
     hits, stored as PxHit objects [http://nusoft.fnal.gov/larsoft/doxsvn/html/classutil_1_1PxHit.html]
     PxHit attributes are wire/time, plane, charge, sumADC, peak.
     wire/time, plane, and charge are the ones used here.
     wire/time are stored in centimeters
     time w.r.t. trigger time.
   */
  std::vector<::util::PxHit>  _hits;

  /**
     To Do: Add Plane variable
  */
  unsigned int _plane;

  /**
     2D cluster start point, stored as PxHit
   */
    ::util::PxHit _start;

  /**
     2D cluster end point, stored as PxHit
   */
    ::util::PxHit _end;

  /**
     Opening angle of cluster, in radians
  */
  float _opening_angle;

  /**
     2D angle of cluster (in wire,time coordinate)
     Angle is measured w.r.t. wire-axis rotating counter-clockwise
   */
  float _angle_2d;

  

};

}

#endif
