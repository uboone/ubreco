#ifndef PROTOSHOWER_H
#define PROTOSHOWER_H

#include <TVector3.h>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "Cluster2D.h"

// Data Products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"

namespace protoshower {


class ProtoShower
{
  friend class ProtoShowerHelper;

public:
  ProtoShower() {};
  ~ProtoShower() {};

  TVector3 vertex() const { return _vertex; }
  const std::vector<::cluster2d::Cluster2D> & clusters() const { return _clusters; }

  // getters
  bool hasCluster2D() const {return _hasCluster2D;}
  bool hasCluster3D() const {return _hasCluster3D;}
  bool hasVertex()    const {return _hasVertex;}
  // setters
  void hasCluster2D(bool on) { _hasCluster2D = on; }
  void hasCluster3D(bool on) { _hasCluster3D = on; }
  void hasVertex   (bool on) { _hasVertex    = on; }

  // 3D vertex associated to this protoshower
  TVector3  _vertex;
  
  // list 2D clusters
  std::vector<::cluster2d::Cluster2D> _clusters;

  /**
     Reset function : clear everything and set bools to false
   */
  void Reset() { 
    _clusters.clear(); 
    _hasCluster2D = false;
    _hasCluster3D = false;
    _hasVertex    = false;
    return;
  }

protected:

  bool _hasCluster2D;
  bool _hasCluster3D;
  bool _hasVertex;



  // Not sure what to do with vertexes yet

};

}

#endif
